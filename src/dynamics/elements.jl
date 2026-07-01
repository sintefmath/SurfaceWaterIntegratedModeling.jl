import Graphs

export DynObject, DynFlowPath, DynTrap, DynCulvert, DynNetwork, setup_network

# Make generic baseclass for dynamic objects
abstract type DynObject end

"""
        DynFlowPath(cells, target_trap, culvert_inlets, culvert_outlets, merges)

Represent a flow path over the terrain, represented by a sequence of grid cells.
The flow path may lead into a trap (target_trap), terminate in an intersecting
flow path (target_trap=0), or flow out of the domain (target_trap=0).

The flow path may also include culverts along the way.  Culvert inlets would
subtract water from the flow path, whereas culvert outlets would add water to it.
The infiltration capacity of each cell in the path is represented externally.

"""
struct DynFlowPath <: DynObject
    cells::Vector{CartesianIndex{2}} # cells along the flow path

    # Target trap index (0 for out-of-domain or intersection with another flow path)
    target_trap::Int

    # culvert inlets and outlets on this path: each (culvert_index, cell_position),
    # where cell_position is the 1-based index of the culvert's inlet/outlet cell in
    # *this* path's `cells`.  The position lets routing charge the infiltration up to
    # the abstraction/addition point, exactly like a `merges` junction.
    culvert_inlets::Vector{Tuple{Int,Int}}
    culvert_outlets::Vector{Tuple{Int,Int}}

    # tributary paths that merge into this one: (tributary_path_index, junction_cell_index)
    # where junction_cell_index is the 1-based index of the junction cell in *this* path's cells.
    merges::Vector{Tuple{Int,Int}}
end

DynFlowPath(cells, target_trap) =
    DynFlowPath(cells, target_trap, Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[])
DynFlowPath(cells) =
    DynFlowPath(cells, 0, Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[])

"""
        DynTrap(trap_ix, spill_path, culvert_inlets, culvert_outlets)

Represent a trap in the terrain, identified by its index in the spillanalysis
structure.  The `spill_path` field is the flow path that water would take when the
trap's capacity is exceeded:
- `spill_path > 0`: the (network-local) index of that flow path;
- `spill_path == 0`: the trap is not (yet) full, so it has no spill path;
- `spill_path == -1`: the trap is full and spills straight out of the domain
  (no in-network successor).  The routing and ordering layers treat any
  `spill_path <= 0` as "no successor"; the sentinel only distinguishes a full
  domain-exiting trap from a not-yet-full one for the solver's state validation.

The trap may also have culvert inlets or outlets within its footprint, which
would add or subtract water from the trap respectively, depending on the current
water level in the trap.  The infiltration capacity of each cell in the trap is
represented externally.

"""
struct DynTrap <: DynObject
    trap_ix::Int # index of the trap in the spillanalysis structure

    # spill path
    spill_path::Int

    # culvert inlets and outlets, each represented by a culvert ID
    culvert_inlets::Vector{Int}
    culvert_outlets::Vector{Int}
end

DynTrap(trap_ix, spill_path) = DynTrap(trap_ix, spill_path, Int[], Int[])


"""
       DynCulvert(inlet, outlet, r, Cd, Ke, Kf, Cw)

Represent a culvert connecting two points in the terrain, represented by the
grid cell indices of the inlet and outlet.  The culvert has a limited capacity,
which is determined by its internal parameters as well as the dynamic water levels
at the inlet and outlet.

Internal parameters include the radius (r), discharge coefficient (Cd), entrance
loss coefficient (Ke), friction loss coefficient (Kf), and weir coefficient
(Cw).  Which parameters are used in the computation of the flow through the
culvert depends on the water levels at the inlet and outlet.

"""
struct DynCulvert <: DynObject
    inlet::CartesianIndex{2}  # grid cell index of the culvert inlet
    outlet::CartesianIndex{2} # grid cell index of the culvert outlet

    r::Float64  # radius of the culvert
    Cd::Float64 # discharge coefficient of the culvert
    Ke::Float64 # entrance loss coefficient
    Kf::Float64 # friction loss coefficient
    Cw::Float64 # weir coefficient (for overtopping flow)
end

"""
    DynCulvert(tstruct, inlet, outlet; r, n=0.013, Cd=0.6, Ke=0.5, Cw=1.7)

Convenience constructor that builds a [`DynCulvert`](@ref) from a minimum of
physical data, filling in SI default hydraulic coefficients and deriving the
friction-loss coefficient `Kf` from Manning's roughness `n` and the barrel
length.

The barrel length is the straight-line distance between the `inlet` and `outlet`
cells: horizontal extent plus the elevation drop read from `tstruct.topography`.
All quantities are SI (metres, m^3/s).

Defaults assume a concrete pipe (`n`) with a square-edged entrance (`Cd`, `Ke`)
and the SI weir coefficient (`Cw`); override any of them when better data exists.
"""
function DynCulvert(tstruct, inlet::CartesianIndex{2}, outlet::CartesianIndex{2};
                    r::Real,
                    n::Real  = 0.013,   # Manning roughness (~concrete)
                    Cd::Real = 0.6,     # orifice discharge coef (square edge)
                    Ke::Real = 0.5,     # entrance loss coef (square edge)
                    Cw::Real = 1.7)     # weir coefficient, SI
    D = 2r
    di, dj = Tuple(outlet - inlet)
    # @@@ grid resolution assumed 1 m/cell; use the real cell size once available
    horiz = hypot(float(di), float(dj))
    drop  = tstruct.topography[inlet] - tstruct.topography[outlet]
    L = hypot(horiz, drop)              # full barrel length (m)
    # Manning friction as a dimensionless loss coefficient (HDS-5 eq. 3.5, SI):
    # Kf = Ku * n^2 * L / R^(4/3), with assembled friction constant Ku = 19.63
    # (SI; = 2g/km^2 with Manning unit coef km = 1) and hydraulic radius R = D/4
    # for a full circular pipe.  (R-based form, not D-based -- see §3 of
    # agent/reports/culvert_hydraulics_reference.md.)
    # @@@ R = D/4, not D/2: hydraulic radius is A/P = pi*r^2 / (2*pi*r) = r/2 = D/4.
    R  = D / 4                          # hydraulic radius, full circular pipe
    Kf = 19.63 * n^2 * L / R^(4/3)
    return DynCulvert(inlet, outlet, float(r), float(Cd), float(Ke), float(Kf), float(Cw))
end

# Network of dynamic objects
"""
         DynNetwork(flow_paths, traps, culverts)

Represent the dynamic elements of the terrain as a network of flow paths, traps,
and culverts.  Each flow path may lead into a trap, and each trap has a spill
path that represents the flow path that water would take when the trap's
capacity is exceeded.  Culverts can be present along flow paths or within traps,
and would add or subtract water from the flow path or trap depending on their
location and the current water levels.

The network is tree-like (flow paths can merge but not split, and multiple flow
paths may flow into the same trap, but each flow path can only lead into one
trap and each trap only has one spill path out of it when it overflows).  Loops
are avoided by requiring that water always flows downstreams (whether through
flow paths or culverts).  Overlapping flow paths are not allowed; if multiple
flow paths share cells, one would be truncated as it merges into the other, and
the trap at the end of the truncated path would be connected to the remaining
path instead.
"""
struct DynNetwork
    flow_paths::Vector{DynFlowPath}
    traps::Vector{DynTrap}
    culverts::Vector{DynCulvert}
end

DynNetwork() = DynNetwork(DynFlowPath[], DynTrap[], DynCulvert[])

"""
      setup_network(tstruct, dyn_coords, full_traps; culverts=DynCulvert[])

Given a spillanalysis structure, a set of coordinates representing dynamic points
in the terrain, and a set of traps that are currently filled, build a network
of dynamic flow paths and traps that represent the flow of water through the
terrain.  The flow paths are built starting from the dynamic points, and are
connected to traps as they lead into them.

`culverts` is an optional vector of [`DynCulvert`](@ref)s.  A culvert is included
in the network when either its inlet or its outlet falls on a cell already in the
network (a trap footprint cell or a flow-path cell).  Including a culvert expands
the network to cover the downstream structure of *both* its endpoints (so the
outlet's onward route, and the path the inlet sits on, are present), which may in
turn bring further culverts into the network.  Culverts may connect what would
otherwise be disjoint networks; such networks are merged into one.

Do not pass the (infinite-rate) culverts stored in the `TrapStructure` here —
those are a different concept; `culverts` must be dynamic [`DynCulvert`](@ref)s.
"""
function setup_network(tstruct, dyn_coords, full_traps; culverts::Vector{DynCulvert}=DynCulvert[])
    subnets = DynNetwork[_subnetwork(tstruct, c, full_traps) for c in dyn_coords]
    isempty(culverts) && return _merge_networks(subnets)

    subnets, included = _expand_with_culverts(tstruct, subnets, culverts, full_traps)
    return _merge_networks(subnets, DynCulvert[culverts[ci] for ci in included], tstruct)
end

# All grid cells currently covered by the given subnetworks: every flow-path cell
# and every trap-footprint cell.  Used to decide which culverts touch the network.
function _occupied_cells(tstruct, subnets)
    CI    = CartesianIndices(tstruct.topography)
    cells = Set{CartesianIndex{2}}()
    for net in subnets
        for p in net.flow_paths, c in p.cells
            push!(cells, c)
        end
        for t in net.traps, k in tstruct.footprints[t.trap_ix]
            push!(cells, CI[k])
        end
    end
    return cells
end

# Grow the set of subnetworks by pulling in every culvert that touches the network,
# tracing a fresh downstream subnetwork from *both* endpoints of each included
# culvert (so the outlet's onward route and the inlet's host path/trap are present).
# Iterates to a fixpoint, since a culvert's expansion may expose another culvert.
# Returns the expanded subnetworks and the sorted indices of the included culverts.
function _expand_with_culverts(tstruct, subnets, culverts, full_traps)
    included = Int[]
    incl_set = Set{Int}()
    changed  = true
    while changed
        changed = false
        occ = _occupied_cells(tstruct, subnets)
        for (ci, cv) in enumerate(culverts)
            ci in incl_set && continue
            if cv.inlet in occ || cv.outlet in occ
                push!(included, ci); push!(incl_set, ci); changed = true
                cv.inlet  ∉ occ && push!(subnets, _subnetwork(tstruct, cv.inlet,  full_traps))
                cv.outlet ∉ occ && push!(subnets, _subnetwork(tstruct, cv.outlet, full_traps))
            end
        end
    end
    return subnets, sort(included)
end

"""
        _subnetwork(tstruct, coord, full_traps)

Create a simple line-like network associated with the complete flow path from a
point in the terrain.
"""
function _subnetwork(tstruct, coord::CartesianIndex, full_traps)
    LI = LinearIndices(tstruct.topography)

    paths, ftraps = flow_path_from(tstruct, LI[coord]; full_traps=full_traps)
    @assert abs(length(paths) - length(ftraps)) <= 1

    # The result is a single chain that strictly alternates between path segments
    # and traps in downstream order (a path may be a zero-length connector where one
    # full trap spills directly into the next).  The chain begins with a trap when
    # the start point lies inside the first full trap — i.e. when the start cell is
    # within that trap's footprint.
    starts_with_trap = !isempty(ftraps) && LI[coord] ∈ tstruct.footprints[ftraps[1]]

    # When the chain ends on a path, that path may still discharge into an unfilled
    # trap downstream; include it as the terminating trap.
    traps = collect(ftraps)
    ends_with_path = length(paths) > length(ftraps) ||
                     (starts_with_trap && length(paths) == length(ftraps))

    # Identify the unfilled trap the chain ultimately pools into (0 = exits the domain).
    # Two terminations require subsuming the trailing full descendants into that parent:
    #   * the chain ends on a path discharging into an unfilled trap, or
    #   * the chain ends on a full trap that wrapped back among full siblings of an
    #     unfilled parent (flow_path_from's wrap-around termination).  In that case the
    #     final returning segment was popped, so `ends_with_path` is false and we instead
    #     recover the parent from the hierarchy.  A full trap that spills straight out of
    #     the domain is *not* such a case and is left as the terminal node.
    # A chain that ends on a *full* trap spilling straight out of the domain has that
    # trap as its terminal node (no subsume, no onward path); it is flagged so
    # `_build_network` emits the out-of-domain sentinel spill_path (-1) for it.
    terminal_exits_domain = !ends_with_path && !isempty(ftraps) &&
                            _spills_out_of_domain(tstruct, ftraps[end])

    tix = if ends_with_path
        _unfilled_trap_at(tstruct, paths[end][end], full_traps)
    elseif !isempty(ftraps) && !_spills_out_of_domain(tstruct, ftraps[end])
        _unfilled_parent_of(tstruct, ftraps[end], full_traps)
    else
        0
    end
    if tix > 0
        # Subsume-parents invariant: a parent trap is always a single node that
        # subsumes its whole subtree.  When the terminal unfilled trap `tix` is a
        # parent, its already-full children may have been recorded as separate nodes
        # (they spill into the parent inside its own footprint).  Collapse them into
        # the single parent node, dropping the path segments internal to its
        # footprint.  (A *full* parent is already subsumed by `flow_path_from`, which
        # records the uppermost full supertrap and deletes its footprint cells.)
        paths, traps, starts_with_trap =
            _subsume_terminal_parent(tstruct, paths, ftraps, tix, starts_with_trap)
    end

    return _build_network(paths, traps, starts_with_trap, tstruct;
                          terminal_exits_domain=terminal_exits_domain)
end

# Is `t` a (transitive) subtrap of `P`?  Walk the agglomeration hierarchy upward from
# `t` until `P` is reached (descendant) or the root is passed (not a descendant).
function _is_descendant(tstruct, t::Int, P::Int)
    cur = t
    while true
        par = parentof(tstruct, cur)
        par === nothing && return false
        par == P && return true
        cur = par
    end
end

# Collapse a terminal unfilled parent `P` and its in-chain full descendants into a
# single node.  `ftraps` are the full traps recorded downstream (a prefix of upstream
# non-descendants followed by a contiguous suffix of `P`'s descendants); `paths` are the
# segments between them.  Returns `(paths, traps, starts_with_trap)` for `_build_network`
# with the descendant traps and the path segments internal to `P`'s footprint removed,
# and the last surviving path retargeted into `P`.
function _subsume_terminal_parent(tstruct, paths, ftraps, P::Int, starts_with_trap::Bool)
    d = findfirst(t -> _is_descendant(tstruct, t, P), ftraps)
    if d === nothing
        # `P` has no recorded descendants: keep the original chain, `P` as terminal.
        return paths, push!(collect(ftraps), P), starts_with_trap
    end

    # Keep the upstream non-descendant traps, then the single parent node `P`.
    new_traps = vcat(ftraps[1:d-1], P)

    # Paths that feed into `ftraps[d]` (the first descendant): with a leading trap the
    # path into trap i is path i-1, otherwise it is path i.
    npath_keep = starts_with_trap ? d - 1 : d
    new_paths = [copy(paths[i]) for i in 1:npath_keep]

    # The last kept path enters `P`; drop its cells that lie inside `P`'s footprint
    # (those are internal to the composite and carry no external flow).
    if !isempty(new_paths)
        fpP = Set(tstruct.footprints[P])
        last = new_paths[end]
        cut = findfirst(c -> c ∈ fpP, last)
        cut !== nothing && (new_paths[end] = last[1:cut-1])
    end

    return new_paths, new_traps, starts_with_trap
end

# True when trap `t`'s spillpoint sits on the domain boundary, i.e. `t` spills
# straight out of the domain.  `_add_outer_bounderies!` encodes this as a spillpoint
# whose current and downstream cells coincide (see spillpoints.jl).
_spills_out_of_domain(tstruct, t::Int) =
    (sp = tstruct.spillpoints[t]; sp.current_region_cell == sp.downstream_region_cell)

# Lowest ancestor (supertrap) of `t` that is not full, or 0 if none exists.  This is
# the trap whose basin a wrapped-around chain of full siblings actually pools into.
function _unfilled_parent_of(tstruct, t::Int, full_traps)
    p = parentof(tstruct, t)
    while p !== nothing && p ∈ full_traps
        p = parentof(tstruct, p)
    end
    return p === nothing ? 0 : p
end

# Return the index of the uppermost unfilled trap that `cell` drains into, or 0 if
# the cell drains out of the domain or only into already-full traps.
function _unfilled_trap_at(tstruct, cell::Int, full_traps)
    cur_reg = tstruct.regions[cell]
    cur_reg <= 0 && return 0
    unfilled = filter(x -> x ∉ full_traps, tstruct.supertraps_of[cur_reg])
    isempty(unfilled) ? 0 : minimum(unfilled)
end

# Build the DynNetwork from the alternating path/trap chain.  Both vectors are
# already in downstream order, so their local indices are simply 1:length, and
# `starts_with_trap` fixes the offset between a path and the trap it flows into:
# when the chain starts with a trap, path i flows into trap i+1 and trap i spills
# into path i; otherwise path i flows into trap i and trap i spills into path i+1.
function _build_network(paths, traps, starts_with_trap, tstruct;
                        terminal_exits_domain::Bool=false)
    CI = CartesianIndices(tstruct.topography)
    link(i, off, n) = i + off <= n ? i + off : 0
    pt = starts_with_trap ? 1 : 0  # offset from a path to the trap it flows into

    # spill_path for trap i.  A terminal trap (no in-network downstream path) gets 0
    # by default, which the solver reads as TRANSITORY (not yet full).  When that
    # terminal trap is instead a *full* trap spilling straight out of the domain, it
    # is flagged with the sentinel -1 so the solver can tell the two apart (a FULL
    # trap must spill somewhere; -1 says "out of the domain").  Routing and ordering
    # already treat any spill_path <= 0 as "no in-network successor", so only the
    # FULL/TRANSITORY validation distinguishes -1 from 0.
    function spill_path_of(i)
        sp = link(i, 1 - pt, length(paths))
        (sp == 0 && terminal_exits_domain && i == length(traps)) ? -1 : sp
    end

    dyn_paths = [DynFlowPath([CI[k] for k in paths[i]], link(i, pt, length(traps)))
                 for i in eachindex(paths)]
    dyn_traps = [DynTrap(traps[i], spill_path_of(i))
                 for i in eachindex(traps)]

    return DynNetwork(dyn_paths, dyn_traps, DynCulvert[])
end

"""
        _merge_networks(networks)

Given a vector of flow path networks, merge any that overlap into a single one.
Return a vector of disjoint networks.  Where flow paths share grid cells, the
later path is truncated at the shared cell and registered as a tributary (a
"merge") of the path it runs into, so that no grid cell is shared by multiple
flow paths.

"""
_merge_networks(networks::Vector{DynNetwork}) =
    _merge_networks(networks, DynCulvert[], nothing)

function _merge_networks(networks::Vector{DynNetwork}, cv_objs::Vector{DynCulvert}, tstruct)
    isempty(networks) && return networks

    # flatten the subnets into one pool with globally unique indices
    all_paths, all_traps = _combine_subnets(networks)

    # truncate paths that share cells, registering each truncated path as a tributary
    _resolve_cell_overlaps!(all_paths)

    # collapse duplicate trap entries (the same physical trap reached from several
    # subnetworks), so each trap_ix appears once with merged connectivity
    all_paths, all_traps = _dedup_traps(all_paths, all_traps)

    # resolve each culvert's inlet/outlet to its owning path or trap
    inlet_owner, outlet_owner = _culvert_owners(tstruct, all_paths, all_traps, cv_objs)

    # _components: group paths+traps into disjoint connected components (culverts
    # link the components of their two endpoints); _build_component: rebuild each.
    return [_build_component(all_paths, all_traps, cv_objs, inlet_owner, outlet_owner, pids, tids)
            for (pids, tids) in _components(all_paths, all_traps, inlet_owner, outlet_owner)]
end

# Collapse duplicate trap entries (same `trap_ix` reached from several subnetworks)
# into one.  The surviving entry keeps a non-zero spill_path if any duplicate had
# one.  Path `target_trap` references are remapped to the surviving index;
# path/spill_path (path references) are untouched.  Culverts are assigned later, by
# cell (_culvert_owners / _build_component), so traps carry none at this stage.
function _dedup_traps(all_paths, all_traps)
    canon     = Dict{Int,Int}()              # trap_ix -> surviving index in new_traps
    new_traps = DynTrap[]
    remap     = zeros(Int, length(all_traps)) # old trap index -> new trap index

    for (ti, t) in enumerate(all_traps)
        if haskey(canon, t.trap_ix)
            ni  = canon[t.trap_ix]
            old = new_traps[ni]
            new_traps[ni] = DynTrap(t.trap_ix,
                old.spill_path == 0 ? t.spill_path : old.spill_path)
            remap[ti] = ni
        else
            push!(new_traps, t)
            canon[t.trap_ix] = length(new_traps)
            remap[ti]        = length(new_traps)
        end
    end

    new_paths = [DynFlowPath(p.cells,
                             p.target_trap == 0 ? 0 : remap[p.target_trap],
                             p.culvert_inlets, p.culvert_outlets, p.merges)
                 for p in all_paths]
    return new_paths, new_traps
end

# For each culvert (in `cv_objs` order), resolve its inlet and outlet cell to the
# path or trap that owns it, as a `(:path|:trap, global_index)` pair.  A trap
# footprint cell wins over a flow-path cell (footprints are the trap interior).
# Returns two parallel vectors (inlet owners, outlet owners).  Empty when there
# are no culverts (and then `tstruct` may be `nothing`).
function _culvert_owners(tstruct, all_paths, all_traps, cv_objs)
    isempty(cv_objs) && return (Tuple{Symbol,Int}[], Tuple{Symbol,Int}[])

    CI        = CartesianIndices(tstruct.topography)
    trap_cell = Dict{CartesianIndex{2},Int}()
    for (ti, t) in enumerate(all_traps), k in tstruct.footprints[t.trap_ix]
        trap_cell[CI[k]] = ti
    end
    path_cell = Dict{CartesianIndex{2},Int}()
    for (pi, p) in enumerate(all_paths), c in p.cells
        path_cell[c] = pi
    end
    owner(cell) = haskey(trap_cell, cell) ? (:trap, trap_cell[cell]) :
                  haskey(path_cell, cell) ? (:path, path_cell[cell]) :
                  (:none, 0)

    inlet_owner  = [owner(cv.inlet)  for cv in cv_objs]
    outlet_owner = [owner(cv.outlet) for cv in cv_objs]
    return inlet_owner, outlet_owner
end

# Flatten a set of subnets into a single pool of paths and traps with globally
# unique indices.  Subnets are culvert-free (see `_subnetwork`/`_build_network`),
# so culverts are not handled here: a culvert's owning path/trap is resolved
# later, by cell, in `_culvert_owners` / `_build_component`.
function _combine_subnets(subnets::Vector{DynNetwork})
    path_offsets = [0; cumsum([length(n.flow_paths) for n in subnets])[1:end-1]]
    trap_offsets = [0; cumsum([length(n.traps)      for n in subnets])[1:end-1]]

    # Non-positive references are sentinels, not indices: 0 = none, -1 = spills out
    # of the domain (trap spill_path only); both pass through an offset unchanged.
    remap(idx, off)      = idx <= 0 ? idx : idx + off
    remap_merges(v, off) = [(m + off, j) for (m, j) in v]

    all_paths = DynFlowPath[]
    all_traps = DynTrap[]

    for (ni, net) in enumerate(subnets)
        poff, toff = path_offsets[ni], trap_offsets[ni]
        for p in net.flow_paths
            push!(all_paths, DynFlowPath(copy(p.cells),
                                         remap(p.target_trap, toff),
                                         Tuple{Int,Int}[], Tuple{Int,Int}[],   # culverts added later
                                         remap_merges(p.merges, poff)))
        end
        for t in net.traps
            push!(all_traps, DynTrap(t.trap_ix,
                                     remap(t.spill_path, poff),
                                     Int[], Int[]))             # culverts added later
        end
    end

    return all_paths, all_traps
end

# Truncate any flow path whose cells overlap with a previously-processed path.
# The truncated path is registered as a tributary (a "merge") of the primary
# (earlier) path that owns the shared cell.  Trap connectivity is left untouched.
# Culverts are assigned later, by cell (_culvert_owners / _build_component), so
# paths carry none at this stage.
function _resolve_cell_overlaps!(all_paths)
    cell_owner = Dict{CartesianIndex{2}, Int}()  # grid cell → owning path index

    for pi in 1:length(all_paths)
        path = all_paths[pi]
        merge_pos = findfirst(cell -> haskey(cell_owner, cell), path.cells)

        if merge_pos !== nothing
            merge_into = cell_owner[path.cells[merge_pos]]
            kept       = path.cells[1:merge_pos-1]

            all_paths[pi] = DynFlowPath(kept, 0, path.culvert_inlets, path.culvert_outlets, path.merges)

            primary      = all_paths[merge_into]
            junction_pos = findfirst(==(path.cells[merge_pos]), primary.cells)
            all_paths[merge_into] = DynFlowPath(primary.cells, primary.target_trap,
                primary.culvert_inlets, primary.culvert_outlets,
                [primary.merges; (pi, junction_pos)])

            for cell in kept
                cell_owner[cell] = pi
            end
        else
            for cell in path.cells
                cell_owner[cell] = pi
            end
        end
    end
end

# Group paths and traps into disjoint connected components using union-find over a
# unified node set: paths are nodes 1:np, traps are nodes np+1:np+nt.  Nodes are
# connected when a path targets a trap, a trap spills into a path, a tributary
# merges into a path, or a culvert links its inlet owner to its outlet owner.
# Returns a vector of `(path_ids, trap_ids)` tuples (global indices per component).
# A trap with no connections forms its own singleton component (a lone terminal
# trap, e.g. reached only via a culvert), which is preserved rather than dropped.
function _components(all_paths, all_traps, inlet_owner, outlet_owner)
    np, nt = length(all_paths), length(all_traps)
    parent = collect(1:(np + nt))

    find_root(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    function unite!(x, y)
        rx, ry = find_root(x), find_root(y)
        rx != ry && (parent[rx] = ry)
    end
    node((kind, id)) = kind == :path ? id : np + id

    for (pi, p) in enumerate(all_paths)
        p.target_trap > 0 && unite!(pi, np + p.target_trap)
        for (m, _) in p.merges
            unite!(m, pi)
        end
    end
    for (ti, t) in enumerate(all_traps)
        t.spill_path > 0 && unite!(np + ti, t.spill_path)
    end
    for k in eachindex(inlet_owner)
        io, oo = inlet_owner[k], outlet_owner[k]
        (io[1] == :none || oo[1] == :none) && continue
        unite!(node(io), node(oo))
    end

    paths_of = Dict{Int, Vector{Int}}()
    traps_of = Dict{Int, Vector{Int}}()
    for pi in 1:np
        push!(get!(paths_of, find_root(pi), Int[]), pi)
    end
    for ti in 1:nt
        push!(get!(traps_of, find_root(np + ti), Int[]), ti)
    end

    roots = union(keys(paths_of), keys(traps_of))
    return [(get(paths_of, r, Int[]), get(traps_of, r, Int[])) for r in roots]
end

# Return (sorted_path_ids, sorted_trap_ids) in upstream-to-downstream order,
# using a joint topological sort of the path/trap DAG via Graphs.  `culvert_links`
# is a list of (inlet_owner, outlet_owner) `(:path|:trap, global_id)` pairs whose
# direction must be respected (inlet processed before outlet) on top of the
# terrain flow edges.
function _topological_order(global_path_ids, global_trap_ids, all_paths, all_traps,
                            culvert_links = Tuple{Tuple{Symbol,Int},Tuple{Symbol,Int}}[])
    np = length(global_path_ids)
    path_node = Dict(gpi => i      for (i, gpi) in enumerate(global_path_ids))
    trap_node = Dict(gti => np + i for (i, gti) in enumerate(global_trap_ids))
    node_of((kind, gid)) = kind == :path ? path_node[gid] : trap_node[gid]

    # Every path/trap referenced from within this component is itself part of the
    # component: _components unions paths and traps via merges, path→trap targets,
    # trap→spill_path, and culvert endpoint links.  So a nonzero reference is
    # guaranteed to resolve here (asserted, not tested).
    g = Graphs.SimpleDiGraph(np + length(global_trap_ids))
    for (li, gpi) in enumerate(global_path_ids)
        p = all_paths[gpi]
        if p.target_trap > 0
            @assert haskey(trap_node, p.target_trap)
            Graphs.add_edge!(g, li, trap_node[p.target_trap])
        end
        for (m, _) in p.merges
            @assert haskey(path_node, m)
            Graphs.add_edge!(g, path_node[m], li)
        end
    end
    for (li, gti) in enumerate(global_trap_ids)
        sp = all_traps[gti].spill_path
        if sp > 0
            @assert haskey(path_node, sp)
            Graphs.add_edge!(g, np + li, path_node[sp])
        end
    end
    # Culvert direction: the inlet owner is processed before the outlet owner
    # (downhill assumption), so add an edge inlet_owner -> outlet_owner.  This is
    # what makes a culvert that links two otherwise terrain-disjoint pieces order
    # them correctly.  A culvert running against terrain flow would close a cycle
    # here -- that is the deferred uphill / reverse-culvert case, and
    # topological_sort_by_dfs throws on it (fail loud rather than mis-route).
    for (inlet_owner, outlet_owner) in culvert_links
        (inlet_owner[1] == :none || outlet_owner[1] == :none) && continue
        Graphs.add_edge!(g, node_of(inlet_owner), node_of(outlet_owner))
    end

    # The network must be acyclic to be ordered upstream-to-downstream.  A cycle is
    # either a terrain-flow programming error or, more commonly, an uphill/reverse
    # culvert (inlet downstream of its outlet).  Fail loud rather than mis-route.
    # @@@ uphill / reverse culverts are deferred; revisit when task-2 routing gains
    #     direction-aware handling instead of rejecting the network here.
    if Graphs.is_cyclic(g)
        error("Cyclic dynamic network: flow paths, traps, and culverts form a " *
              "directed loop and cannot be ordered upstream-to-downstream.  The " *
              "usual cause is an uphill/reverse culvert (its inlet lies downstream " *
              "of its outlet), which is not yet supported.")
    end
    order = Graphs.topological_sort_by_dfs(g)
    return ([global_path_ids[i]      for i in order if i <= np],
            [global_trap_ids[i - np] for i in order if i >  np])
end

# Reconstruct one DynNetwork from this component's global path and trap indices,
# remapping all internal references (including culverts) to local 1-based indices.
# Culverts whose inlet/outlet owners lie in this component are attached to the
# owning local path or trap via its culvert_inlets / culvert_outlets list.
function _build_component(all_paths, all_traps, cv_objs, inlet_owner, outlet_owner,
                          global_path_ids, global_trap_ids)
    path_set = Set(global_path_ids)
    trap_set = Set(global_trap_ids)

    # Culverts belonging to this component: both endpoints' owners are necessarily
    # in the same component (the culvert unites them), so testing the inlet owner
    # suffices.  (Sets are order-independent, so this is valid before the sort.)
    in_comp((kind, id)) = kind == :path ? (id in path_set) :
                          kind == :trap ? (id in trap_set) : false
    comp_cv = [ci for ci in eachindex(cv_objs) if in_comp(inlet_owner[ci])]

    # Order paths/traps upstream-to-downstream, honoring culvert direction (inlet
    # owner before outlet owner) on top of terrain flow.
    culvert_links = [(inlet_owner[ci], outlet_owner[ci]) for ci in comp_cv]
    global_path_ids, global_trap_ids =
        _topological_order(global_path_ids, global_trap_ids, all_paths, all_traps,
                           culvert_links)

    path_map = Dict(gpi => lpi for (lpi, gpi) in enumerate(global_path_ids))
    trap_map = Dict(gti => lti for (lti, gti) in enumerate(global_trap_ids))

    # Per owner, the local culvert indices it hosts.  A path also records the
    # 1-based position of the culvert's inlet/outlet cell within its `cells`, so
    # routing can charge infiltration up to that point (like a `merges` junction);
    # a trap has no along-path position, so it stores the bare culvert index.
    culvert_map  = Dict(gci => lci for (lci, gci) in enumerate(comp_cv))
    path_inlets  = Dict{Int,Vector{Tuple{Int,Int}}}();  path_outlets = Dict{Int,Vector{Tuple{Int,Int}}}()
    trap_inlets  = Dict{Int,Vector{Int}}();             trap_outlets = Dict{Int,Vector{Int}}()
    pos_on(gpi, cell) = findfirst(==(cell), all_paths[gpi].cells)
    for gci in comp_cv
        lc = culvert_map[gci]
        cv = cv_objs[gci]
        ik, iid = inlet_owner[gci]
        ok, oid = outlet_owner[gci]
        ik == :path ? push!(get!(path_inlets,  iid, Tuple{Int,Int}[]), (lc, pos_on(iid, cv.inlet))) :
                      push!(get!(trap_inlets,  iid, Int[]), lc)
        ok == :path ? push!(get!(path_outlets, oid, Tuple{Int,Int}[]), (lc, pos_on(oid, cv.outlet))) :
                      push!(get!(trap_outlets, oid, Int[]), lc)
    end

    local_paths = [DynFlowPath(
        all_paths[gpi].cells,
        get(trap_map, all_paths[gpi].target_trap, 0),
        get(path_inlets,  gpi, Tuple{Int,Int}[]),
        get(path_outlets, gpi, Tuple{Int,Int}[]),
        [(path_map[m], j) for (m, j) in all_paths[gpi].merges if m ∈ path_set]
    ) for gpi in global_path_ids]

    # spill_path < 0 is the out-of-domain sentinel (-1): keep it as-is rather than
    # remapping it as a path index.
    remap_spill(sp) = sp < 0 ? sp : get(path_map, sp, 0)
    local_traps = [DynTrap(
        all_traps[gti].trap_ix,
        remap_spill(all_traps[gti].spill_path),
        get(trap_inlets,  gti, Int[]),
        get(trap_outlets, gti, Int[])
    ) for gti in global_trap_ids]

    return DynNetwork(local_paths, local_traps, [cv_objs[gci] for gci in comp_cv])
end
