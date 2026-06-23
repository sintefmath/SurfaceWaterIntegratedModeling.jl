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

    # culvert inlets and outlets, each represented by a culvert ID
    culvert_inlets::Vector{Int}
    culvert_outlets::Vector{Int}

    # tributary paths that merge into this one: (tributary_path_index, junction_cell_index)
    # where junction_cell_index is the 1-based index of the junction cell in *this* path's cells.
    merges::Vector{Tuple{Int,Int}}
end

DynFlowPath(cells, target_trap) = DynFlowPath(cells, target_trap, Int[], Int[], Tuple{Int,Int}[])
DynFlowPath(cells) = DynFlowPath(cells, 0, Int[], Int[], Tuple{Int,Int}[])

"""
        DynTrap(trap_ix, spill_path, culvert_inlets, culvert_outlets)

Represent a trap in the terrain, identified by its index in the spillanalysis
structure.  Every trap must have a spill path (spill_path) that represents the
flow path that water would take when the trap's capacity is exceeded.

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
    # Manning friction recast as a dimensionless loss coefficient (SI form);
    # @@@ SI constant 19.6 (= 2g/Ku^2) -- verify the grouping against HDS-5
    Kf = 19.6 * n^2 * L / D^(4/3)
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
                push!(subnets, _subnetwork(tstruct, cv.inlet,  full_traps))
                push!(subnets, _subnetwork(tstruct, cv.outlet, full_traps))
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
    if ends_with_path
        tix = _unfilled_trap_at(tstruct, paths[end][end], full_traps)
        tix > 0 && push!(traps, tix)
    end

    return _build_network(paths, traps, starts_with_trap, tstruct)
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
function _build_network(paths, traps, starts_with_trap, tstruct)
    CI = CartesianIndices(tstruct.topography)
    link(i, off, n) = i + off <= n ? i + off : 0
    pt = starts_with_trap ? 1 : 0  # offset from a path to the trap it flows into

    dyn_paths = [DynFlowPath([CI[k] for k in paths[i]], link(i, pt, length(traps)))
                 for i in eachindex(paths)]
    dyn_traps = [DynTrap(traps[i], link(i, 1 - pt, length(paths)))
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

    # flatten into one pool with globally unique indices
    all_paths, all_traps, _ = _combine_networks(networks)

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

# Merge all networks into a single flat pool with globally unique indices
function _combine_networks(networks::Vector{DynNetwork})
    path_offsets    = [0; cumsum([length(n.flow_paths) for n in networks])[1:end-1]]
    trap_offsets    = [0; cumsum([length(n.traps)      for n in networks])[1:end-1]]
    culvert_offsets = [0; cumsum([length(n.culverts)   for n in networks])[1:end-1]]

    remap(idx, off)        = idx == 0 ? 0 : idx + off
    remap_vec(v, off)      = isempty(v) ? copy(v) : v .+ off
    remap_merges(v, off)   = [(m + off, j) for (m, j) in v]

    all_paths    = DynFlowPath[]
    all_traps    = DynTrap[]
    all_culverts = DynCulvert[]

    for (ni, net) in enumerate(networks)
        poff, toff, coff = path_offsets[ni], trap_offsets[ni], culvert_offsets[ni]
        for p in net.flow_paths
            push!(all_paths, DynFlowPath(copy(p.cells),
                                         remap(p.target_trap, toff),
                                         remap_vec(p.culvert_inlets, coff),
                                         remap_vec(p.culvert_outlets, coff),
                                         remap_merges(p.merges, poff)))
        end
        for t in net.traps
            push!(all_traps, DynTrap(t.trap_ix,
                                     remap(t.spill_path, poff),
                                     remap_vec(t.culvert_inlets, coff),
                                     remap_vec(t.culvert_outlets, coff)))
        end
        append!(all_culverts, net.culverts)
    end

    return all_paths, all_traps, all_culverts
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
# using a joint topological sort of the path/trap DAG via Graphs.
function _topological_order(global_path_ids, global_trap_ids, all_paths, all_traps)
    np = length(global_path_ids)
    path_node = Dict(gpi => i      for (i, gpi) in enumerate(global_path_ids))
    trap_node = Dict(gti => np + i for (i, gti) in enumerate(global_trap_ids))

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

    # The flow graph is assumed acyclic (water flows strictly downstream); a cycle
    # is a programming error, and topological_sort_by_dfs throws if one is present.
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

    global_path_ids, global_trap_ids =
        _topological_order(global_path_ids, global_trap_ids, all_paths, all_traps)

    path_map = Dict(gpi => lpi for (lpi, gpi) in enumerate(global_path_ids))
    trap_map = Dict(gti => lti for (lti, gti) in enumerate(global_trap_ids))

    # Culverts belonging to this component: both endpoints' owners are necessarily
    # in the same component (the culvert unites them), so testing the inlet owner
    # suffices.  Build the local culvert list and, per owner, its inlet/outlet hits.
    in_comp((kind, id)) = kind == :path ? (id in path_set) :
                          kind == :trap ? (id in keys(trap_map)) : false

    comp_cv      = [ci for ci in eachindex(cv_objs) if in_comp(inlet_owner[ci])]
    culvert_map  = Dict(gci => lci for (lci, gci) in enumerate(comp_cv))
    path_inlets  = Dict{Int,Vector{Int}}();  path_outlets = Dict{Int,Vector{Int}}()
    trap_inlets  = Dict{Int,Vector{Int}}();  trap_outlets = Dict{Int,Vector{Int}}()
    register!(d, owner_id, lc) = push!(get!(d, owner_id, Int[]), lc)
    for gci in comp_cv
        lc = culvert_map[gci]
        ik, iid = inlet_owner[gci]
        ok, oid = outlet_owner[gci]
        register!(ik == :path ? path_inlets  : trap_inlets,  iid, lc)
        register!(ok == :path ? path_outlets : trap_outlets, oid, lc)
    end

    local_paths = [DynFlowPath(
        all_paths[gpi].cells,
        get(trap_map, all_paths[gpi].target_trap, 0),
        get(path_inlets,  gpi, Int[]),
        get(path_outlets, gpi, Int[]),
        [(path_map[m], j) for (m, j) in all_paths[gpi].merges if m ∈ path_set]
    ) for gpi in global_path_ids]

    local_traps = [DynTrap(
        all_traps[gti].trap_ix,
        get(path_map, all_traps[gti].spill_path, 0),
        get(trap_inlets,  gti, Int[]),
        get(trap_outlets, gti, Int[])
    ) for gti in global_trap_ids]

    return DynNetwork(local_paths, local_traps, [cv_objs[gci] for gci in comp_cv])
end
