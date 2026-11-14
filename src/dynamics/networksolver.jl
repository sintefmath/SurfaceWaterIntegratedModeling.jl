import Interpolations
import Graphs
using DifferentialEquations   # brings DiscreteCallback, CallbackSet, SciMLBase into scope

export dynNetworkRateFunction!, solveDynNetwork!

# ============================================================================
# Geometry helpers for the dynamic network solver.
#
# A `DynNetwork` carries only topology; the geometry lives in the `TrapStructure`.
# Precomputed once per solve, these answer the two questions the rate function asks
# of each trap at a stored volume V:
#   * water_level(g, V)         -- water surface elevation
#   * wetted_infiltration(g, V) -- infiltration over the currently-submerged footprint
# Reuses `_compute_z_vol_tables` (volume<->level) and the `infilfun(bottom .<= z)`
# wetted-area pattern from fill_sequence.jl.
# ============================================================================

"""
    TrapGeometry

Precomputed geometry of one trap in a [`DynNetwork`](@ref), used by the rate function.
Volumes are *net of subtraps*: `capacity` is the trap's own storage above its subtraps,
and the state `V` ranges over `0 .. capacity`.

# Fields
- `trap_ix`: index of the trap in the [`TrapStructure`](@ref)
- `capacity`: own storage volume, net of subtraps
- `footprint`: linear indices of the footprint cells
- `bottom`: terrain bottom per footprint cell; parent-trap cells over a subtrap are raised
  to the subtrap's spillpoint so a submerged subtrap reads as wetted
- `infil`: infiltration rate per footprint cell (0 where impermeable)
- `zmin`: lowest bottom elevation (level of an empty trap)
- `v2z`: stored-volume -> water-level interpolation (see [`water_level`](@ref))
"""
struct TrapGeometry
    trap_ix::Int            # index of the trap in the TrapStructure
    capacity::Float64       # own storage volume, net of subtraps
    footprint::Vector{Int}  # linear indices of the footprint cells
    bottom::Vector{Float64} # terrain bottom per footprint cell
    infil::Vector{Float64}  # infiltration rate per footprint cell (0 where impermeable)
    zmin::Float64           # lowest bottom elevation (level of an empty trap)
    v2z   # volume -> level map; untyped (perf is a later concern, correctness first)
end

# ----------------------------------------------------------------------------
"""
    _build_trap_geometry(tstruct, net, infiltration; zvt=nothing) -> Vector{TrapGeometry}

Precompute the geometry the rate function needs for every trap in `net`.

# Arguments
- `tstruct`: supplies each trap's footprint, topography, volumes and spillpoint.
- `net`: the network; its `traps` fix the output order.
- `infiltration`: per-cell rate grid (0 = impermeable).
- `zvt`: cached volume↔level tables ([`_compute_z_vol_tables`](@ref)); computed here if
  `nothing`.

# Returns
One [`TrapGeometry`](@ref) per trap, in `net.traps` order.  A parent's `bottom` is raised to
its subtraps' spillpoint and its `capacity` is net of them; footprint cells at or above the
spillpoint get `infil = 0`, since they never pond.  Nothing is mutated.
"""
function _build_trap_geometry(tstruct::TrapStructure,
                              net::DynNetwork,
                              infiltration::AbstractMatrix{<:Real};
                              zvt = nothing)

    zvt === nothing && (zvt = _compute_z_vol_tables(tstruct))
    nreg = numregions(tstruct)

    geom = Vector{TrapGeometry}(undef, length(net.traps))
    for (i, dtrap) in enumerate(net.traps)
        tix       = dtrap.trap_ix
        footprint = tstruct.footprints[tix]
        bottom    = Float64.(tstruct.topography[footprint])

        # parent trap: raise the bottom to the upper subtrap's spillpoint, capacity net of
        # subtraps.  leaf trap: no subtraps, capacity is the full trap volume.
        if tix > nreg
            children = subtrapsof(tstruct, tix)
            bottom   = max.(bottom, tstruct.spillpoints[children[1]].elevation)
            capacity = Float64(tstruct.trapvolumes[tix] - tstruct.subvolumes[tix])
        else
            capacity = Float64(tstruct.trapvolumes[tix])
        end

        zmin    = minimum(bottom)
        zvtable = zvt[tix]
        # a single-row table is a degenerate flat trap with no storage: level is always zmin
        v2z = length(zvtable[2]) == 1 ?
            (_ -> zmin) :
            Interpolations.linear_interpolation(zvtable[2], zvtable[1],
                                                extrapolation_bc = Interpolations.Line())

        infil = Float64.(infiltration[footprint])

        # Cells at/above the spillpoint never pond, so drop their infiltration.  This keeps the
        # wetted-area loss continuous at V = capacity — else a trap pinned at its spillpoint by a
        # balancing through-flow chatters full<->not-full.  Test the actual terrain, not raised
        # `bottom` (same rule as `_ponding_infiltration` / `fill_trap_until`).
        _sp = Float64(tstruct.spillpoints[tix].elevation)   # concrete: Spillpoint.elevation is ::Real
        @inbounds for k in eachindex(infil)
            tstruct.topography[footprint[k]] >= _sp && (infil[k] = 0.0)
        end

        geom[i] = TrapGeometry(tix, capacity, footprint, bottom, infil, zmin, v2z)
    end
    return geom
end

# ----------------------------------------------------------------------------
"""
    water_level(g::TrapGeometry, V) -> Float64

Water surface elevation of a trap at a given stored volume.

# Arguments
- `g`: the trap's geometry.
- `V`: stored volume, net of subtraps.

# Returns
The surface elevation: `g.zmin` when empty (`V <= 0`), else `g.v2z(V)`, held at the spillpoint
once `V` reaches `g.capacity`.
"""
function water_level(g::TrapGeometry, V::Real)
    return V <= 0.0 ? g.zmin : Float64(g.v2z(min(V, g.capacity)))
end

# ----------------------------------------------------------------------------
"""
    wetted_infiltration(g::TrapGeometry, V) -> Float64

Infiltration loss of a trap at a given volume — only the submerged part of its footprint loses
water.

# Arguments
- `g`: the trap's geometry.
- `V`: stored volume, net of subtraps.

# Returns
The per-cell rate summed over cells whose `bottom` lies at or below [`water_level`](@ref)`(g, V)`.
A step function of `V` (cells wet in whole units), which is why the steady-state detector must be
discrete.  At capacity the level reaches the spillpoint and the whole footprint infiltrates,
matching `_footprint_infiltration`.
"""
function wetted_infiltration(g::TrapGeometry, V::Real)
    z = water_level(g, V)
    s = 0.0
    @inbounds for k in eachindex(g.bottom)
        if g.bottom[k] <= z
            s += g.infil[k]
        end
    end
    return s
end

# ============================================================================
# Flow routing through the network (lossy paths).
#
# Flow paths are stateless: they route water downstream instantaneously but lossily —
# a path cell with infiltration removes water in transit, capped by the flow passing.
# A full trap emits `inflow - footprint infiltration` down its spill path; an
# accumulating trap keeps its inflow.  Tributaries (`DynFlowPath.merges`) deliver into
# the path they join.  The router's output is the total inflow arriving at each trap.
# ============================================================================

"""
    _network_order(net) -> (order, np)

Visiting order for the router: each element comes after everything upstream of it.

# Arguments
- `net`: the network to order.

# Returns
`(order, np)` — `order` a topological sort over the combined nodes, `np = length(net.flow_paths)`
the offset that splits them: node `p <= np` is flow path `p`, node `q > np` is trap `q - np`.

Edges follow the flow: trap → its spill path, path → its target trap, tributary → the path it
joins, culvert inlet owner → outlet owner.  `net.traps` order alone won't do, since it doesn't
fix the path/trap interleaving once tributaries exist.  Nothing is mutated.

!!! note "Why culverts get an edge and NBS placements do not"
    An edge is needed only for a read-after-write *within* one routing pass.  A culvert has
    one: its inlet draws `min(capacity, current)` — what is actually available there — and its
    outlet then delivers exactly that, so the inlet must be visited first.

    An NBS has none.  Its output depends only on the stored layer state, so
    [`_nbs_injections`](@ref) fixes it before the pass, and the draws the pass collects feed
    only the layers' rates.  A carrier path may therefore be visited before the inflow path
    feeding the same placement.

    NBS placements *do* get edges in `_coupling_graph` (`network_utils.jl`), but that is an
    undirected graph deciding which elements must be solved **together**, not in what order.
"""
function _network_order(net::DynNetwork)
    np = length(net.flow_paths)
    nt = length(net.traps)
    g  = Graphs.SimpleDiGraph(np + nt)

    for (p, path) in enumerate(net.flow_paths)
        path.target_trap > 0 && Graphs.add_edge!(g, p, np + path.target_trap)
        for (m, _) in path.merges
            Graphs.add_edge!(g, m, p)            # tributary m feeds path p
        end
    end
    for (t, trap) in enumerate(net.traps)
        trap.spill_path > 0 && Graphs.add_edge!(g, np + t, trap.spill_path)
    end
    _add_culvert_edges!(g, net, np)        # inlet owner before outlet owner

    return Graphs.topological_sort_by_dfs(g), np
end

# ----------------------------------------------------------------------------
"""
    _merge_targets(net) -> Vector{Int}

Which path each tributary feeds into.

# Arguments
- `net`: the network to scan.

# Returns
A vector over `net.flow_paths`: `merge_target[p]` is the path tributary `p` joins, or `0` if `p`
is not a tributary.  Well defined because a path is a tributary of at most one other.  Nothing
is mutated.
"""
function _merge_targets(net::DynNetwork)
    merge_target = zeros(Int, length(net.flow_paths))
    for (a, path) in enumerate(net.flow_paths), (m, _) in path.merges
        merge_target[m] = a
    end
    return merge_target
end

# ----------------------------------------------------------------------------
"""
    CulvertPlan

Static culvert routing data for one solve — everything [`_culvert_flow`](@ref) needs to price a
culvert without re-deriving it each step.  Every field is a vector indexed by culvert.

# Fields
- `diam`: barrel diameter (`2r`).
- `inlet_invert`, `outlet_invert`: terrain elevation at each end's cell.
- `inlet_is_path`, `outlet_is_path`: whether that end is owned by a flow path (else a trap).
- `inlet_owner`, `outlet_owner`: the owning path's or trap's local index.
"""
struct CulvertPlan
    diam::Vector{Float64}         # barrel diameter (2r) per culvert
    inlet_invert::Vector{Float64} # terrain elevation at each culvert's inlet cell
    outlet_invert::Vector{Float64}# terrain elevation at each culvert's outlet cell
    inlet_is_path::BitVector
    inlet_owner::Vector{Int}      # local path or trap index owning the inlet
    outlet_is_path::BitVector
    outlet_owner::Vector{Int}     # local path or trap index owning the outlet
end

"""
    _build_culvert_plan(net, tstruct) -> CulvertPlan

Resolve each culvert's ends to their owners, once per solve.

# Arguments
- `net`: supplies the culverts and the paths/traps that host their endpoints.
- `tstruct`: supplies the terrain the inverts are read from.

# Returns
The [`CulvertPlan`](@ref) for `net.culverts`.  Nothing is mutated.
"""
function _build_culvert_plan(net::DynNetwork, tstruct)
    nc   = length(net.culverts)
    plan = CulvertPlan([2 * cv.r for cv in net.culverts],
                       Float64[tstruct.topography[cv.inlet]  for cv in net.culverts],
                       Float64[tstruct.topography[cv.outlet] for cv in net.culverts],
                       falses(nc), zeros(Int, nc), falses(nc), zeros(Int, nc))
    for (p, path) in enumerate(net.flow_paths)
        for (ci, _) in path.culvert_inlets;  plan.inlet_is_path[ci]  = true; plan.inlet_owner[ci]  = p; end
        for (ci, _) in path.culvert_outlets; plan.outlet_is_path[ci] = true; plan.outlet_owner[ci] = p; end
    end
    for (t, trap) in enumerate(net.traps)
        for ci in trap.culvert_inlets;  plan.inlet_owner[ci]  = t; end   # is_path stays false
        for ci in trap.culvert_outlets; plan.outlet_owner[ci] = t; end
    end
    return plan
end

"""
    _add_culvert_edges!(g, net, np) -> g

Order each culvert's inlet before its outlet, so the drawn flow is known when it is delivered.

# Arguments
- `g`: the routing-order graph, **mutated** — one inlet-owner → outlet-owner edge per culvert.
- `net`: supplies the culverts and their endpoint owners.
- `np`: the path/trap node offset (see [`_network_order`](@ref)).

# Returns
`g`, mutated.

Assumes downhill culverts: the network is verified acyclic with these edges at construction, so
one running against the flow order would fail to sort.
"""
function _add_culvert_edges!(g, net::DynNetwork, np::Int)
    inlet_node  = zeros(Int, length(net.culverts))
    outlet_node = zeros(Int, length(net.culverts))
    for (p, path) in enumerate(net.flow_paths)
        for (ci, _) in path.culvert_inlets;  inlet_node[ci]  = p; end
        for (ci, _) in path.culvert_outlets; outlet_node[ci] = p; end
    end
    for (t, trap) in enumerate(net.traps)
        for ci in trap.culvert_inlets;  inlet_node[ci]  = np + t; end
        for ci in trap.culvert_outlets; outlet_node[ci] = np + t; end
    end
    for ci in eachindex(net.culverts)
        (inlet_node[ci] > 0 && outlet_node[ci] > 0) &&
            Graphs.add_edge!(g, inlet_node[ci], outlet_node[ci])
    end
    return g
end

"""
    _path_event_templates(net) -> Vector{Vector{Tuple{Int,Symbol,Int}}}

The ordered stops along each flow path — everything the router must act on mid-path.

# Arguments
- `net`: supplies each path's merges, culvert endpoints and NBS inlets.

# Returns
One entry per flow path: a vector of `(position, kind, idx)` sorted by `position` (the 1-based
cell index within that path), `kind` one of `:trib`, `:cvin`, `:cvout`, `:nbsin`, and `idx` the
tributary path / culvert / NBS it refers to.

Walking these in order is what lets [`_path_delivered!`](@ref) charge each segment only the
cells it actually crosses.  Static for a solve, so built once and reused every step.  Nothing
is mutated.
"""
function _path_event_templates(net::DynNetwork)
    return [begin
                ev = Tuple{Int,Symbol,Int}[]
                for (m, j)   in fp.merges;          push!(ev, (j,   :trib,  m));  end
                for (ci, pos) in fp.culvert_inlets;  push!(ev, (pos, :cvin,  ci)); end
                for (ci, pos) in fp.culvert_outlets; push!(ev, (pos, :cvout, ci)); end
                for (ni, pos) in fp.nbs_inlets;      push!(ev, (pos, :nbsin, ni)); end
                sort!(ev; by = first)
                ev
            end
            for fp in net.flow_paths]
end

"""
    _culvert_flow(plan, net, ci, trap_level) -> Float64

Volumetric capacity of one culvert at the current water levels — what it *could* carry, before
the inlet's available flow caps it (see [`_path_delivered!`](@ref)).

# Arguments
- `plan`: the [`CulvertPlan`](@ref) giving this culvert's owners, inverts and diameter.
- `net`: supplies the `DynCulvert` itself.
- `ci`: culvert index.
- `trap_level`: current water-surface elevation per trap, `net.traps` order.

# Returns
Capacity in m^3/s, always `>= 0`.

A trap endpoint uses its live level above the culvert's invert.  A flow-path endpoint has no
pool, so it is treated as not-submerged with head fixed at the diameter — its inlet-control
capacity is then the weir capacity at head `D`.  Downhill-only (`allow_reverse = false`): a
higher outlet pool yields 0, not reverse flow.  Nothing is mutated.
"""
function _culvert_flow(plan::CulvertPlan, net::DynNetwork, ci::Int,
                       trap_level::AbstractVector{<:Real})
    cv = net.culverts[ci]
    D  = plan.diam[ci]
    if plan.inlet_is_path[ci]
        ih, isub = D, false
    else
        ih   = max(trap_level[plan.inlet_owner[ci]] - plan.inlet_invert[ci], 0.0)
        isub = ih >= D
    end
    if plan.outlet_is_path[ci]
        oh, osub = D, false
    else
        oh   = max(trap_level[plan.outlet_owner[ci]] - plan.outlet_invert[ci], 0.0)
        osub = oh >= D
    end
    return culvert_rate(cv;
                        inlet_submerged = isub,  inlet_head  = ih,  inlet_invert  = plan.inlet_invert[ci],
                        outlet_submerged = osub, outlet_head = oh,  outlet_invert = plan.outlet_invert[ci],
                        allow_reverse = false)
end

# ----------------------------------------------------------------------------
"""
    _path_delivered!(path_runoff, head_flow, events, trib_output, culvert_actual, cvplan, net,
                     trap_level, nbs_draw) -> Float64

Walk one flow path, applying every stop along it, and return what arrives at the far end.

# Arguments
- `path_runoff`: this path's per-cell oblivious residual — the reference the flow is attenuated
  against ([`_attenuate_range`](@ref)).  Read only.
- `head_flow`: signed flow entering at the path's head.
- `events`: this path's ordered stops, from [`_path_event_templates`](@ref).
- `trib_output`: delivered flow per tributary path, read at each `:trib` stop.
- `culvert_actual`: **mutated** — set at each `:cvin` to the amount actually drawn, read at each
  `:cvout` to deliver it.
- `cvplan`, `net`, `trap_level`: passed through to [`_culvert_flow`](@ref) to price a culvert.
- `nbs_draw`: **mutated** — each `:nbsin` adds the whole passing flow to its placement's entry.

# Returns
The signed flow delivered at the end of the path.

At each stop: a tributary adds its signed output; a culvert outlet adds what its inlet drew; a
culvert inlet draws `min(capacity, available)`; an NBS inflow cell intercepts everything leaving
it, i.e. after its own residual is charged, which is what enters the footprint (and what `O_0`
counts there), and zeroes the path.  Cell `k`'s residual is charged on the segment *leaving* `k`.

Writing the drawn amount to `culvert_actual` and delivering exactly that at the outlet is the
mass-conservation invariant — never the requested capacity.  `nbs_draw` is signed: a negative
draw is an upstream NBS's deficit.
"""
function _path_delivered!(path_runoff, head_flow, events,
                          trib_output, culvert_actual, cvplan, net, trap_level, nbs_draw)
    current = head_flow
    lo = 1
    for (pos, kind, idx) in events
        current = _attenuate_range(path_runoff, lo, pos - 1, current)
        if kind === :trib
            current += trib_output[idx]
        elseif kind === :cvout
            current += culvert_actual[idx]                # deliver
        elseif kind === :cvin
            a = min(_culvert_flow(cvplan, net, idx, trap_level), max(current, 0.0))
            culvert_actual[idx] = a                       # drawn == delivered
            current -= a
        else                                              # :nbsin — intercept ALL flow leaving the
            # inflow cell (charged by its own residual first, as O_0 counts it); signed: an
            # upstream NBS storing shows up as a deficit in this NBS's input Q_in
            nbs_draw[idx] += _attenuate_range(path_runoff, pos, pos, current)
            current = 0.0
        end
        lo = pos
    end
    return _attenuate_range(path_runoff, lo, length(path_runoff), current)
end

"""
    _route_trap_node!(i, net, trap_inflow, path_flow, footprint_infil, spilling, cvplan,
                      trap_level, culvert_actual) -> nothing

Route one trap node: apply the culverts touching it, then spill its surplus if it is full.

# Arguments
- `i`: local trap index.
- `net`: supplies the trap's culvert lists and spill path.
- `trap_inflow`: **mutated** at `i` — culvert deliveries added, culvert draws subtracted.
- `path_flow`: **mutated** at the trap's spill path — a full trap adds its surplus there.
- `footprint_infil`: whole-footprint loss per trap; a full trap's spill is net of it.
- `spilling`: whether each trap spills (is FULL in this network).
- `cvplan`, `trap_level`: passed to [`_culvert_flow`](@ref); `nothing` cvplan skips culverts.
- `culvert_actual`: **mutated** — records what each culvert leaving this trap draws.
- `trap_volume`: per-trap volume, or `nothing`.  An empty trap's culverts draw at most what
  reaches it, whatever their head (a parent's floor sits at its children's spill level).

# Returns
`nothing`; everything lands in the mutated accumulators.

A full trap emits `max(trap_inflow[i] - footprint_infil[i], 0)`.  If its `spill_path` is 0 it
spills out of the domain, so the surplus simply leaves.
"""
function _route_trap_node!(i, net, trap_inflow, path_flow, footprint_infil,
                           spilling, cvplan, trap_level, culvert_actual, trap_volume = nothing)
    if cvplan !== nothing
        trap = net.traps[i]
        for ci in trap.culvert_outlets                # deliver (source drew earlier)
            trap_inflow[i] += culvert_actual[ci]
        end
        empty = trap_volume !== nothing && trap_volume[i] <= 0.0
        for ci in trap.culvert_inlets                 # drain (delivered at its outlet)
            q = _culvert_flow(cvplan, net, ci, trap_level)
            # an empty trap only passes on what reaches it (mass conservation)
            empty && (q = min(q, max(trap_inflow[i], 0.0)))
            culvert_actual[ci] = q
            trap_inflow[i]    -= q
        end
    end
    if spilling[i]
        spill = max(trap_inflow[i] - footprint_infil[i], 0.0)
        sp = net.traps[i].spill_path
        sp > 0 && (path_flow[sp] += spill)            # sp == 0: spills out of domain
    end
    return nothing
end

"""
    _route_flow(net, external_inflow, spilling, footprint_infil, path_runoff,
                path_events, order, merge_target,
                cvplan = nothing, trap_level = nothing,
                path_inject = nothing, trap_inject = nothing; trap_volume = nothing)
                -> (trap_inflow, nbs_draw)

The router: how much water arrives at each trap, given who is spilling right now.

All static routing data is pre-supplied — [`_build_rate_params`](@ref) precomputes it once per
solve onto [`DynNetworkRateParams`](@ref) — so the per-step cost here is just the walk.

# Arguments
- `net`: the network being routed.
- `external_inflow`: constant inflow per trap (terrain runoff); seeds the result.
- `spilling`: whether each trap spills (is FULL in this network), hence emitting.
- `footprint_infil`: whole-footprint loss per trap, deducted from a full trap's spill.
- `path_runoff`: per-path oblivious residual, the attenuation reference.  Read only.
- `path_events`: per-path ordered stops ([`_path_event_templates`](@ref)).
- `order`, `merge_target`: the routing plan ([`_network_order`](@ref), [`_merge_targets`](@ref)).
- `cvplan`: [`CulvertPlan`](@ref), or `nothing` to ignore culverts entirely.
- `trap_level`: per-trap water surface, needed only for culvert heads.
- `path_inject`, `trap_inject`: NBS output added at path heads / into traps
  ([`_nbs_injections`](@ref)), or `nothing` with no NBS.
- `trap_volume` (keyword): per-trap volume, so an empty trap's culverts draw no more than
  reaches it; `nothing` skips that cap.

# Returns
Total inflow at each trap (`net.traps` order), and the flow each NBS placement intercepts at its
inflow cells (`net.nbs` order; signed, as an upstream NBS that is storing arrives as a deficit).

Walks `order`, so everything upstream of a node is final when it is visited.  Nothing the caller
owns is mutated.  Segmented routing
charges flow only the cells it actually travels, so a tributary joining at junction `j` is not
charged the main path's pre-junction losses.
"""
function _route_flow(net::DynNetwork,
                     external_inflow::AbstractVector{<:Real},
                     spilling::AbstractVector{Bool},
                     footprint_infil::AbstractVector{<:Real},
                     path_runoff::AbstractVector{<:AbstractVector{<:Real}},
                     path_events::AbstractVector{<:AbstractVector},
                     order::AbstractVector{<:Integer},
                     merge_target::AbstractVector{<:Integer},
                     cvplan = nothing,
                     trap_level = nothing,
                     path_inject = nothing,
                     trap_inject = nothing;
                     trap_volume = nothing)

    np = length(net.flow_paths)
    nt = length(net.traps)
    @assert length(external_inflow) == length(spilling) == length(footprint_infil) == nt
    @assert length(path_runoff) == length(merge_target) == length(path_events) == np

    # every path's water arrives from an upstream trap's spill or an NBS injection at its head
    trap_inflow = Float64.(external_inflow)
    trap_inject === nothing || (trap_inflow .+= trap_inject)
    path_flow   = zeros(Float64, np)
    trib_output = zeros(Float64, np)
    culvert_actual = cvplan === nothing ? Float64[] : zeros(Float64, length(net.culverts))
    nbs_draw = zeros(Float64, length(net.nbs))

    for node in order
        if node <= np                               # a flow path
            p = node
            head = path_inject === nothing ? path_flow[p] : path_flow[p] + path_inject[p]
            delivered = _path_delivered!(path_runoff[p], head,
                                         path_events[p], trib_output, culvert_actual,
                                         cvplan, net, trap_level, nbs_draw)
            tt = net.flow_paths[p].target_trap
            if tt > 0
                trap_inflow[tt] += delivered        # into the downstream trap
            elseif merge_target[p] > 0
                trib_output[p] = delivered          # buffer for the main path's routing
            end                                      # else: exits the domain
        else                                        # a trap
            _route_trap_node!(node - np, net, trap_inflow, path_flow, footprint_infil,
                              spilling, cvplan, trap_level, culvert_actual, trap_volume)
        end
    end
    return trap_inflow, nbs_draw
end

# ----------------------------------------------------------------------------
"""
    _footprint_infiltration(geom) -> Vector{Float64}

The loss a full trap incurs, when its whole footprint is submerged.

# Arguments
- `geom`: per-trap geometry, from [`_build_trap_geometry`](@ref).

# Returns
Per-cell `infil` summed over each trap's footprint, in `geom` order — the ceiling that
[`wetted_infiltration`](@ref) approaches as `V` reaches capacity.  Nothing is mutated.
"""
function _footprint_infiltration(geom::Vector{TrapGeometry})
    return [sum(g.infil) for g in geom]
end

# ----------------------------------------------------------------------------
"""
    _path_cell_values(net, grid) -> Vector{Vector{Float64}}

Sample a terrain grid along each flow path's cells.

# Arguments
- `net`: supplies the paths and their cells.
- `grid`: any per-cell terrain grid.

# Returns
One vector per flow path, in `net.flow_paths` order, holding `grid` at that path's cells; empty
for a zero-length connector.  Nothing is mutated.

The caller decides what the values mean.  `_build_rate_params` passes the oblivious runoff grid
to get the residual [`_attenuate_range`](@ref) charges all flow against, or — with no such grid —
negated infiltration, which puts every cell at its capacity floor and reproduces the same rule.
"""
function _path_cell_values(net::DynNetwork, grid::AbstractMatrix{<:Real})
    return [isempty(p.cells) ? Float64[] : Float64[grid[c] for c in p.cells]
            for p in net.flow_paths]
end

# ============================================================================
# Network rate function.
#
# `dynNetworkRateFunction!` is the in-place ODE RHS over the trap-volume state.  Its
# `DynNetworkRateParams` bundle everything static for a solve (geometry, external inflow,
# infiltration sums, routing plan), so per-step work is just routing + a loop over traps.
# ============================================================================

# ============================================================================
# NBS in the solver (model and invariants: agent/NBS_INVARIANTS.md).  Each placement adds its
# layer states after the trap volumes, injects (O_t − O_0)·share at its exits and ponds and its
# piped overflow at its outlets (`_nbs_injections`), and is fed O_0 plus the `:nbsin` draws the
# routing pass collects.
# ============================================================================

"""
    NBSElement

One NBS placement resolved into router terms for a solve, built by [`_build_nbs_plan`](@ref).
`O_0_total` and the shares come from [`_nbs_baseline`](@ref).
"""
struct NBSElement
    system    ::NBSSystem
    A         ::Float64                        # footprint area (m^2)
    state_base::Int                            # layers occupy state nt+state_base+1 ... (top first)
    n_terrain ::Int                            # top layers overflowing onto the terrain
    O_0_total ::Float64                        # oblivious throughput
    terrain_paths::Vector{Tuple{Int,Float64}}  # (carrier path, share) per exit
    terrain_traps::Vector{Tuple{Int,Float64}}  # (trap, share) per pond
    piped_paths  ::Vector{Tuple{Int,Int}}      # (carrier path, layer) per piped outlet
end

"""
    NBSPlan

Every NBS placement of a solve, resolved into router terms.

# Fields
- `elems`: one [`NBSElement`](@ref) per placement, in `net.nbs` order.
- `nlayer_total`: layer states appended after the trap volumes — how much longer the ODE state
  vector is than the trap count.
"""
struct NBSPlan
    elems       ::Vector{NBSElement}
    nlayer_total::Int
end

"""
    _nbs_state_count(p) -> Int

How many NBS layer states the ODE state vector carries beyond the trap volumes.

# Arguments
- `p`: the [`DynNetworkRateParams`](@ref) for this solve.

# Returns
`p.nbsplan.nlayer_total`, or `0` when the net has no NBS.
"""
_nbs_state_count(p) = p.nbsplan === nothing ? 0 : p.nbsplan.nlayer_total

"""
    _attenuate_range(runoff, lo, hi, d) -> Float64

The one flow-tracking rule, applied to every flow the router carries — network spills,
tributary and culvert additions, NBS output diffs alike.

# Arguments
- `runoff`: per-cell oblivious residual for this path.  Read only.
- `lo`, `hi`: the inclusive cell range to cross (`hi < lo` crosses nothing).
- `d`: the signed flow entering the range.

# Returns
What survives past cell `hi`, each cell applying `max(V+d, 0) - max(V, 0)`.

A cell with spare capacity (`V < 0`) infiltrates the flow until it vanishes; a cell already
carrying background flow (`V >= 0`) passes it unchanged; negative `d` is symmetric.  Charging
the *residual* rather than the raw infiltration is what avoids re-charging a cell the background
already saturated.  Mirrors `_track_flow!` (`flow.jl`), which builds the very grid it reads.
"""
function _attenuate_range(runoff, lo::Int, hi::Int, d::Float64)
    @inbounds for k in lo:hi
        V = runoff[k]
        d = max(V + d, 0.0) - max(V, 0.0)
    end
    return d
end

"""
    _attenuate_diff(runoff, d) -> Float64

[`_attenuate_range`](@ref) over a whole path.

# Arguments
- `runoff`: the path's per-cell oblivious residual.  Read only.
- `d`: the signed flow to attenuate.

# Returns
What survives to the end of the path.  Convenience form for full-vector callers and the NBS
unit tests.
"""
_attenuate_diff(runoff::AbstractVector{<:Real}, d::Float64) =
    _attenuate_range(runoff, 1, length(runoff), d)

"""
    _build_nbs_plan(net, tstruct, runoff, precipitation) -> NBSPlan or nothing

Resolve every NBS placement in `net` into router terms (`nothing` if there are none): its
baseline throughput and shares ([`_nbs_baseline`](@ref)), the carrier path leaving each exit and
piped outlet, and the trap under each pond.  Nothing is mutated.

Every exit and outlet is a network seed, so a path departs from it; that is asserted.
"""
function _build_nbs_plan(net::DynNetwork, tstruct, runoff::AbstractMatrix{<:Real},
                         precipitation::Union{AbstractMatrix{<:Real},Real})
    isempty(net.nbs) && return nothing
    departing = Dict{CartesianIndex{2},Int}()           # cell -> first path departing from it
    for (p, fp) in enumerate(net.flow_paths)
        get!(departing, fp.departure_point, p)
    end
    function carrier(cell)
        p = get(departing, cell, 0)
        @assert p != 0 "NBS exit/outlet $cell has no carrier path (expected it to be a network seed)"
        return p
    end

    elems = NBSElement[]
    base  = 0
    for nb in net.nbs
        O_0, exit_shares, pond_shares = _nbs_baseline(nb, tstruct, runoff, precipitation)
        push!(elems, NBSElement(nb.system,
                                Float64(length(nb.footprint)),   # @@@ 1 m^2/cell; use the real cell area
                                base, nb.n_terrain, O_0,
                                [(carrier(c), s) for (c, s) in zip(nb.link.exits, exit_shares)],
                                [(_pond_trap(net, tstruct, c), s) for (c, s) in zip(nb.link.ponds, pond_shares)],
                                [(carrier(o), l) for (o, l) in zip(nb.outlets, _piped_layers(nb.system, nb.n_terrain))]))
        base += length(nb.system.layers)
    end
    return NBSPlan(elems, base)
end

"""
    _nbs_injections(V, plan, nt, np) -> (path_inject, trap_inject)

What every NBS placement adds to the routing at state `V`, computed from the layer storage
alone, before the routing pass.  Each exit's carrier path gets `(O_terrain - O_0)·share` at its
head and each pond's trap the same share term; each piped outlet's carrier path gets that layer's
overflow.  Adding the *difference* to `O_0` is what keeps `external_inflow`, which already holds
`O_0`, from being counted twice.
"""
function _nbs_injections(V, plan::NBSPlan, nt::Int, np::Int)
    path_inject = zeros(Float64, np)
    trap_inject = zeros(Float64, nt)
    for el in plan.elems
        Vs   = view(V, _layer_range(el, nt))
        diff = _nbs_terrain_outflow(el.system, Vs, el.A, el.n_terrain) - el.O_0_total
        for (pth, sh) in el.terrain_paths; path_inject[pth] += diff * sh; end
        for (tr,  sh) in el.terrain_traps; trap_inject[tr]  += diff * sh; end
        for (pth, l) in el.piped_paths
            path_inject[pth] += first(_layer_fluxes(el.system.layers[l], Vs[l], el.A))
        end
    end
    return path_inject, trap_inject
end

# ----------------------------------------------------------------------------
"""
    DynNetworkRateParams

Static precomputed inputs to [`dynNetworkRateFunction!`](@ref) for one solve, built by
[`_build_rate_params`](@ref) and passed as the `ODEProblem` parameter.  Indexed
network-locally (`net.traps` / `net.flow_paths` order).

# Fields
- `net`: the [`DynNetwork`](@ref) being solved
- `geom`: per-trap [`TrapGeometry`](@ref)
- `external_inflow`: constant inflow per trap.  This already accounts for rain landing on the
  flow paths: a path cell lies in its downstream trap's spill region, so `watercourses`
  accumulates it into that trap's inflow.  Paths therefore carry no inflow of their own.
- `footprint_infil`: whole-footprint infiltration per trap
- `path_runoff`: per-path oblivious residual runoff, the read-only reference the router
  attenuates all flow against (real grid where available, else `-infiltration`)
- `order`, `merge_target`: the static routing plan ([`_network_order`](@ref), [`_merge_targets`](@ref))
- `path_events`: per-path in-order `(cell_pos, kind, idx)` stops — tributary junctions, plus
  culvert inlet/outlet and NBS-inlet positions when the net has culverts / NBS
- `cvplan`: culvert routing data (`nothing` if no culverts)
- `nbsplan`: NBS signed-diff routing data (`nothing` if no NBS)
"""
struct DynNetworkRateParams
    net::DynNetwork                       # The network being solved
    geom::Vector{TrapGeometry}            # struct holding info about trap's geometry
    external_inflow::Vector{Float64}      # constant inflow from terrain/weather per trap
    footprint_infil::Vector{Float64}      # total footprint infiltration per trap (full trap loss)
    path_runoff::Vector{Vector{Float64}}  # per-path oblivious residual runoff, the read-only reference
                                          # the router attenuates all flow against (real grid where
                                          # available, else -infiltration = cells at their capacity floor)
    order::Vector{Int}                    # topological sort of the flowgraph
    merge_target::Vector{Int}             # path that each tributary feeds into (0 if none)
    cvplan::Union{CulvertPlan,Nothing}    # key information per culvert for routing, or nothing if none
    path_events::Vector{Vector{Tuple{Int,Symbol,Int}}}  # in-order tributary/culvert/NBS-inlet stops per path
    nbsplan::Union{NBSPlan,Nothing}      # NBS signed-diff routing data, or nothing if none
end

# ----------------------------------------------------------------------------
"""
    _build_rate_params(tstruct, net, infiltration, external_inflow;
                       runoff=nothing, precipitation=nothing, zvt=nothing) -> DynNetworkRateParams

Precompute everything static for one solve, so per-step work is just routing plus a loop over
traps.

# Arguments
- `tstruct`, `net`: the terrain and the network to solve.
- `infiltration`: per-cell infiltration-rate grid (0 = impermeable).
- `external_inflow`: constant inflow per trap, `net.traps` order.
- `runoff`: the oblivious runoff grid.  Used for `path_runoff`, and required for the NBS plan.
- `precipitation`: the period's rain rate; required only when the net has NBS.
- `zvt`: cached volume↔level tables ([`_compute_z_vol_tables`](@ref)).

# Returns
The [`DynNetworkRateParams`](@ref) the ODE reads as its `p`.  Nothing is mutated.

Errors if `net` has NBS but `runoff` or `precipitation` is missing.  A net with no culverts gets
`cvplan = nothing`; with no NBS, `nbsplan = nothing`.
"""
function _build_rate_params(tstruct::TrapStructure,
                            net::DynNetwork,
                            infiltration::AbstractMatrix{<:Real},
                            external_inflow::AbstractVector{<:Real};
                            runoff::Union{AbstractMatrix{<:Real},Nothing} = nothing,
                            precipitation::Union{AbstractMatrix{<:Real},Real,Nothing} = nothing,
                            zvt = nothing)
    tgeom = _build_trap_geometry(tstruct, net, infiltration; zvt=zvt)
    structure = DynNetworkRateParams(net, tgeom, Float64[], _footprint_infiltration(tgeom),
                                     Vector{Float64}[], first(_network_order(net)),
                                     _merge_targets(net),
                                     isempty(net.culverts) ? nothing : _build_culvert_plan(net, tstruct),
                                     _path_event_templates(net), nothing)
    return _with_flows(structure, tstruct, infiltration, external_inflow;
                       runoff = runoff, precipitation = precipitation)
end

# `p` with its flow-dependent parts recomputed: external inflow, path residuals and NBS plan.  The
# rest depends only on the network and the period's infiltration, so a context whose network did
# not change keeps it across a touch.
function _with_flows(p::DynNetworkRateParams, tstruct, infiltration, external_inflow;
                     runoff = nothing, precipitation = nothing)
    net = p.net
    @assert length(external_inflow) == length(net.traps)
    # The reference the router attenuates all flow against.  With a real grid it is
    # background-aware (a cell the background saturated charges the network flow nothing); without
    # one (hand-built nets) cells sit at their capacity floor (-infiltration).
    path_runoff = runoff === nothing ?
        [Float64[-c for c in ci] for ci in _path_cell_values(net, infiltration)] :
        _path_cell_values(net, runoff)
    if !isempty(net.nbs)
        runoff === nothing &&
            error("_build_rate_params: net has NBS but no runoff grid was supplied")
        precipitation === nothing &&
            error("_build_rate_params: net has NBS but no precipitation was supplied " *
                  "(needed for each placement's throughput: rain on the footprint + inflow)")
    end
    nbsplan = isempty(net.nbs) ? nothing : _build_nbs_plan(net, tstruct, runoff, precipitation)
    return DynNetworkRateParams(net, p.geom, Float64.(external_inflow), p.footprint_infil,
                                path_runoff, p.order, p.merge_target, p.cvplan, p.path_events,
                                nbsplan)
end

# ----------------------------------------------------------------------------
"""
    _routed_inflow(V, p) -> (inflow, spilling, nbs_draw)

Route the network at a given state — the single point where "what is arriving where" is decided.

# Arguments
- `V`: the ODE state — trap volumes, then any NBS layer states.
- `p`: the [`DynNetworkRateParams`](@ref) for this solve.

# Returns
`(inflow, spilling, nbs_draw)` — total inflow per trap, whether each trap spills, and the flow
each NBS placement intercepts at its inflow cells.  A trap spills when it is FULL in this
network (it has a spill path), not when `V` reaches capacity: a filling trap keeps its own rate
through capacity until the `:fill` event ends the solve, so the rate has no kink there for the
event location to fall short of (a stage trying `V > C` must not stop the trap filling).

Shared by the rate function and the `:unspill` condition, so the two cannot disagree about who
is spilling.  `nbs_draw` is complete only because every path crossing a placement's inflow cells
is in the placement's component (`_nbs_coupled_nodes`; invariant C1 in NBS_INVARIANTS.md).
"""
function _routed_inflow(V, p::DynNetworkRateParams)
    geom = p.geom
    nt   = length(geom)
    spilling = Bool[p.net.traps[i].spill_path != 0 for i in 1:nt]
    # culverts need each trap's water-surface elevation for the head calc
    trap_level = p.cvplan === nothing ? nothing :
                 Float64[water_level(geom[i], V[i]) for i in 1:nt]
    # NBS output depends on the layer state only, so it is fixed before the routing pass
    inject = p.nbsplan === nothing ? (nothing, nothing) :
             _nbs_injections(V, p.nbsplan, nt, length(p.net.flow_paths))
    inflow, nbs_draw = _route_flow(p.net, p.external_inflow, spilling,
                                   p.footprint_infil, p.path_runoff,
                                   p.path_events, p.order, p.merge_target,
                                   p.cvplan, trap_level, inject...; trap_volume = V)
    return inflow, spilling, nbs_draw
end

"""
    dynNetworkRateFunction!(dV, V, p::DynNetworkRateParams, t=0.0)

In-place ODE rate of the trap-volume state `V` (plus appended NBS layer states).

A *full* trap (`V >= capacity`) passes excess downstream at steady volume:
`dV = inflow - footprint_infil - max(inflow - footprint_infil, 0)` (0 while well fed, negative
once inflow drops below its losses).  An *accumulating* trap fills at its wetted-area rate:
`dV = inflow - wetted_infiltration(V)`.  `inflow` is the external inflow plus everything routed
in from upstream (see [`_route_flow`](@ref)).

# Arguments
- `dV`: **written in place** — the rate for each trap, then for each NBS layer.  The only
  argument mutated; `V` and `p` are read only.
- `V`: current state — one volume per trap in `p.net.traps` order, then the NBS layer states.
- `p`: the static parameters for this solve (geometry, external inflow, routing plan).
- `t`: unused — the rate is autonomous; present for the `DifferentialEquations.jl` RHS
  signature.

# Returns
`nothing`; the rates land in `dV`.
"""
function dynNetworkRateFunction!(dV, V, p::DynNetworkRateParams, t = 0.0)
    geom = p.geom
    nt   = length(geom)

    inflow, spilling, nbs_draw = _routed_inflow(V, p)

    @assert length(dV) == length(V) == nt + _nbs_state_count(p)
    @inbounds for i in 1:nt
        if spilling[i]
            # full: whole footprint wetted; excess overflows, so dV is ~0 while fed
            loss  = p.footprint_infil[i]
            spill = max(inflow[i] - loss, 0.0)
        else
            # accumulating: infiltrate only through the wetted footprint, no spill
            loss  = wetted_infiltration(geom[i], V[i])
            spill = 0.0
        end
        dV[i] = inflow[i] - loss - spill
        # Physical floor: volume cannot go below zero.  When an EMPTY trap has negative
        # computed dV (e.g. a culvert trying to draw from it), clamp to zero rather than
        # letting the ODE push V < 0 and blow up the v2z interpolation.
        V[i] <= 0.0 && dV[i] < 0.0 && (dV[i] = 0.0)
    end

    # NBS layers, fed at the top by Q_in = O_0_total + the :nbsin draws.  Per inflow cell c the two
    # add up to max(R_c + d_c, 0), so Q_in >= 0 exactly; max() only absorbs rounding.
    if p.nbsplan !== nothing
        for (k, el) in enumerate(p.nbsplan.elems)
            r = _layer_range(el, nt)
            _nbs_layer_rates!(view(dV, r), el.system, view(V, r), el.A,
                              max(el.O_0_total + nbs_draw[k], 0.0))
        end
    end
    return nothing
end

# State indices of `el`'s layers (top first), after the `nt` trap volumes.
_layer_range(el::NBSElement, nt::Int) = (nt + el.state_base + 1):(nt + el.state_base + length(el.system.layers))

# ============================================================================
# Event detection.  Two callbacks stop the integration at the first event:
#
# 1. cb_topo (VectorContinuousCallback, LeftRootFind): topology events on downward
#    zero-crossings —
#      :fill    capacity - V   trap fills, starts spilling
#      :empty   V              trap empties, exposes its subtraps
#      :unspill inflow - loss  a full trap's net inflow goes negative
#    Evolving traps get :fill (parents also :empty); full traps get :unspill.  LeftRootFind lets
#    conditions starting at 0 (full traps, zero initial net inflow) wait for a genuine crossing
#    rather than firing degenerately.
#
# 2. cb_ss (DiscreteCallback): steady state, checked at accepted steps — every evolving trap at
#    rate exactly 0 or past a sign change (stalled at a wetted-area step), every NBS layer below a
#    relative rate tolerance.  Discrete, not Continuous: wetted_infiltration is a step function of
#    V, so interpolated mid-step states would give spurious crossings.  Vetoes when any evolving
#    trap is at V >= C (:fill wins).  Only when predicting: a solve that must reach `tmax`
#    (commit, sampling) never stops early, or the remaining flow would be lost.
# ============================================================================

"""
    EventCondition

One topology event the solve watches for.  The `VectorContinuousCallback` monitors one scalar
condition per entry.

# Fields
- `kind`: `:fill`, `:empty` or `:unspill`.
- `trap`: the network-local trap index it concerns.
"""
struct EventCondition
    kind::Symbol
    trap::Int
end

"""
    DynNetworkEvent

The event a solve terminated on.  Mutable, so the callback can record the hit and
[`solveDynNetwork!`](@ref) read it back after integrating.

# Fields
- `kind`: `:fill`, `:empty`, `:unspill`, or `:none` (steady state or `tmax` cutoff).
- `trap`: the network-local trap it concerns, `0` if none.  `solveDynNetwork!` maps this to a
  global index before returning.
"""
mutable struct DynNetworkEvent
    kind::Symbol
    trap::Int
end
DynNetworkEvent() = DynNetworkEvent(:none, 0)

# ----------------------------------------------------------------------------
"""
    _event_conditions(p, evolving, nreg) -> Vector{EventCondition}

The topology events to monitor for this solve.

# Arguments
- `p`: the [`DynNetworkRateParams`](@ref), for the traps.
- `evolving`: local indices of the traps below capacity; the rest are full.
- `nreg`: region count — a trap index above it is a parent, at or below it a leaf.

# Returns
A `Vector{EventCondition}`: `:fill` per evolving trap, `:empty` per evolving *parent*, and
`:unspill` per full trap.  The order fixes the callback's condition indices.

A leaf reaching `V = 0` is just sitting at its floor rather than changing topology, so only
parents get `:empty`.  Nothing is mutated.
"""
function _event_conditions(p::DynNetworkRateParams,
                           evolving::AbstractVector{<:Integer},
                           nreg::Int)
    conds = EventCondition[]
    for e in evolving
        push!(conds, EventCondition(:fill,  e))
        p.net.traps[e].trap_ix > nreg &&
            push!(conds, EventCondition(:empty, e))
    end
    evolving_set = Set(evolving)
    for i in eachindex(p.net.traps)
        i in evolving_set && continue
        push!(conds, EventCondition(:unspill, i))
    end
    return conds
end

# ----------------------------------------------------------------------------
"""
    _build_steadystate_callback(p, evolving, du0) -> DiscreteCallback

Steady-state detector for an NBS-free network.

# Arguments
- `p`: the [`DynNetworkRateParams`](@ref) for this solve.
- `evolving`: local indices of the traps to watch.
- `du0`: each trap's rate at `t = 0`, for the sign-change test.

# Returns
A `DiscreteCallback` that terminates at the first accepted step where every evolving trap has
settled — its rate is exactly 0, or has crossed the sign it started at, meaning a trap
oscillating across a wetted-infiltration step has reached a sub-capacity equilibrium (the
multi-trap analogue of `fill_trap_until`'s stagnation stop).  A small but nonzero rate is not
settled: the trap is still filling or draining, however slowly.

Discrete rather than continuous because `wetted_infiltration` is a step function of `V`, so
interpolated mid-step states would give spurious crossings.  Vetoes while any evolving trap sits
at `V >= C`, letting `:fill` terminate instead.
"""
function _build_steadystate_callback(p::DynNetworkRateParams,
                                     evolving::AbstractVector{<:Integer},
                                     du0::AbstractVector{<:Real})
    function condition(u, t, integrator)
        for e in evolving
            u[e] >= p.geom[e].capacity && return false   # veto: :fill wins
        end
        inflow, _ = _routed_inflow(u, p)
        for e in evolving
            rate = inflow[e] - wetted_infiltration(p.geom[e], u[e])
            # Match the rate function's physical floor: an empty trap cannot drain
            # further, so its actual dV is 0 (otherwise cb_ss never sees a leaf that has
            # drained to V == 0 settle — its unclamped rate stays negative and the solve
            # marches to t → ∞).
            u[e] <= 0.0 && rate < 0.0 && (rate = 0.0)
            # Settled if the actual rate is 0, OR it has crossed its starting sign — a
            # trap oscillating across a wetted-infiltration step has reached a
            # sub-capacity equilibrium and makes no further net progress (the multi-trap
            # analogue of `fill_trap_until`'s stagnation/sign-change stop).
            _trap_settled(rate, du0[e]) || return false
        end
        return true
    end
    return DiscreteCallback(condition, terminate!; save_positions = (false, false))
end

"""
    _build_steadystate_callback_nbs(p, ss_indices, layer_tol, du0) -> DiscreteCallback

Steady-state detector for a network with NBS.

# Arguments
- `p`: the [`DynNetworkRateParams`](@ref) for this solve.
- `ss_indices`: state indices to watch — the evolving traps plus the NBS layer states.
- `layer_tol`: the `|dS/dt|` below which an NBS layer counts as settled.
- `du0`: each state's rate at `t = 0`, for the traps' sign-change test.

# Returns
A `DiscreteCallback`, as [`_build_steadystate_callback`](@ref) but evaluating the *full* rate
function: the cheap `_routed_inflow` check cannot express the layer dynamics.  Traps settle as
there.  A layer approaches its equilibrium exponentially, never reaching a zero rate, so it
settles below `layer_tol`; a sign change does not settle it, since a layer can drain and then
refill.  A trap at capacity still vetoes, letting `:fill` terminate.
"""
function _build_steadystate_callback_nbs(p::DynNetworkRateParams,
                                         ss_indices::AbstractVector{<:Integer},
                                         layer_tol::Real,
                                         du0::AbstractVector{<:Real})
    nt        = length(p.geom)
    trap_test = Int[i for i in ss_indices if i <= nt]
    dbuf      = zeros(Float64, length(du0))
    function condition(u, t, integrator)
        for e in trap_test
            u[e] >= p.geom[e].capacity && return false   # veto: :fill wins
        end
        dynNetworkRateFunction!(dbuf, u, p, t)
        for i in ss_indices
            settled = i <= nt ? _trap_settled(dbuf[i], du0[i]) : abs(dbuf[i]) < layer_tol
            settled || return false
        end
        return true
    end
    return DiscreteCallback(condition, terminate!; save_positions = (false, false))
end

# A trap's rate is settled at exactly 0, or once it has crossed the sign `rate0` it started at.
_trap_settled(rate, rate0) = rate == 0.0 || (rate0 != 0.0 && signbit(rate) != signbit(rate0))

# Relative rate tolerance for NBS layers in the steady-state test, scaled by the largest
# starting rate of the solve (`abstol` is a volume, not a rate).
const _NBS_SETTLED_RTOL = 1e-6

# ----------------------------------------------------------------------------
"""
    _build_event_callback(p, evolving, V0, nreg) -> (callback, event)

The callback that halts integration at the first topology-changing event.

# Arguments
- `p`: the [`DynNetworkRateParams`](@ref) for this solve.
- `evolving`: local indices of the non-FULL traps; the rest get `:unspill` instead.
- `V0`: the entry state, for the conditions' initial values.
- `nreg`: region count, restricting `:empty` to parent traps.

# Returns
`(callback, event)` — a `VectorContinuousCallback` (LeftRootFind), and the
[`DynNetworkEvent`](@ref) it writes the hit into.  The caller reads `event` after `solve`
returns; it stays `:none` if nothing fired.

LeftRootFind lets a condition that starts *at* zero — a full trap with zero net inflow — wait
for a genuine crossing instead of firing degenerately.
"""
function _build_event_callback(p::DynNetworkRateParams,
                               evolving::AbstractVector{<:Integer},
                               V0::AbstractVector{<:Real},
                               nreg::Int)

    conds = _event_conditions(p, evolving, nreg)
    event = DynNetworkEvent()

    function condition(out, V, t, integrator)
        # `:unspill` conditions need the routed inflow; `:fill`/`:empty` read V directly.
        # (Route unconditionally — a big fused component almost always has a full trap.)
        inflow, _ = _routed_inflow(V, p)
        @inbounds for (k, ec) in enumerate(conds)
            if ec.kind == :fill
                out[k] = p.geom[ec.trap].capacity - V[ec.trap]
            elseif ec.kind == :empty
                out[k] = V[ec.trap]
            else   # :unspill: fire when the full trap's net inflow drops below its losses
                out[k] = inflow[ec.trap] - p.footprint_infil[ec.trap]
            end
        end
    end

    function affect!(integrator, ix)
        k  = isa(ix, AbstractVector) ? findfirst(!iszero, ix) : ix
        ec = conds[k]
        event.kind = ec.kind
        event.trap = ec.trap
        terminate!(integrator)
    end

    return VectorContinuousCallback(condition, affect!, length(conds);
                                    rootfind = SciMLBase.LeftRootFind), event
end

# ----------------------------------------------------------------------------
"""
    _validate_network(tstruct, net, state, geom) -> nothing

Enforce the three-state caller contract (FULL / TRANSITORY / EMPTY — see
[`solveDynNetwork!`](@ref)) before a solve.

# Arguments
- `tstruct`: supplies the trap hierarchy the parent/child checks walk.
- `net`: the network whose traps are checked.
- `state`: the entry volumes, `net.traps` order.
- `geom`: per-trap geometry, supplying each capacity.

# Returns
`nothing`, or throws an informative error on the first violation.  Nothing is mutated.

Checks that a FULL trap has a spill path and that a FULL parent's in-network children are FULL
too, and that a TRANSITORY trap has *no* spill path.  A parent at `V == 0` is TRANSITORY, not
EMPTY: it subsumes its full children, so zero is its own-volume floor with them submerged.
"""
function _validate_network(tstruct::TrapStructure,
                           net::DynNetwork,
                           state::AbstractVector{<:Real},
                           geom::Vector{TrapGeometry})
    nreg     = numregions(tstruct)
    nt       = length(net.traps)
    trap_map = Dict(net.traps[i].trap_ix => i for i in 1:nt)
    for i in 1:nt
        trap = net.traps[i]
        V    = state[i]
        C    = geom[i].capacity
        if V == C                                     # FULL
            trap.spill_path != 0 ||
                error("solveDynNetwork!: FULL trap $(trap.trap_ix) has spill_path == 0. " *
                      "Rebuild the network (via setup_network) with this trap in full_traps " *
                      "before the next call so it obtains a downstream spill path " *
                      "(spill_path > 0), or the out-of-domain sentinel (spill_path == -1) " *
                      "if it spills straight out of the domain.")
            if trap.trap_ix > nreg                    # parent: in-network children must be FULL
                for child in subtrapsof(tstruct, trap.trap_ix)
                    ci = get(trap_map, child, 0)
                    ci == 0 && continue               # child not in network; assumed FULL externally
                    state[ci] == geom[ci].capacity ||
                        error("solveDynNetwork!: FULL parent trap $(trap.trap_ix) has " *
                              "in-network child $child with V < C. " *
                              "A parent trap can only hold water above both children's spillpoints, " *
                              "so both children must be FULL whenever the parent is FULL.")
                end
            end
        elseif V == 0.0 && trap.trap_ix <= nreg       # EMPTY (leaf trap at its floor)
            # a lowest-level trap may legitimately sit empty; nothing to check.
        else                                          # TRANSITORY
            # also covers a parent at its floor (V==0): it subsumes its full children, so V==0 is
            # the parent's own-volume floor with children submerged — TRANSITORY, not EMPTY.
            trap.spill_path == 0 ||
                error("solveDynNetwork!: TRANSITORY trap $(trap.trap_ix) has spill_path > 0. " *
                      "A trap that is not yet full should not have a downstream spill path. " *
                      "Rebuild the network without this trap in full_traps.")
        end
    end
end

"""
    _t0_fast_path(V0, du0, p, nreg) -> (; time, trap, kind) or nothing

Catch an entry state that cannot be held, so the caller can skip integration entirely.

# Arguments
- `V0`: the entry state.
- `du0`: the rate at `V0`, from one rate-function evaluation.
- `p`: the [`DynNetworkRateParams`](@ref) for this solve.
- `nreg`: region count, distinguishing a parent from a leaf.

# Returns
`(; time = 0.0, trap, kind)` with `trap` a **global** index, or `nothing` when the state is
viable and the ODE should run.  Nothing is mutated.

A FULL trap already draining (`du0 < 0`) unspills at once; a parent at its floor (`V == 0`) with
negative net inflow must expose its children now.  The parent check recomputes the *unclamped*
net rate, since `du0`'s physical floor guard has already zeroed it at `V == 0`.
"""
function _t0_fast_path(V0, du0, p::DynNetworkRateParams, nreg::Int)
    nt = length(p.geom)
    for i in 1:nt
        V0[i] == p.geom[i].capacity && du0[i] < 0.0 || continue
        return (time = 0.0, trap = p.net.traps[i].trap_ix, kind = :unspill)
    end
    inflow0, _ = _routed_inflow(V0, p)
    for i in 1:nt
        (V0[i] == 0.0 && p.net.traps[i].trap_ix > nreg) || continue
        inflow0[i] - wetted_infiltration(p.geom[i], 0.0) < 0.0 || continue
        return (time = 0.0, trap = p.net.traps[i].trap_ix, kind = :empty)
    end
    return nothing
end

# ============================================================================
# Driver: solveDynNetwork!.
#
# Wires the pieces together: validate the caller contract, build the static rate
# parameters, classify traps into evolving / full, run the ODE, and stop at the
# first topology-changing event.  The `state` vector is updated in place.
# ============================================================================

"""
    solveDynNetwork!(state, tstruct, net, infiltration, inflow;
                     tmax=Inf, runoff=nothing, precipitation=nothing,
                     abstol=1e-6, reltol=1e-4, stop_at_steady=true, stop_at=nothing,
                     zvt=nothing)
                     -> (; time, trap, kind)

Evolve the trap volumes of a [`DynNetwork`](@ref) forward until the first topology-changing
event, steady state, or `tmax`.  `state` is updated in place.

# Arguments
- `state`: stored volume (net of subtraps) per trap, `net.traps` order.  **Mutated in place**;
  must satisfy the three-state contract (below) on entry.
- `tstruct`: the [`TrapStructure`](@ref) supplying `net`'s geometry
- `net`: the network to solve
- `infiltration`: per-cell infiltration-rate grid (0 = impermeable)
- `inflow`: constant inflow per trap, `net.traps` order.  Flow paths take no inflow of their
  own: rain landing on a path cell is already in its downstream trap's entry here, because the
  cell lies in that trap's spill region and `watercourses` accumulates it there.  Everything
  else reaching a path — an upstream trap's spill, a culvert delivery, an NBS output — is
  routed by the solver itself.

# Keyword arguments
- `tmax`: integrate at most this far (elapsed).  With no event by then, stop at `tmax` and
  return `kind = :none` with `state` at `tmax` — lets a caller read the volumes at an
  externally-fixed time.  Default `Inf` (run to the first event or steady state).
- `runoff`: oblivious runoff grid, needed only when `net` has NBS
- `precipitation`: the period's rain rate (scalar or per-cell), needed only when `net` has NBS —
  each placement's throughput is the rain on its footprint plus the flow entering it
- `zvt`: cached volume↔level tables ([`_compute_z_vol_tables`](@ref))
- `abstol`, `reltol`: ODE tolerances
- `stop_at_steady`: stop early at steady state, returning `time = Inf` (default `true`, for
  predicting the next event).  Pass `false` when the state at `tmax` is needed (a commit, a
  timepoint sample): the solve then always runs to `tmax`, so slow flow is not dropped.
- `stop_at`: stop at this elapsed time (`<= tmax`) instead of `tmax`, returning `time = stop_at`.
  The step sequence is then that of a solve to `tmax` (the solver's steps depend on the window),
  so a commit can reproduce its prediction exactly.  Default `nothing` (stop at `tmax`).

# Returns
`(; time, trap, kind)` (state is in place): `time` the event time (`Inf` = steady state,
`tmax` = window elapsed); `trap` the global trap index or 0; `kind` one of `:fill`, `:empty`,
`:unspill`, `:none` (`:none` covers steady state and a `tmax` cutoff — tell apart by `time`).

# Three-state contract (validated at entry)
- **FULL** (V == C): `spill_path != 0` (`> 0` into an in-network path, `-1` out of domain); a
  parent's in-network children are all FULL too.
- **TRANSITORY** (0 < V < C): `spill_path == 0`.
- **EMPTY** (V == 0): a leaf trap (`trap_ix <= numregions(tstruct)`).

# Caller protocol after each event

| Event       | Action before next call                                               |
|-------------|-----------------------------------------------------------------------|
| `:fill`     | `state[trap]` already clamped to C.  Add trap to `full_traps`,       |
|             | rebuild network via `setup_network`.                                  |
| `:unspill`  | `state` unchanged (trap is still at C).  Remove trap from            |
|             | `full_traps`, set `state[trap] = prevfloat(C)`, rebuild network.     |
| `:empty`    | `state[trap]` clamped to 0.  Remove parent + all its children from   |
|             | `full_traps`.  Set `state[child] = prevfloat(C_child)`.  Rebuild.   |
| `:none`     | Steady state; no further event.                                       |

A FULL trap that stays adequately fed keeps `dV == 0` exactly (the spill term cancels
`inflow - loss`), so it remains at exactly C — no drift, nothing for the caller to clamp.
A FULL trap whose inflow falls below its losses fires `:unspill` instead.
"""
function solveDynNetwork!(state::AbstractVector{Float64},
                          tstruct::TrapStructure,
                          net::DynNetwork,
                          infiltration::AbstractMatrix{<:Real},
                          inflow::AbstractVector{<:Real};
                          tmax = Inf,
                          runoff::Union{AbstractMatrix{<:Real},Nothing} = nothing,
                          precipitation::Union{AbstractMatrix{<:Real},Real,Nothing} = nothing,
                          # 1e-6 m^3 (~mL) / 1e-4: enough for physical accuracy, and ~halves the
                          # step count vs 1e-8 on the culvert worst case (see PARITY_TOL in tests)
                          abstol = 1e-6, reltol = 1e-4, stop_at_steady = true,
                          stop_at = nothing, zvt = nothing)
    # Returns an (immutable) `DynNetworkRateParams` with all the static data the ODE needs
    p = _build_rate_params(tstruct, net, infiltration, inflow;
                           runoff=runoff, precipitation=precipitation, zvt=zvt)
    return solveDynNetwork!(state, tstruct, p; tmax=tmax, abstol=abstol, reltol=reltol,
                            stop_at_steady=stop_at_steady, stop_at=stop_at)
end

"""
    solveDynNetwork!(state, tstruct, p::DynNetworkRateParams;
                     tmax=Inf, abstol=1e-6, reltol=1e-4, stop_at_steady=true, stop_at=nothing)
                     -> (; time, trap, kind)

As above, but from prebuilt rate params `p` (see [`_build_rate_params`](@ref)), which already
hold the network, its inflow and everything read from the runoff and precipitation grids.  A
[`DynNetworkContext`](@ref) builds `p` once and solves from it, so the rates stay those of the
context's creation even when the live grids change before a commit.
"""
function solveDynNetwork!(state::AbstractVector{Float64},
                          tstruct::TrapStructure,
                          p::DynNetworkRateParams;
                          tmax = Inf, abstol = 1e-6, reltol = 1e-4, stop_at_steady = true,
                          stop_at = nothing)
    net = p.net
    nt  = length(net.traps)
    @assert length(state) == nt + _nbs_state_count(p) """
        state must have one entry per trap in net.traps plus one per NBS layer \
        ($(nt) traps + $(_nbs_state_count(p)) NBS layer states)"""
    V0 = copy(state)   # immutable snapshot; ODE evolves this, state is updated at the end

    _validate_network(tstruct, net, V0, p.geom)

    # initial rates, for the t=0 fast-path checks
    du0 = similar(V0, Float64)
    dynNetworkRateFunction!(du0, V0, p, 0.0)
    nreg = numregions(tstruct)

    # return an immediate event if the entry state cannot be held (state left for the caller)
    ev0 = _t0_fast_path(V0, du0, p, nreg)
    ev0 === nothing || return ev0

    # where a solve that meets no event ends (elapsed)
    t_stop = stop_at === nothing ? tmax : stop_at
    @assert t_stop <= tmax "stop_at must lie within the window"
    # empty window: nothing to advance
    t_stop <= 0.0 && return (time = max(t_stop, 0.0), trap = 0, kind = :none)

    # evolving = every non-FULL trap.  Empty/zero-inflow traps are included: their inflow can
    # start mid-integration (culverts, NBS) and without a :fill callback V could exceed capacity.
    evolving = [i for i in 1:nt if V0[i] != p.geom[i].capacity]

    # NBS layer states evolve but never fire a topology event, so they join the steady test only
    ss_indices = p.nbsplan === nothing ? evolving :
        vcat(evolving, collect((nt + 1):(nt + p.nbsplan.nlayer_total)))

    # nothing evolving: nothing can change
    isempty(ss_indices) && return (time = stop_at_steady ? Inf : t_stop, trap = 0, kind = :none)
    # already steady (only when predicting; a solve to `tmax` must integrate even slow flow)
    layer_tol = _NBS_SETTLED_RTOL * maximum(i -> abs(du0[i]), ss_indices)
    stop_at_steady &&
        all(i <= nt ? du0[i] == 0.0 : abs(du0[i]) < layer_tol for i in ss_indices) &&
        return (time = Inf, trap = 0, kind = :none)

    # cb_topo: topology events (:fill/:empty/:unspill).  cb_ss: steady state — a DiscreteCallback
    # because step-function infil must be checked at accepted steps, not interpolated ones.  A
    # trap-free NBS-only net has no topology events, so cb_topo is skipped.
    cbs = Any[]
    if nt == 0 # no traps
        event = DynNetworkEvent()
    else
        cb_topo, event = _build_event_callback(p, evolving, V0, nreg)
        push!(cbs, cb_topo)
    end
    stop_at_steady && push!(cbs, p.nbsplan === nothing ?
        _build_steadystate_callback(p, evolving, du0) :
        _build_steadystate_callback_nbs(p, ss_indices, layer_tol, du0))
    # stop at `stop_at`: a step lands on it (tstops) and the solve ends there
    stop_at === nothing ||
        push!(cbs, DiscreteCallback((u, t, integrator) -> t >= stop_at, terminate!;
                                    save_positions = (false, false)))
    callback = CallbackSet(cbs...)
    # finite cap so DiffEq never gets Inf as a tspan endpoint (its internal range(t,Inf,n) rejects it)
    tmax_ode = min(tmax, 1e12)
    # Auto-switching Tsit5/Rodas5P. A NBS layer's overflow k*(S-Smax) is a fast decaying mode
    # (eigenvalue ~ -k), so an explicit method is stability-capped at dt ~ 3.3/k regardless of
    # accuracy — over a long tmax (a last-period solve with no next event marching to a distant
    # trap fill) that means ~1e6 tiny steps. The stiff switch integrates that mode A-stably, so
    # dt is accuracy-bound and grows; Tsit5 still handles the easy segments. Finite-diff Jacobian
    # because the RHS (_nbs_injections) allocates Float64 and isn't Dual-safe for ForwardDiff.
    sol = solve(ODEProblem(dynNetworkRateFunction!, V0, (0.0, tmax_ode), p),
                AutoTsit5(Rodas5P(autodiff = AutoFiniteDiff()));
                callback = callback, tstops = stop_at === nothing ? Float64[] : [stop_at],
                abstol = abstol, reltol = reltol,
                # only the end point is read: don't keep the trajectory
                save_everystep = false, save_start = false, dense = false)

    state .= sol.u[end]        # write ODE result back in place

    # no topology event: tmax cutoff (stopped at tmax_ode) vs genuine steady state (earlier)
    if event.kind == :none
        stop_at === nothing || sol.t[end] < stop_at ||
            return (time = stop_at, trap = 0, kind = :none)
        return sol.t[end] >= tmax_ode ? (time = tmax, trap = 0, kind = :none) :
                                        (time = Inf,  trap = 0, kind = :none)
    end

    # clamp the trigger to its exact threshold so the next call's validator passes
    if event.kind == :fill
        state[event.trap] = p.geom[event.trap].capacity
    elseif event.kind == :empty
        state[event.trap] = 0.0
    end

    return (time = sol.t[end],
            trap = net.traps[event.trap].trap_ix,
            kind = event.kind)
end
