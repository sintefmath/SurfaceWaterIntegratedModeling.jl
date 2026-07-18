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

    An NBS has none.  Its output is `O_terrain(V)` / `E_l(V)`, a function of the **stored layer
    state**, so [`_nbs_routing`](@ref) fixes `path_diff` / `trap_diff` before the walk begins;
    the router only adds them, exactly as it adds `external_inflow`.  In the other direction
    `nbs_draw` is write-only during the walk and is not read until after
    [`_routed_inflow`](@ref) returns.  Neither direction reads what the other writes, so no
    ordering constraint exists to encode — a carrier path may legitimately be visited before
    the inflow path feeding the same placement.

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
culvert inlet draws `min(capacity, available)`; an NBS inlet intercepts everything passing and
zeroes the carrier.  Cell `k`'s residual is charged on the segment *leaving* `k`.

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
        else                                              # :nbsin — intercept ALL passing flow,
            nbs_draw[idx] += current                      # signed: an upstream NBS storing shows up
            current = 0.0                                 # as a deficit in this NBS's live input I_1
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
- `spilling`: whether each trap is at capacity.
- `cvplan`, `trap_level`: passed to [`_culvert_flow`](@ref); `nothing` cvplan skips culverts.
- `culvert_actual`: **mutated** — records what each culvert leaving this trap draws.

# Returns
`nothing`; everything lands in the mutated accumulators.

A full trap emits `max(trap_inflow[i] - footprint_infil[i], 0)`.  If its `spill_path` is 0 it
spills out of the domain, so the surplus simply leaves.
"""
function _route_trap_node!(i, net, trap_inflow, path_flow, footprint_infil,
                           spilling, cvplan, trap_level, culvert_actual)
    if cvplan !== nothing
        trap = net.traps[i]
        for ci in trap.culvert_outlets                # deliver (source drew earlier)
            trap_inflow[i] += culvert_actual[ci]
        end
        for ci in trap.culvert_inlets                 # drain (delivered at its outlet)
            q = _culvert_flow(cvplan, net, ci, trap_level)
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
                cvplan = nothing, trap_level = nothing, nbsrt = nothing) -> Vector{Float64}

The router: how much water arrives at each trap, given who is spilling right now.

All static routing data is pre-supplied — [`_build_rate_params`](@ref) precomputes it once per
solve onto [`DynNetworkRateParams`](@ref) — so the per-step cost here is just the walk.

# Arguments
- `net`: the network being routed.
- `external_inflow`: constant inflow per trap (terrain runoff); seeds the result.
- `spilling`: whether each trap is at capacity, hence emitting.
- `footprint_infil`: whole-footprint loss per trap, deducted from a full trap's spill.
- `path_runoff`: per-path oblivious residual, the attenuation reference.  Read only.
- `path_events`: per-path ordered stops ([`_path_event_templates`](@ref)).
- `order`, `merge_target`: the routing plan ([`_network_order`](@ref), [`_merge_targets`](@ref)).
- `cvplan`: [`CulvertPlan`](@ref), or `nothing` to ignore culverts entirely.
- `trap_level`: per-trap water surface, needed only for culvert heads.
- `nbsrt`: [`NBSRouting`](@ref), or `nothing` with no NBS.  Its `path_diff` / `trap_diff` are
  read in; its `nbs_draw` is **written through** with the intercepted flow.

# Returns
Total inflow at each trap, in `net.traps` order.

Walks `order`, so everything upstream of a node is final when it is visited.  Nothing the caller
owns is mutated apart from `nbsrt`: the accumulators are built fresh here.  Segmented routing
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
                     nbsrt = nothing)

    np = length(net.flow_paths)
    nt = length(net.traps)
    @assert length(external_inflow) == length(spilling) == length(footprint_infil) == nt
    @assert length(path_runoff) == length(merge_target) == length(path_events) == np

    # Working accumulators: `trap_inflow` seeded with external trap inflow (plus any NBS
    # internal-depression diff); `path_flow`, `trib_output` and `culvert_actual` zeroed —
    # every path's water arrives from an upstream trap's spill or an NBS head-injection.
    trap_inflow = Float64.(external_inflow)
    nbsrt === nothing || (trap_inflow .+= nbsrt.trap_diff)
    path_flow   = zeros(Float64, np)
    trib_output = zeros(Float64, np)
    culvert_actual = cvplan === nothing ? Float64[] : zeros(Float64, length(net.culverts))
    nbs_draw = nbsrt === nothing ? Float64[] : nbsrt.nbs_draw

    for node in order
        if node <= np                               # a flow path
            p = node
            # NBS output diff (if any) is head-injected here and rides the same signed tracker.
            head = nbsrt === nothing ? path_flow[p] : path_flow[p] + nbsrt.path_diff[p]
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
                              spilling, cvplan, trap_level, culvert_actual)
        end
    end
    return trap_inflow
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
# NBS as a signed-diff router element (see agent/NBS_ROUTING_NODE_PLAN.md).
#
# `watercourses` is NBS-oblivious: the baseline oblivious runoff already passes through the
# footprint and routes downstream (its output is `O_0`).  The NBS captures its *live* total
# input `I_1` at the footprint-inflow cells (`:nbsin` router draws + the oblivious background
# throughput `O_0_total`), runs its layered storage ODE (state appended after the `nt` trap
# volumes), and emits only the **diff** `O_1 - O_0` of its live output over the oblivious
# baseline.  Net downstream is `O_0 + diff = O_1`; `external_inflow` keeps the `O_0` baseline
# (no double-count).
#
# Output diff placement.  Each terrain endpoint gets `(O_terrain(V) - O_0_total)*ratio_e`,
# `ratio_e = O_0[e]/O_0_total` (even split when `O_0_total ≈ 0`); boundary exits head-inject
# it on the carrier path departing from the exit cell, internal depressions deposit straight
# into the accumulation trap.  Each piped outlet head-injects `+E_l(V)` on the path departing
# from its outlet cell.  A head-injected diff is attenuated along the carrier path's oblivious
# runoff by `_attenuate_diff` and delivered to the path's target trap, cascading onward via
# the normal trap spill.  Fluxes are power-law in layer storage (`compute_outflow`), mm↔m^3
# via `S_mm = V*1000/A`.

"""
    NBSElement

One NBS placement, resolved into router terms for a solve (see the section header above for the
signed-diff model).

# Fields
- `system`: the layered storage model ([`NBSSystem`](@ref))
- `A`: footprint area (m^2), for `S_mm = V*1000/A`
- `state_base`: 0-based offset of this element's layer block, after the `nt` trap states
- `n_terrain`: number of top layers re-emitting at terrain rather than through a pipe
- `O_0_total`: the placement's oblivious throughput — read off the **input** boundary, as
  `Σ runoff` over `footprint_inflow_cells` plus the rain landing on the footprint.  The
  footprint has zero infiltration by contract, so nothing is lost in transit and throughput
  equals what enters; it therefore also seeds the ODE's live input.
- `terrain_paths`: `(carrier path, ratio_e)` per boundary outlet
- `terrain_traps`: `(accumulation trap, ratio_e)` per internal depression
- `piped_paths`: `(carrier path, layer index)` per piped outlet

`ratio_e` is a *normalised share*, weighted by the oblivious `runoff` at each outlet cell (and
at each ponding cell) — the receiving watercourse, not the footprint cells feeding it.  It is
independent of `O_0_total`: several footprint cells may feed one outlet, and the outlet's own
rain and any flow converging on it from outside the footprint are counted in its weight but
must never be counted in the throughput.
"""
struct NBSElement
    system    ::NBSSystem
    A         ::Float64                     # footprint area (m^2) for S_mm = V*1000/A
    state_base::Int                         # 0-based offset of its layer block after the nt trap states
    n_terrain ::Int                         # top layers re-emitting at terrain
    O_0_total ::Float64                     # oblivious throughput, from the input boundary:
                                            # Σ runoff over footprint_inflow_cells + rain on footprint
                                            # (exact: zero footprint infiltration => output == input)
    terrain_paths::Vector{Tuple{Int,Float64}}  # (carrier path, ratio_e) per boundary outlet
    terrain_traps::Vector{Tuple{Int,Float64}}  # (accumulation trap, ratio_e) per internal depression
    piped_paths  ::Vector{Tuple{Int,Int}}      # (carrier path, layer index) per piped outlet
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

# Below this, the endpoint weights count as zero and the terrain diff is split evenly across
# endpoints rather than by ratio (which would be 0/0).
const _NBS_O0_EPS = 1e-12

"""
    _footprint_rain(precipitation, footprint) -> Float64

Total rain landing on a footprint.

# Arguments
- `precipitation`: the period's rain rate — a uniform scalar, or a per-cell grid.
- `footprint`: the placement's cells (linear indices).

# Returns
The rain summed over `footprint`, in the same units as the `runoff` grid: `watercourses` seeds
runoff as `precipitation - infiltration` per cell and only ever sums it downstream, applying no
area or mm conversion, so the two are directly comparable.  `O(1)` for a uniform rate.

Infiltration is not subtracted — an NBS footprint has none by contract (`fill_sequence` forces
it to zero), so every drop of rain on it becomes flow.
"""
_footprint_rain(precipitation::Real, footprint::Vector{Int}) =
    Float64(precipitation) * length(footprint)
_footprint_rain(precipitation::AbstractMatrix{<:Real}, footprint::Vector{Int}) =
    sum(Float64[Float64(precipitation[f]) for f in footprint]; init = 0.0)

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
    NBSRouting

Per-step NBS scratch, rebuilt each rate-function call and threaded through the router.  It
carries traffic in **both** directions, and each field is indexed in a **different space**.

# Fields

| Field | Indexed by | Direction |
|---|---|---|
| `path_diff` | flow path (`1:np`) | in — router reads |
| `trap_diff` | trap (`1:nt`) | in — router reads |
| `nbs_draw` | NBS placement (`nbsplan.elems`) | out — router writes |

**In.** `path_diff` and `trap_diff` are the placements' signed output diffs, computed from the
current layer state by [`_nbs_routing`](@ref).  Both are the same quantity — a share
`(O_terrain - O_0_total) * ratio_e` — differing only in where it lands: `path_diff` is injected
at a carrier path's head and then attenuated along it (and also collects each piped outlet's
`+E_l`), while `trap_diff` is deposited straight into an internal depression's trap.

**Out.** `nbs_draw` is the *network-routed* flow intercepted at each placement's
footprint-inflow cells (the `:nbsin` events).  Layer 1 is fed by the live input
`I_1 = O_0_total + nbs_draw`: the oblivious half is static and known when the plan is built,
but the network half depends on which traps are spilling right now, so it can only be collected
while routing — hence a router output rather than a plan field.  It is signed: an upstream NBS
that is storing water rides a negative diff, which arrives here as a *deficit* in this
placement's input rather than as flow that quietly bypassed it.

**The two directions do not meet within a step, and must not.** The `in` fields come from the
layer state `V`; the `out` field feeds `dV`.  So a placement's output does *not* reflect the
draw collected alongside it — outflow follows stored volume, not what is arriving this instant,
and the draw reaches the output only through integration.  That is what keeps the RHS free of an
algebraic loop, and why routing needs no NBS ordering (see [`_network_order`](@ref)).  Sanity
check: at steady state `dS/dt = 0`, so `O_terrain = I_1` and `diffbase = nbs_draw` — the
placement passes the network draw straight through.

The attenuation reference `path_runoff` is not here — it is static and lives on the rate params.
"""
struct NBSRouting
    path_diff  ::Vector{Float64}   # per flow path — signed diff injected at its head (in)
    trap_diff  ::Vector{Float64}   # per trap — signed diff deposited into it (in)
    nbs_draw   ::Vector{Float64}   # per NBS placement — :nbsin flow captured while routing (out)
end

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

Resolve every NBS placement in `net` into router terms, once per solve.

# Arguments
- `net`: supplies the placements (`net.nbs`), and the paths/traps their outputs target.
- `tstruct`: supplies trap footprints, for mapping a ponding cell to its covering trap.
- `runoff`: the oblivious runoff grid.  Read only.
- `precipitation`: the period's rain rate (scalar or per-cell).

# Returns
The [`NBSPlan`](@ref), or `nothing` when `net` has no NBS.

Per placement this resolves two independent quantities, from `runoff`, `precipitation`, and the
cell lists `setup_network` already cached on the [`DynNBSPlacement`](@ref):

- **`O_0_total`**, the throughput, off the *input* boundary — `Σ runoff` over
  `footprint_inflow_cells` plus the rain on the footprint.  Exact, because the footprint has
  zero infiltration by contract, so output equals input.
- **`ratio_e`**, a normalised share per endpoint, weighted by the `runoff` at each outlet cell
  in `footprint_outflow_cells` and each cell in `internal_accumulation_cells`.

Keeping them separate matters: an outlet's weight legitimately includes its own rain and any
flow converging on it from outside the footprint, none of which is the placement's throughput.
The layer models and ODE-state layout come straight off `net.nbs`.

Nothing is mutated.  Asserts the coupling invariants: an outlet must have a carrier path
departing it (it is a network seed), and an internal depression must be covered by a network
trap.
"""
function _build_nbs_plan(net::DynNetwork, tstruct, runoff::AbstractMatrix{<:Real},
                         precipitation::Union{AbstractMatrix{<:Real},Real})
    isempty(net.nbs) && return nothing
    LI = LinearIndices(tstruct.topography)
    trap_of_cell = Dict{Int,Int}()             # linear cell -> local trap index
    for (ti, t) in enumerate(net.traps), c in tstruct.footprints[t.trap_ix]
        trap_of_cell[c] = ti
    end
    # carrier path per seed cell: the path departing from an output (outlet) cell
    dep_path = Dict{CartesianIndex{2},Int}()
    for (p, fp) in enumerate(net.flow_paths)
        haskey(dep_path, fp.departure_point) || (dep_path[fp.departure_point] = p)
    end

    elems = NBSElement[]
    base  = 0
    for nb in net.nbs
        A_foot = Float64(length(nb.footprint))   # @@@ 1 m^2/cell; use real cell area when available

        # Throughput, off the INPUT boundary.  The footprint has zero infiltration by contract,
        # so nothing is lost crossing it and output == input == (flow in) + (rain on footprint).
        # The flowgraph is single-successor, so an inflow cell's whole runoff enters; a cell with
        # negative runoff (spare infiltration capacity) never propagates, hence max(., 0).
        O_0_total = _footprint_rain(precipitation, nb.footprint)
        for ic in nb.footprint_inflow_cells
            O_0_total += max(Float64(runoff[ic]), 0.0)
        end

        # Weights: the oblivious flow at each *outlet* cell (outside the footprint) and at each
        # ponding cell.  Only ever used as a normalised share, so the outlet's own rain and any
        # flow converging on it from outside are harmless here (they must not, and do not, reach
        # O_0_total above).
        weights = Tuple{Float64,Symbol,Int}[]
        for oc in nb.footprint_outflow_cells
            carrier = get(dep_path, oc, 0)       # the outlet cell is a seed, so a path must depart it
            @assert carrier != 0 "NBS outlet cell $oc has no carrier path " *
                "(expected it to be a network seed)"
            push!(weights, (max(Float64(runoff[oc]), 0.0), :path, carrier))
        end
        for pc in nb.internal_accumulation_cells
            tr = get(trap_of_cell, LI[pc], 0)
            @assert tr != 0 "NBS internal-depression cell $pc (region $(tstruct.regions[LI[pc]])) " *
                "is covered by no network trap — coupling invariant broken"
            push!(weights, (max(Float64(runoff[pc]), 0.0), :trap, tr))
        end
        W_total = sum(Float64[w for (w, _, _) in weights]; init = 0.0)
        nend    = length(weights)
        ratio(w) = W_total > _NBS_O0_EPS ? w / W_total : (nend > 0 ? 1.0 / nend : 0.0)

        terrain_paths = Tuple{Int,Float64}[]
        terrain_traps = Tuple{Int,Float64}[]
        for (w, kind, tgt) in weights
            kind === :path ? push!(terrain_paths, (tgt, ratio(w))) :
                             push!(terrain_traps, (tgt, ratio(w)))
        end

        piped_paths = Tuple{Int,Int}[]
        piped = 0
        for (l, L) in enumerate(nb.system.layers)
            (l > nb.n_terrain && L.Kout > 0.0) || continue
            piped += 1
            carrier = get(dep_path, nb.outlets[piped], 0)   # the outlet cell is a seed too
            @assert carrier != 0 "NBS piped outlet $(nb.outlets[piped]) has no carrier path " *
                "(expected it to be a network seed)"
            push!(piped_paths, (carrier, l))
        end

        push!(elems, NBSElement(nb.system, A_foot, base, nb.n_terrain, O_0_total,
                                terrain_paths, terrain_traps, piped_paths))
        base += length(nb.system.layers)
    end
    return NBSPlan(elems, base)
end

"""
    _nbs_routing(V, p, nt, np) -> NBSRouting

The signed output diffs of every NBS placement at the current layer state.

# Arguments
- `V`: the ODE state — trap volumes, then the NBS layer states.
- `p`: the [`DynNetworkRateParams`](@ref); `p.nbsplan` must not be `nothing`.
- `nt`, `np`: trap and flow-path counts, sizing the two output vectors.

# Returns
A fresh [`NBSRouting`](@ref) with `path_diff` and `trap_diff` filled and `nbs_draw` zeroed for
the router to fill on the way back.  Nothing is mutated; the scratch is allocated per call.

Terrain re-emission contributes `(O_terrain(V) - O_0_total) * ratio_e` to `path_diff` (injected
at the carrier path's head) or to `trap_diff` (deposited into the accumulation trap) — the same
share either way, differing only in where it lands.  Each piped outlet adds `+E_l(V)` to its
carrier path.

Emitting the *diff* over the oblivious baseline, rather than the output itself, is what keeps
`external_inflow`'s `O_0` from being counted twice.
"""
function _nbs_routing(V, p, nt::Int, np::Int)
    path_diff  = zeros(Float64, np)
    trap_diff  = zeros(Float64, nt)
    for el in p.nbsplan.elems
        base = el.state_base
        O_terrain = 0.0
        @inbounds for l in 1:el.n_terrain
            L = el.system.layers[l]
            O_terrain += compute_outflow(L.Kout, L.nout, L.Smax,
                                         V[nt + base + l] * 1000.0 / el.A) * el.A * 1e-3
        end
        diffbase = O_terrain - el.O_0_total
        for (pth, r) in el.terrain_paths; path_diff[pth]  += diffbase * r; end
        for (tr,  r) in el.terrain_traps; trap_diff[tr]   += diffbase * r; end
        for (pth, l) in el.piped_paths
            L = el.system.layers[l]
            E = compute_outflow(L.Kout, L.nout, L.Smax,
                                V[nt + base + l] * 1000.0 / el.A) * el.A * 1e-3
            path_diff[pth] += E
        end
    end
    return NBSRouting(path_diff, trap_diff, zeros(Float64, length(p.nbsplan.elems)))
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
    nt = length(net.traps)
    np = length(net.flow_paths)
    @assert length(external_inflow) == nt
    order, _        = _network_order(net)
    # Per-path oblivious residual runoff — the reference the router attenuates all flow against.
    # With a real grid it is background-aware (a cell the background saturated charges the network
    # flow nothing); without one (hand-built nets) cells sit at their capacity floor (-infiltration),
    # which reproduces the plain infiltration-prefix behaviour.
    path_runoff = runoff === nothing ?
        [Float64[-c for c in ci] for ci in _path_cell_values(net, infiltration)] :
        _path_cell_values(net, runoff)
    cvplan  = isempty(net.culverts) ? nothing : _build_culvert_plan(net, tstruct)
    if isempty(net.nbs)
        nbsplan = nothing
    else
        runoff === nothing &&
            error("_build_rate_params: net has NBS but no runoff grid was supplied")
        precipitation === nothing &&
            error("_build_rate_params: net has NBS but no precipitation was supplied " *
                  "(needed for each placement's throughput: rain on the footprint + inflow)")
        nbsplan = _build_nbs_plan(net, tstruct, runoff, precipitation)
    end
    # in-order tributary/culvert stops per path (tributaries only when there are no culverts)
    events = _path_event_templates(net)
    tgeom = _build_trap_geometry(tstruct, net, infiltration; zvt=zvt)
    footprint_infil = _footprint_infiltration(tgeom)
    return DynNetworkRateParams(net,
                                tgeom,
                                Float64.(external_inflow),
                                footprint_infil,
                                path_runoff,
                                order,
                                _merge_targets(net),
                                cvplan,
                                events,
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
`(inflow, spilling, nbs_draw)` — total inflow per trap, whether each trap is at capacity (and so
spilling), and the NBS live-input draws (`nothing` with no NBS).

Shared by the rate function and the `:unspill` condition, so the two cannot disagree about who
is spilling.  With NBS the signed output diffs are folded into the routing via `nbsrt`, and the
`:nbsin` flow captured during routing comes back as `nbs_draw`.  Nothing the caller owns is
mutated; scratch is allocated per call.

!!! note "`nbs_draw` is complete on return"
    A placement's draw is summed across every path crossing its footprint-inflow cells, so it is
    only trustworthy once *all* of them have been walked.  Three things make that hold, and a
    caller must not read `nbs_draw` before this function returns:

    1. `_nbs_coupled_nodes` (`network_utils.jl`) forces every path crossing a footprint-inflow
       cell into that placement's component, so no contributor can sit in another `net`.
    2. [`_route_flow`](@ref) walks every node of `order`, so every `:nbsin` fires.
    3. `nbs_draw` is write-only during the walk — order of contributions is irrelevant to a sum.

    Point 1 is the load-bearing one: were the coupling incomplete, a path's flow would simply
    never reach the placement.  It cannot go unnoticed, though — `_localize_path` resolves each
    `nbs_inlet` through `nbsmap`, so a path landing in a component without its NBS throws at
    build rather than under-counting silently.
"""
function _routed_inflow(V, p::DynNetworkRateParams)
    geom = p.geom
    nt   = length(geom)
    spilling = Bool[V[i] >= geom[i].capacity for i in 1:nt]
    # culverts need each trap's water-surface elevation for the head calc
    trap_level = p.cvplan === nothing ? nothing :
                 Float64[water_level(geom[i], V[i]) for i in 1:nt]
    nbsrt = p.nbsplan === nothing ? nothing :
            _nbs_routing(V, p, nt, length(p.net.flow_paths))
    inflow = _route_flow(p.net, p.external_inflow, spilling,
                         p.footprint_infil, p.path_runoff,
                         p.path_events, p.order, p.merge_target,
                         p.cvplan, trap_level, nbsrt)
    return inflow, spilling, (nbsrt === nothing ? nothing : nbsrt.nbs_draw)
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

    # NBS layer-state derivatives (m^3/time).  Layer 1 (top) is fed by the live total input
    # `I_1` (oblivious background throughput `O_0_total` + the :nbsin network draws captured
    # during routing); each lower layer by the layer above's infiltration; the bottom layer's
    # infiltration leaves the system to ground.  Each layer loses its own overflow and its
    # infiltration.  Fluxes are power-law in the layer storage, in mm, converted to m^3
    # (`S_mm = V*1000/A`, `Q_m3 = Q_mm*1e-3`).
    if p.nbsplan !== nothing
        @inbounds for (k, el) in enumerate(p.nbsplan.elems)
            base    = el.state_base
            # layer-1 inflow = live input I_1 = O_0_total + the :nbsin draws.  Physically >= 0
            # (a negative draw is an upstream NBS's deficit, capped by attenuation at the flow
            # present, which is part of this footprint's own O_0_total); max() guards float noise.
            prev_qi = max(el.O_0_total + nbs_draw[k], 0.0)
            for (l, L) in enumerate(el.system.layers)
                S_mm = V[nt + base + l] * 1000.0 / el.A
                # depth-rate (mm/time) -> volume-rate (m^3/time): * A / 1000
                qo   = compute_outflow(L.Kout, L.nout, L.Smax, S_mm) * el.A * 1e-3
                qi   = compute_outflow(L.Kinf, L.ninf, L.Smin, S_mm) * el.A * 1e-3
                # @@@ evapotranspiration deferred: ET is an explicit 0.0 placeholder, so
                #     wiring it in (from EVCoeff/EVS11) is a one-line change per layer.
                ET = 0.0
                dV[nt + base + l] = prev_qi - qo - qi - ET
                # Physical floor: a layer cannot drain below empty.
                V[nt + base + l] <= 0.0 && dV[nt + base + l] < 0.0 && (dV[nt + base + l] = 0.0)
                prev_qi = qi
            end
        end
    end
    return nothing
end

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
# 2. cb_ss (DiscreteCallback): steady state — max|dV/dt| < abstol at an accepted step.  Discrete,
#    not Continuous: wetted_infiltration is a step function of V, so interpolated mid-step states
#    would give spurious crossings.  Vetoes when any evolving trap is at V >= C (:fill wins).
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
    _build_steadystate_callback(p, evolving, abstol, du0) -> DiscreteCallback

Steady-state detector for an NBS-free network.

# Arguments
- `p`: the [`DynNetworkRateParams`](@ref) for this solve.
- `evolving`: local indices of the traps to watch.
- `abstol`: the `|dV/dt|` threshold below which a trap counts as settled.
- `du0`: each trap's rate at `t = 0`, for the sign-change test.

# Returns
A `DiscreteCallback` that terminates at the first accepted step where every evolving trap has
settled — `|dV/dt| < abstol`, or its rate has crossed the sign it started at, meaning a trap
oscillating across a wetted-infiltration step has reached a sub-capacity equilibrium (the
multi-trap analogue of `fill_trap_until`'s stagnation stop).

Discrete rather than continuous because `wetted_infiltration` is a step function of `V`, so
interpolated mid-step states would give spurious crossings.  Vetoes while any evolving trap sits
at `V >= C`, letting `:fill` terminate instead.
"""
function _build_steadystate_callback(p::DynNetworkRateParams,
                                     evolving::AbstractVector{<:Integer},
                                     abstol::Real,
                                     du0::AbstractVector{<:Real})
    function condition(u, t, integrator)
        inflow, spilling = _routed_inflow(u, p)
        for e in evolving
            spilling[e] && return false   # veto: :fill wins
        end
        for e in evolving
            rate = inflow[e] - wetted_infiltration(p.geom[e], u[e])
            # Match the rate function's physical floor: an empty trap cannot drain
            # further, so its actual dV is 0 (otherwise cb_ss never sees a leaf that has
            # drained to V == 0 settle — its unclamped rate stays negative and the solve
            # marches to t → ∞).
            u[e] <= 0.0 && rate < 0.0 && (rate = 0.0)
            # Settled if the actual rate is ~0, OR it has crossed its starting sign — a
            # trap oscillating across a wetted-infiltration step has reached a
            # sub-capacity equilibrium and makes no further net progress (the multi-trap
            # analogue of `fill_trap_until`'s stagnation/sign-change stop).
            settled = abs(rate) < abstol ||
                      (du0[e] != 0.0 && signbit(rate) != signbit(du0[e]))
            settled || return false
        end
        return true
    end
    return DiscreteCallback(condition, terminate!; save_positions = (false, false))
end

"""
    _build_steadystate_callback_nbs(p, ss_indices, abstol, du0) -> DiscreteCallback

Steady-state detector for a network with NBS.

# Arguments
- `p`: the [`DynNetworkRateParams`](@ref) for this solve.
- `ss_indices`: state indices to watch — the evolving traps plus the NBS layer states.
- `abstol`: the `|dV/dt|` threshold below which an index counts as settled.
- `du0`: each state's rate at `t = 0`, for the sign-change test.

# Returns
A `DiscreteCallback`, as [`_build_steadystate_callback`](@ref) but evaluating the *full* rate
function: the cheap `_routed_inflow` check cannot express the layer dynamics.  A trap at
capacity still vetoes, letting `:fill` terminate.
"""
function _build_steadystate_callback_nbs(p::DynNetworkRateParams,
                                         ss_indices::AbstractVector{<:Integer},
                                         abstol::Real,
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
            rate    = dbuf[i]
            settled = abs(rate) < abstol ||
                      (du0[i] != 0.0 && signbit(rate) != signbit(du0[i]))
            settled || return false
        end
        return true
    end
    return DiscreteCallback(condition, terminate!; save_positions = (false, false))
end

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
                     abstol=1e-6, reltol=1e-4, zvt=nothing) -> (; time, trap, kind)

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
- `abstol`, `reltol`: ODE tolerances; `abstol` also sets the steady-state |dV/dt| threshold

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
                          abstol = 1e-6, reltol = 1e-4,
                          zvt = nothing)
    nt = length(net.traps)
    # Returns an (immutable) `DynNetworkRateParams` with all the static data the ODE needs
    p  = _build_rate_params(tstruct, net, infiltration, inflow;
                            runoff=runoff, precipitation=precipitation, zvt=zvt)
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

    # empty window: nothing to advance
    tmax <= 0.0 && return (time = max(tmax, 0.0), trap = 0, kind = :none)

    # evolving = every non-FULL trap.  Empty/zero-inflow traps are included: their inflow can
    # start mid-integration (culverts, NBS) and without a :fill callback V could exceed capacity.
    evolving = [i for i in 1:nt if V0[i] != p.geom[i].capacity]

    # NBS layer states evolve but never fire a topology event, so they join the steady test only
    ss_indices = p.nbsplan === nothing ? evolving :
        vcat(evolving, collect((nt + 1):(nt + p.nbsplan.nlayer_total)))

    # already steady: nothing evolving, or every evolving rate ~0
    isempty(ss_indices) && return (time = Inf, trap = 0, kind = :none)
    all(abs(du0[i]) <= abstol for i in ss_indices) &&
        return (time = Inf, trap = 0, kind = :none)

    # cb_topo: topology events (:fill/:empty/:unspill).  cb_ss: steady state — a DiscreteCallback
    # because step-function infil must be checked at accepted steps, not interpolated ones.  A
    # trap-free NBS-only net has no topology events, so cb_topo is skipped.
    cb_ss = p.nbsplan === nothing ?
        _build_steadystate_callback(p, evolving, abstol, du0) :
        _build_steadystate_callback_nbs(p, ss_indices, abstol, du0)
    if nt == 0 # no traps
        event    = DynNetworkEvent()
        callback = cb_ss
    else
        cb_topo, event = _build_event_callback(p, evolving, V0, nreg)
        callback = CallbackSet(cb_topo, cb_ss)
    end
    # finite cap so DiffEq never gets Inf as a tspan endpoint (its internal range(t,Inf,n) rejects it)
    tmax_ode = min(tmax, 1e12)
    # Explicit Tsit5: the RHS has C0 kinks (culvert regime switches, routing `min`) but isn't
    # stiff.  The auto-switcher misreads the kinks as stiffness and pays for dense Jacobians
    # (~4x slower, no accuracy gain).
    # @@@ revisit with an auto-switching solver + sparse Jacobian if a genuinely stiff config appears
    sol = solve(ODEProblem(dynNetworkRateFunction!, V0, (0.0, tmax_ode), p), Tsit5();
                callback = callback,
                abstol = abstol, reltol = reltol)

    state .= sol.u[end]        # write ODE result back in place

    # no topology event: tmax cutoff (stopped at tmax_ode) vs genuine steady state (earlier)
    if event.kind == :none
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
