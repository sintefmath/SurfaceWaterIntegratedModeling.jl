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

A [`TrapGeometry`](@ref) for every trap in `net`, in `net.traps` order.  `infiltration` is a
per-cell rate grid (0 = impermeable).  `zvt` is the cached volume↔level table
([`_compute_z_vol_tables`](@ref)); computed here if `nothing`.
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

Water surface elevation of trap `g` at stored volume `V`.  Empty (`V <= 0`) → `zmin`;
full (`V >= capacity`) → `Inf` (whole footprint reads as submerged); else interpolated.
"""
function water_level(g::TrapGeometry, V::Real)
    return V <= 0.0       ? g.zmin :
           V >= g.capacity ? Inf    :
                            Float64(g.v2z(V))
end

# Actual water-surface elevation of trap `g` at volume `V`, with no `Inf` sentinel:
# a full trap sits at its spill level (the table top), not `Inf`.  Used for culvert
# heads, where the real surface elevation matters and `Inf - Inf` would be `NaN`.
function _surface_level(g::TrapGeometry, V::Real)
    return V <= 0.0 ? g.zmin : Float64(g.v2z(min(V, g.capacity)))
end

# ----------------------------------------------------------------------------
"""
    wetted_infiltration(g::TrapGeometry, V) -> Float64

Infiltration loss of trap `g` at volume `V`: the per-cell rate summed over the
currently-submerged footprint cells (bottom at or below the water level).  A full trap
(level `Inf`) infiltrates over its whole footprint.
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

# Combined topological order over path nodes (1:np) and trap nodes (np+1:np+nt).  Edges follow
# the flow (trap→spill path, path→target trap, tributary→joined path).  Returns (order, np).
# `net.traps` order alone doesn't fix the path/trap interleaving once tributaries exist, so the
# combined graph is sorted.
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
# `merge_target[p]` = the path that tributary `p` feeds into, or 0 if `p` is not a
# tributary.  Each path is a tributary of at most one other, so this is well defined.
function _merge_targets(net::DynNetwork)
    merge_target = zeros(Int, length(net.flow_paths))
    for (a, path) in enumerate(net.flow_paths), (m, _) in path.merges
        merge_target[m] = a
    end
    return merge_target
end

# ----------------------------------------------------------------------------
# Culvert routing data.  Each culvert's inlet/outlet is owned by a flow path or a trap
# (resolved at construction).  `CulvertPlan` precomputes, per culvert, which kind owns each
# endpoint and its local index, the barrel diameter, and the two invert elevations.

struct CulvertPlan
    diam::Vector{Float64}         # barrel diameter (2r) per culvert
    inlet_invert::Vector{Float64} # terrain elevation at each culvert's inlet cell
    outlet_invert::Vector{Float64}# terrain elevation at each culvert's outlet cell
    inlet_is_path::BitVector
    inlet_owner::Vector{Int}      # local path or trap index owning the inlet
    outlet_is_path::BitVector
    outlet_owner::Vector{Int}     # local path or trap index owning the outlet
end

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

# Add an inlet-owner -> outlet-owner edge per culvert to the routing-order graph, so
# the inlet is routed before the outlet.  Construction already verified the network
# is acyclic with these edges (see `_topological_order`).
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

# Per-path event stream for the culvert-aware router: tributary junctions and culvert
# inlet/outlet positions, sorted by cell position, as `(position, :trib|:cvin|:cvout, idx)`.
# Static for a solve, so built once and reused every rate-function call.
function _path_event_templates(net::DynNetwork)
    return [begin
                ev = Tuple{Int,Symbol,Int}[]
                for (m, j)   in fp.merges;          push!(ev, (j,   :trib,  m));  end
                for (ci, pos) in fp.culvert_inlets;  push!(ev, (pos, :cvin,  ci)); end
                for (ci, pos) in fp.culvert_outlets; push!(ev, (pos, :cvout, ci)); end
                sort!(ev; by = first)
                ev
            end
            for fp in net.flow_paths]
end

# Volumetric flow through culvert `ci` given the current trap water levels.  A flow-
# path endpoint is treated as not-submerged with head fixed at the diameter (so its
# inlet-control capacity is the weir capacity at head D); a trap endpoint uses its
# live water level above the culvert's invert.  Downhill-only (no reverse flow).
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
    _route_flow(net, external_inflow, spilling, footprint_infil, path_cell_infil;
                path_inflow=zeros(np), cvplan=nothing, trap_level=nothing) -> Vector{Float64}

Total inflow arriving at each trap of `net` (in `net.traps` order): `external_inflow`
(per trap) and `path_inflow` (per path) plus everything routed in from upstream spills.

A spilling trap emits `max(inflow - footprint_infil, 0)` into its spill path.  Path flow is
routed in segments between tributary junctions so each tributary is charged only the
infiltration of the cells it actually travels (see the core method below).

This convenience form computes the static routing data on the fly (tests / hand-built nets);
the solver uses the precomputed core method via [`DynNetworkRateParams`](@ref).
"""
function _route_flow(net::DynNetwork,
                     external_inflow::AbstractVector{<:Real},
                     spilling::AbstractVector{Bool},
                     footprint_infil::AbstractVector{<:Real},
                     path_cell_infil::AbstractVector{<:AbstractVector{<:Real}};
                     path_inflow = zeros(length(net.flow_paths)),
                     cvplan = nothing, trap_level = nothing)
    order, _      = _network_order(net)
    prefix        = [_infil_prefix(ci) for ci in path_cell_infil]
    path_events   = _path_event_templates(net)
    return _route_flow(net, external_inflow, spilling, footprint_infil,
                       prefix, path_inflow, path_events, order, _merge_targets(net),
                       cvplan, trap_level)
end

# Flow delivered at the end of a path: head flow travels the cells losing per-segment
# infiltration; at each `events` stop a tributary adds its output, a culvert draws/adds at its
# cell.  `events` is the in-order `(pos, kind, idx)` stream.  Mutates `culvert_actual` with what
# each culvert inlet here drew.
function _path_delivered!(prefix, head_flow, events,
                          trib_output, culvert_actual, cvplan, net, trap_level)
    current  = head_flow
    prev_pfx = 0.0
    for (pos, kind, idx) in events
        current = max(current - (prefix[pos] - prev_pfx), 0.0)
        if kind === :trib
            current += trib_output[idx]
        elseif kind === :cvout
            current += culvert_actual[idx]                # deliver
        else                                              # :cvin
            a = min(_culvert_flow(cvplan, net, idx, trap_level), current)
            culvert_actual[idx] = a                       # drawn == delivered
            current -= a
        end
        prev_pfx = prefix[pos]
    end
    return max(current - (prefix[end] - prev_pfx), 0.0)
end

# Trap node: apply the culvert flows touching this trap (deliver what culverts
# routed in drew; drain what culverts draining it carry), then, if it is full,
# spill its surplus into its spill path.  Mutates `trap_inflow`, `path_flow`, and
# `culvert_actual`.
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

# Core router, all static routing data pre-supplied.  Segmented routing charges flow
# infiltration only over the cells it travels, so a tributary joining at junction j isn't
# charged the main path's pre-junction infiltration.  No tributaries → one lump deduction.
function _route_flow(net::DynNetwork,
                     external_inflow::AbstractVector{<:Real},
                     spilling::AbstractVector{Bool},
                     footprint_infil::AbstractVector{<:Real},
                     path_infil_prefix::AbstractVector{<:AbstractVector{<:Real}},
                     external_path_inflow::AbstractVector{<:Real},
                     path_events::AbstractVector{<:AbstractVector},
                     order::AbstractVector{<:Integer},
                     merge_target::AbstractVector{<:Integer},
                     cvplan = nothing,
                     trap_level = nothing)

    np = length(net.flow_paths)
    nt = length(net.traps)
    @assert length(external_inflow) == length(spilling) == length(footprint_infil) == nt
    @assert length(path_infil_prefix) == length(merge_target) ==
            length(external_path_inflow) == length(path_events) == np

    # Working accumulators: `trap_inflow` seeded with external trap inflow, `path_flow` with
    # per-path head inflow, `trib_output`/`culvert_actual` zeroed.
    trap_inflow = Float64.(external_inflow)
    path_flow   = Float64.(external_path_inflow)
    trib_output = zeros(Float64, np)
    culvert_actual = cvplan === nothing ? Float64[] : zeros(Float64, length(net.culverts))

    for node in order
        if node <= np                               # a flow path
            p = node
            delivered = _path_delivered!(path_infil_prefix[p], path_flow[p],
                                         path_events[p], trib_output, culvert_actual,
                                         cvplan, net, trap_level)
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

Whole-footprint infiltration of each trap (the loss a full, submerged trap incurs): the
per-cell `infil` summed over the footprint.  Used by [`_route_flow`](@ref).
"""
function _footprint_infiltration(geom::Vector{TrapGeometry})
    return [sum(g.infil) for g in geom]
end

# ----------------------------------------------------------------------------
"""
    _path_cell_infiltration(net, infiltration) -> Vector{Vector{Float64}}

Per-cell infiltration of each flow path (in `net.flow_paths` order), the grid sampled along
each path's cells.  [`_build_rate_params`](@ref) turns these into prefix sums for the router.
"""
function _path_cell_infiltration(net::DynNetwork, infiltration::AbstractMatrix{<:Real})
    return [isempty(p.cells) ? Float64[] : Float64[infiltration[c] for c in p.cells]
            for p in net.flow_paths]
end

# Prefix sum of a per-cell infiltration vector: prefix[k] = sum of cells 1..k-1,
# so prefix[1]=0 and prefix[end]=total.  Used for O(1) segment-infil lookups.
function _infil_prefix(cell_infil::AbstractVector{<:Real})
    n      = length(cell_infil)
    prefix = Vector{Float64}(undef, n + 1)
    prefix[1] = 0.0
    for k in 1:n
        prefix[k+1] = prefix[k] + cell_infil[k]
    end
    return prefix
end

# ============================================================================
# Network rate function.
#
# `dynNetworkRateFunction!` is the in-place ODE RHS over the trap-volume state.  Its
# `DynNetworkRateParams` bundle everything static for a solve (geometry, external inflow,
# infiltration sums, routing plan), so per-step work is just routing + a loop over traps.
# ============================================================================

"""
    DynNetworkRateParams

Static precomputed inputs to [`dynNetworkRateFunction!`](@ref) for one solve, built by
[`_build_rate_params`](@ref) and passed as the `ODEProblem` parameter.  Indexed
network-locally (`net.traps` / `net.flow_paths` order).

# Fields
- `net`: the [`DynNetwork`](@ref) being solved
- `geom`: per-trap [`TrapGeometry`](@ref)
- `external_inflow`: constant inflow per trap
- `external_path_inflow`: constant inflow onto each flow path (e.g. rain on path cells)
- `footprint_infil`: whole-footprint infiltration per trap
- `path_infil_prefix`: per-path prefix sums of per-cell infiltration (`prefix[1]=0`,
  `prefix[end]`=total)
- `order`, `merge_target`: the static routing plan ([`_network_order`](@ref), [`_merge_targets`](@ref))
- `path_events`: per-path in-order `(cell_pos, kind, idx)` stops — tributary junctions, plus
  culvert inlet/outlet positions when the net has culverts
- `cvplan`: culvert routing data (`nothing` if no culverts)
- `nbsplan`: NBS signed-correction data (`nothing` if no NBS)
"""

# ============================================================================
# NBS signed-correction model (see agent/NBS_ROUTING_REDESIGN.md).
#
# `watercourses` is NBS-oblivious, so the baseline inflow already routes the full footprint
# pass-through.  An NBS is a *correction* to that baseline.  Its layered storage lives in the
# ODE state after the `nt` trap volumes: the top layer is fed by the oblivious capture `I`,
# each lower layer by the one above's infiltration, the bottom layer's infiltration to ground.
#
# Terrain re-emit correction at each footprint drainage endpoint `c`: `-X*Q_c`, with
# `X = 1 - O_terrain/I` (`O_terrain` = top-`n_terrain` overflow, `Q_c` = oblivious drainage
# there).  Propagated to the endpoint's first downstream trap (`_propagate_correction`) and
# added to its inflow.  Piped outlets inject their discharge as a positive correction.  Fluxes
# are power-law in layer storage (`compute_outflow`), mm↔m^3 via `S_mm = V*1000/A`.
#
# @@@ Deferred: NBS→NBS series coupling (upstream correction lowering a downstream capture).
struct NBSLayerParams
    Kout::Float64; nout::Float64; Smax_mm::Float64
    Kinf::Float64; ninf::Float64; Smin_mm::Float64
    A::Float64                            # layer area (== footprint area, m^2)
end

# Where a signed correction lands: the network-local trap it is added to (0 = it leaves the
# network and is dropped), and the oblivious runoff `V` along the flow-path cells it travels
# through to get there (empty for a correction deposited straight into a trap).
struct NBSDelivery
    trap  ::Int
    path_V::Vector{Float64}
end

# One NBS element's static routing data.
struct NBSElement
    layers    ::Vector{NBSLayerParams}
    state_base::Int                            # 0-based offset of its layer block after the nt trap states
    n_terrain ::Int                            # top layers re-emitting at terrain
    capture   ::Float64                        # I: oblivious footprint capture (Σ Q over endpoints)
    endpoints ::Vector{Tuple{Float64,NBSDelivery}}   # (Q_c, delivery) per terrain drainage endpoint
    outlets   ::Vector{Tuple{Int,NBSDelivery}}       # (piped layer index, delivery) per piped outlet
end

struct NBSPlan
    elems       ::Vector{NBSElement}
    nlayer_total::Int
end

const _NBS_CAPTURE_EPS = 1e-12

# Number of appended NBS layer states for a rate-params object (0 when no NBS).
_nbs_state_count(p) = p.nbsplan === nothing ? 0 : p.nbsplan.nlayer_total

# Attenuate a signed correction `c` along `path_V` (oblivious runoff per cell) with the rule
# `max(V+c,0) - max(V,0)`: infiltration reduces the correction only where the baseline flow
# was there to reduce, so a negative correction is not eaten by empty-cell clamping and a
# positive one re-fills spare capacity first.  Empty `path_V` returns `c` unchanged.
function _propagate_correction(c::Float64, path_V::Vector{Float64})
    @inbounds for V in path_V
        c = max(V + c, 0.0) - max(V, 0.0)
    end
    return c
end

# Walk the flowgraph downstream from `start`, collecting the oblivious runoff `V` at each
# cell, until reaching a cell inside a network trap's footprint (the delivery target) or the
# domain edge (target 0, correction dropped).  Returns `(target_trap, path_V)`.
function _walk_to_trap(start::Int, tstruct, runoff, trap_of_cell::Dict{Int,Int})
    g     = tstruct.flowgraph
    pathV = Float64[]
    cur   = start
    while true
        haskey(trap_of_cell, cur) && return trap_of_cell[cur], pathV
        push!(pathV, Float64(runoff[cur]))
        ds = Graphs.outneighbors(g, cur)
        isempty(ds) && return 0, pathV
        cur = ds[1]
    end
end

# The footprint's terrain drainage endpoints, read from the oblivious runoff `runoff`.  A
# footprint cell whose flow leaves the footprint is a boundary exit; a footprint trap-bottom
# (flowgraph sink) is an internal depression.  Each yields `(Q_c, delivery)`: `Q_c` is the
# oblivious NBS drainage there, `delivery` is where `-X*Q_c` lands — the first downstream
# network trap (boundary, via `_walk_to_trap`) or the containing trap (internal, deposited
# directly).  Interior cells (draining to another footprint cell) are not endpoints.
function _nbs_endpoints(nb::DynNBSPlacement, tstruct, runoff, trap_of_cell::Dict{Int,Int})
    g       = tstruct.flowgraph
    footset = Set(nb.footprint)
    endpoints = Tuple{Float64,NBSDelivery}[]
    for f in nb.footprint
        Qf = max(Float64(runoff[f]), 0.0)
        Qf > 0.0 || continue
        ds = Graphs.outneighbors(g, f)
        if isempty(ds)                             # internal depression: ponds into its trap
            push!(endpoints, (Qf, NBSDelivery(get(trap_of_cell, f, 0), Float64[])))
        elseif ds[1] ∉ footset                     # boundary exit: flow crosses to ds[1]
            trap, pathV = _walk_to_trap(ds[1], tstruct, runoff, trap_of_cell)
            push!(endpoints, (Qf, NBSDelivery(trap, pathV)))
        end
    end
    return endpoints
end

# Piped outlets: each outflowing layer below the top `n_terrain` discharges at one explicit
# outlet cell.  Returns `(layer_index, delivery)` per outlet, the delivery walking from the
# outlet cell to its first downstream network trap.
function _nbs_outlets(nb::DynNBSPlacement, tstruct, runoff, trap_of_cell::Dict{Int,Int}, LI)
    outlets = Tuple{Int,NBSDelivery}[]
    piped = 0
    for (l, L) in enumerate(nb.system.layers)
        (l > nb.n_terrain && L.Kout > 0.0) || continue
        piped += 1
        trap, pathV = _walk_to_trap(LI[nb.outlets[piped]], tstruct, runoff, trap_of_cell)
        push!(outlets, (l, NBSDelivery(trap, pathV)))
    end
    return outlets
end

# Build the NBS plan for `net` (nothing when it has no NBS), reading the oblivious `runoff`
# grid to resolve each element's capture, drainage endpoints and outlet deliveries.  The
# layer models and their ODE-state layout come straight off `net.nbs`.
function _build_nbs_plan(net::DynNetwork, tstruct, runoff::AbstractMatrix{<:Real})
    isempty(net.nbs) && return nothing
    LI = LinearIndices(tstruct.topography)
    trap_of_cell = Dict{Int,Int}()             # linear cell -> local trap index
    for (ti, t) in enumerate(net.traps), c in tstruct.footprints[t.trap_ix]
        trap_of_cell[c] = ti
    end

    elems = NBSElement[]
    base  = 0
    for nb in net.nbs
        lyrs   = nb.system.layers
        A_foot = Float64(length(nb.footprint))   # @@@ 1 m^2/cell; use real cell area when available
        layers = NBSLayerParams[
            NBSLayerParams(float(L.Kout), float(L.nout), float(L.Smax),
                           float(L.Kinf), float(L.ninf), float(L.Smin), A_foot) for L in lyrs]
        endpoints = _nbs_endpoints(nb, tstruct, runoff, trap_of_cell)
        outlets   = _nbs_outlets(nb, tstruct, runoff, trap_of_cell, LI)
        capture   = sum(Float64[Qc for (Qc, _) in endpoints]; init = 0.0)
        push!(elems, NBSElement(layers, base, nb.n_terrain, capture, endpoints, outlets))
        base += length(lyrs)
    end
    return NBSPlan(elems, base)
end

# Add each NBS element's signed correction to the routed trap inflow.  The terrain re-emit
# reduction `-X*Q_c` (`X = 1 - O_terrain/I`) is propagated to each endpoint's target trap;
# the piped discharge is injected positively at each outlet's target.  `I = 0` (pure
# drainage) falls back to an equal split of `O_terrain` over the endpoints.
function _apply_nbs_corrections!(inflow, V, p, nt::Int)
    for el in p.nbsplan.elems
        base = el.state_base
        O_terrain = 0.0
        @inbounds for l in 1:el.n_terrain
            lp = el.layers[l]
            O_terrain += compute_outflow(lp.Kout, lp.nout, lp.Smax_mm,
                                         V[nt + base + l] * 1000.0 / lp.A) * 1e-3
        end
        I    = el.capture
        nend = length(el.endpoints)
        for (Qc, del) in el.endpoints
            del.trap == 0 && continue
            c = I > _NBS_CAPTURE_EPS ? (O_terrain / I - 1.0) * Qc : O_terrain / nend
            inflow[del.trap] += _propagate_correction(c, del.path_V)
        end
        for (l, del) in el.outlets
            del.trap == 0 && continue
            lp = el.layers[l]
            e  = compute_outflow(lp.Kout, lp.nout, lp.Smax_mm,
                                 V[nt + base + l] * 1000.0 / lp.A) * 1e-3
            inflow[del.trap] += _propagate_correction(e, del.path_V)
        end
    end
    return inflow
end

# ----------------------------------------------------------------------------
struct DynNetworkRateParams
    net::DynNetwork                       # The network being solved
    geom::Vector{TrapGeometry}            # struct holding info about trap's geometry
    external_inflow::Vector{Float64}      # constant inflow from terrain/weather per trap 
    external_path_inflow::Vector{Float64} # constant inflow (rain) on flow path cells
    footprint_infil::Vector{Float64}      # total footprint infiltration per trap (full trap loss)   
    path_infil_prefix::Vector{Vector{Float64}} # `cumsum` of infiltration along each flow path's cells
    order::Vector{Int}                    # topological sort of the flowgraph
    merge_target::Vector{Int}             # path that each tributary feeds into (0 if none)
    cvplan::Union{CulvertPlan,Nothing}    # key information per culvert for routing, or nothing if none
    path_events::Vector{Vector{Tuple{Int,Symbol,Int}}}  # in-order tributary/culvert stops per path
    nbsplan::Union{NBSPlan,Nothing}      # NBS signed-correction data, or nothing if none
end

# ----------------------------------------------------------------------------
"""
    _build_rate_params(tstruct, net, infiltration, external_inflow;
                       path_inflow=nothing, zvt=nothing) -> DynNetworkRateParams

Precompute the static [`DynNetworkRateParams`](@ref) for a solve.  `path_inflow` is constant
inflow per flow path (zeros if `nothing`); `zvt` a cached volume↔level table; `runoff` the
oblivious runoff grid, read only for the NBS plan (omit when the net has no NBS).
"""
function _build_rate_params(tstruct::TrapStructure,
                            net::DynNetwork,
                            infiltration::AbstractMatrix{<:Real},
                            external_inflow::AbstractVector{<:Real};
                            path_inflow = nothing,
                            runoff::Union{AbstractMatrix{<:Real},Nothing} = nothing,
                            zvt = nothing)
    nt = length(net.traps)
    np = length(net.flow_paths)
    @assert length(external_inflow) == nt
    path_inflow_vec = path_inflow === nothing ? zeros(Float64, np) : Float64.(path_inflow)
    @assert length(path_inflow_vec) == np
    order, _        = _network_order(net)
    cell_infil      = _path_cell_infiltration(net, infiltration)
    prefix          = [_infil_prefix(ci) for ci in cell_infil]
    cvplan  = isempty(net.culverts) ? nothing : _build_culvert_plan(net, tstruct)
    if isempty(net.nbs)
        nbsplan = nothing
    else
        runoff === nothing &&
            error("_build_rate_params: net has NBS but no runoff grid was supplied")
        nbsplan = _build_nbs_plan(net, tstruct, runoff)
    end
    # in-order tributary/culvert stops per path (tributaries only when there are no culverts)
    events = _path_event_templates(net)
    tgeom = _build_trap_geometry(tstruct, net, infiltration; zvt=zvt)
    footprint_infil = _footprint_infiltration(tgeom)
    return DynNetworkRateParams(net,
                                tgeom,
                                Float64.(external_inflow),
                                path_inflow_vec,
                                footprint_infil,
                                prefix,
                                order,
                                _merge_targets(net),
                                cvplan,
                                events,
                                nbsplan)
end

# ----------------------------------------------------------------------------
"""
    dynNetworkRateFunction!(dV, V, p::DynNetworkRateParams, t=0.0)

In-place ODE rate of the trap-volume state `V` (plus appended NBS layer states).

A *full* trap (`V >= capacity`) passes excess downstream at steady volume:
`dV = inflow - footprint_infil - max(inflow - footprint_infil, 0)` (0 while well fed, negative
once inflow drops below its losses).  An *accumulating* trap fills at its wetted-area rate:
`dV = inflow - wetted_infiltration(V)`.  `inflow` is the external inflow plus everything routed
in from upstream (see [`_route_flow`](@ref)).
"""
# Routed inflow into every trap at state `V`, plus which traps are full (spilling).  Shared by
# the rate function and the :unspill condition (which needs the raw net the rate fn clamps).
function _routed_inflow(V, p::DynNetworkRateParams)
    geom = p.geom
    nt   = length(geom)
    spilling = Bool[V[i] >= geom[i].capacity for i in 1:nt]
    # culverts need each trap's real water-surface elevation (not the Inf sentinel)
    trap_level = p.cvplan === nothing ? nothing :
                 Float64[_surface_level(geom[i], V[i]) for i in 1:nt]
    inflow = _route_flow(p.net, p.external_inflow, spilling,
                         p.footprint_infil, p.path_infil_prefix, p.external_path_inflow,
                         p.path_events, p.order, p.merge_target,
                         p.cvplan, trap_level)
    # NBS signed correction to the routed inflow (terrain re-emit reduction + piped discharge)
    p.nbsplan === nothing || _apply_nbs_corrections!(inflow, V, p, nt)
    return inflow, spilling
end

function dynNetworkRateFunction!(dV, V, p::DynNetworkRateParams, t = 0.0)
    geom = p.geom
    nt   = length(geom)

    inflow, spilling = _routed_inflow(V, p)

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

    # NBS layer-state derivatives (m^3/time).  Layer 1 (top) is fed by the oblivious
    # footprint capture `I`; each lower layer by the layer above's infiltration; the bottom
    # layer's infiltration leaves the system to ground.  Each layer loses its own overflow
    # and its infiltration.  Fluxes are power-law in the layer storage, in mm, converted to
    # m^3 (`S_mm = V*1000/A`, `Q_m3 = Q_mm*1e-3`).
    if p.nbsplan !== nothing
        @inbounds for el in p.nbsplan.elems
            base    = el.state_base
            prev_qi = el.capture                          # layer-1 inflow = capture I
            for (l, lp) in enumerate(el.layers)
                S_mm = V[nt + base + l] * 1000.0 / lp.A
                qo   = compute_outflow(lp.Kout, lp.nout, lp.Smax_mm, S_mm) * 1e-3
                qi   = compute_outflow(lp.Kinf, lp.ninf, lp.Smin_mm, S_mm) * 1e-3
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

# One monitored topology event condition (:fill/:empty/:unspill) on a network-local `trap`.
struct EventCondition
    kind::Symbol
    trap::Int
end

"""
    DynNetworkEvent

The event a solve terminated on: its `kind` (`:fill`, `:empty`, `:unspill`,
or `:none`) and the network-local `trap` index it concerns (0 if none).
The driver maps `trap` to a global trap index for its return value.
"""
mutable struct DynNetworkEvent
    kind::Symbol
    trap::Int
end
DynNetworkEvent() = DynNetworkEvent(:none, 0)

# ----------------------------------------------------------------------------
# The monitored conditions: :fill per evolving trap, :empty per evolving parent (a leaf at V=0
# is just its floor, no topology change), :unspill per full trap.  (LeftRootFind rationale: see
# the section header.)
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
# Steady-state detector (NBS-free): fires at the first accepted step where max|dV/dt| across
# the evolving traps drops below abstol.  Discrete because wetted_infiltration is a step
# function of V (interpolated mid-step states give spurious crossings).  Vetoes if any evolving
# trap is at V >= C (the :fill VCB terminates first).
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

# Steady-state detector with NBS: the `_routed_inflow` check can't express the layer dynamics,
# so evaluate the full rate function and settle when |dV/dt| < abstol (or a sign flip) across
# every `ss_indices` (evolving traps + NBS layer states).  A spilling trap still vetoes.
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

The `VectorContinuousCallback` (LeftRootFind) halting integration at the first
topology-changing event, and the [`DynNetworkEvent`](@ref) it records into.  `evolving` lists
the non-FULL traps; `nreg` restricts `:empty` to parent traps.
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
# Enforce the three-state caller contract (see [`solveDynNetwork!`](@ref)) before a solve;
# throws an informative error on the first violation.
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

# t=0 fast paths: an immediate event (no integration) when the entry state cannot be held,
# else `nothing`.  A FULL trap already draining (du0<0) unspills; a parent at its floor (V==0)
# with negative net inflow must expose its children now.  `du0`'s floor guard zeroes the rate
# at V==0, so the parent check recomputes the unclamped net rate.
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
                     tmax=Inf, path_inflow=nothing, runoff=nothing,
                     abstol=1e-6, reltol=1e-4, zvt=nothing) -> (; time, trap, kind)

Evolve the trap volumes of a [`DynNetwork`](@ref) forward until the first topology-changing
event, steady state, or `tmax`.  `state` is updated in place.

# Arguments
- `state`: stored volume (net of subtraps) per trap, `net.traps` order.  **Mutated in place**;
  must satisfy the three-state contract (below) on entry.
- `tstruct`: the [`TrapStructure`](@ref) supplying `net`'s geometry
- `net`: the network to solve
- `infiltration`: per-cell infiltration-rate grid (0 = impermeable)
- `inflow`: constant inflow per trap, `net.traps` order

# Keyword arguments
- `tmax`: integrate at most this far (elapsed).  With no event by then, stop at `tmax` and
  return `kind = :none` with `state` at `tmax` — lets a caller read the volumes at an
  externally-fixed time.  Default `Inf` (run to the first event or steady state).
- `path_inflow`: constant inflow onto each flow path (rain on path cells); zeros if `nothing`
- `runoff`: oblivious runoff grid, needed only when `net` has NBS
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

FULL traps that were neither the event trigger nor changed topology may have tiny
ODE-induced drift away from exact C.  The caller should clamp them to exactly C
(consistent with its authoritative `full_traps` list) before the next call.
"""
function solveDynNetwork!(state::AbstractVector{Float64},
                          tstruct::TrapStructure,
                          net::DynNetwork,
                          infiltration::AbstractMatrix{<:Real},
                          inflow::AbstractVector{<:Real};
                          tmax = Inf,
                          path_inflow = nothing,
                          runoff::Union{AbstractMatrix{<:Real},Nothing} = nothing,
                          # 1e-6 m^3 (~mL) / 1e-4: enough for physical accuracy, and ~halves the
                          # step count vs 1e-8 on the culvert worst case (see PARITY_TOL in tests)
                          abstol = 1e-6, reltol = 1e-4,
                          zvt = nothing)
    nt = length(net.traps)
    # Returns an (immutable) `DynNetworkRateParams` with all the static data the ODE needs
    p  = _build_rate_params(tstruct, net, infiltration, inflow;
                            path_inflow=path_inflow, runoff=runoff, zvt=zvt)
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
