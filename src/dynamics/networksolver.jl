import Interpolations
import Graphs

export dynNetworkRateFunction!, solveDynNetwork

# ============================================================================
# Geometry helpers for the dynamic network solver.
#
# A `DynNetwork` carries only topology; the geometry it refers to lives in the
# `TrapStructure`.  These helpers precompute, once per solve, the per-trap
# geometry the rate function needs: storage capacity, the footprint cells, their
# terrain bottoms, the per-cell infiltration rate, and the stored-volume -> water-
# level relationship.  From those we can answer the two questions the rate
# function repeatedly asks of each trap:
#
#   * water_level(g, V)          -- the water surface elevation at stored volume V
#   * wetted_infiltration(g, V)  -- the infiltration loss at stored volume V, i.e.
#                                   the infiltration rate summed over only the
#                                   currently-submerged footprint cells.
#
# This reuses the existing machinery: `_compute_z_vol_tables` for the volume<->
# level tables and the `infilfun(trap_bottom .<= z)` wetted-area pattern from
# `_setup_dvdt` / `fill_trap_until` in fill_sequence.jl.
# ============================================================================

"""
    TrapGeometry

Precomputed geometry of a single trap participating in a [`DynNetwork`](@ref),
used by the dynamic network rate function.

Volumes are *net of subtraps* (the convention used throughout fill_sequence and
by [`waterheight`](@ref)): `capacity` is the trap's own storage above its
subtraps, and the stored-volume state `V` ranges over `0 .. capacity`.

# Fields
- `trap_ix`: index of the trap in the [`TrapStructure`](@ref)
- `capacity`: own storage volume, net of subtraps
- `footprint`: linear indices of the trap's footprint cells
- `bottom`: terrain bottom at each footprint cell.  For a parent trap the cells
  lying over a subtrap are raised to the subtrap's spillpoint elevation, so that
  a submerged subtrap reads as wetted (matching `_compute_z_vol_tables`).
- `infil`: infiltration rate at each footprint cell (0 where impermeable)
- `zmin`: lowest bottom elevation (water level of an empty trap)
- `v2z`: stored-volume -> water-level interpolation (see [`water_level`](@ref))
"""
struct TrapGeometry
    trap_ix::Int
    capacity::Float64
    footprint::Vector{Int}
    bottom::Vector{Float64}
    infil::Vector{Float64}
    zmin::Float64
    v2z   # volume -> level map; untyped (perf is a later concern, correctness first)
end

# ----------------------------------------------------------------------------
"""
    _build_trap_geometry(tstruct, net, infiltration; zvt=nothing) -> Vector{TrapGeometry}

Precompute a [`TrapGeometry`](@ref) for every trap in the [`DynNetwork`](@ref)
`net`, in the same order as `net.traps` (i.e. indexed by the network-local trap
index, not the global trap index).

`infiltration` is a grid the same shape as `tstruct.topography`, giving a
per-cell infiltration rate (0 for impermeable / fully saturated cells).

`zvt` is the pre-computed volume↔level table from [`_compute_z_vol_tables`](@ref).
If `nothing` (default) it is computed here; pass a cached value when calling
repeatedly for the same [`TrapStructure`](@ref) to avoid redundant work.
"""
function _build_trap_geometry(tstruct::TrapStructure,
                              net::DynNetwork,
                              infiltration::AbstractMatrix{<:Real};
                              zvt = nothing)

    # Reuse the existing per-trap volume<->level tables.  These are computed for
    # every trap in the structure; we only index the ones present in `net`.
    zvt === nothing && (zvt = _compute_z_vol_tables(tstruct))
    nreg = numregions(tstruct)

    geom = Vector{TrapGeometry}(undef, length(net.traps))
    for (i, dtrap) in enumerate(net.traps)
        tix       = dtrap.trap_ix
        footprint = tstruct.footprints[tix]
        bottom    = Float64.(tstruct.topography[footprint])

        # Mirror fill_trap_until / _compute_z_vol_tables: for a parent trap the
        # bottom is raised to the upper subtrap's spillpoint, and the capacity is
        # the trap volume net of its subtraps.  A lowest-level trap (region) has
        # no subtraps, so its capacity is just its trap volume.
        if tix > nreg
            children = subtrapsof(tstruct, tix)
            bottom   = max.(bottom, tstruct.spillpoints[children[1]].elevation)
            capacity = Float64(tstruct.trapvolumes[tix] - tstruct.subvolumes[tix])
        else
            capacity = Float64(tstruct.trapvolumes[tix])
        end

        zmin    = minimum(bottom)
        zvtable = zvt[tix]
        # Volume -> level map (mirrors _setup_dvdt).  A single-row table means a
        # degenerate flat trap with no storage; level is then always the bottom.
        v2z = length(zvtable[2]) == 1 ?
            (_ -> zmin) :
            Interpolations.linear_interpolation(zvtable[2], zvtable[1],
                                                extrapolation_bc = Interpolations.Line())

        infil = Float64.(infiltration[footprint])

        geom[i] = TrapGeometry(tix, capacity, footprint, bottom, infil, zmin, v2z)
    end
    return geom
end

# ----------------------------------------------------------------------------
"""
    water_level(g::TrapGeometry, V) -> Float64

Water surface elevation of trap `g` holding stored volume `V` (net of subtraps).

An empty trap (`V <= 0`) sits at its bottom `zmin`; a full trap (`V >= capacity`)
is reported as `Inf` so that its entire footprint reads as submerged.  In between,
the level is interpolated from the volume<->level table.  Mirrors the `z`
computation in `_setup_dvdt`.
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

Infiltration loss rate of trap `g` holding stored volume `V`: the per-cell
infiltration rate summed over only the currently-submerged footprint cells
(those whose terrain bottom lies at or below the current water level).

This realises the wetted-area dependence the prompt calls for: a partially-filled
trap infiltrates only through its wetted footprint, while a full trap (`V >=
capacity`, level `Inf`) infiltrates over its whole footprint.  Mirrors the
`infilfun(trap_bottom .<= z)` term in `_setup_dvdt`.
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
# Flow paths are not stateful: they route water downstream instantaneously, but
# they are lossy -- a path cell with an infiltration value removes water in
# transit (capped by the flow actually passing through).  This mirrors the
# throughput-capped logic of `_update_runoff!`, lifted from the grid to the
# `DynNetwork` abstraction.
#
# A full (spilling) trap emits an output = inflow - footprint infiltration, which
# is routed down its spill path; a not-yet-full (leaf / accumulating) trap keeps
# its inflow and passes nothing on.  Tributaries (`DynFlowPath.merges`) deliver
# their flow into the path they merge into.  The single quantity the rate
# function needs out of all this is the total inflow arriving at each trap.
# ============================================================================

# Combined topological order over the network's path and trap nodes.  Nodes
# 1:np are flow paths, np+1:np+nt are traps.  Edges follow the flow: a trap to
# its spill path, a path to its target trap, and a tributary path to the path it
# merges into.  Returns (order, np); `node <= np` is path `node`, otherwise it is
# trap `node - np`.
#
# `net.traps` alone is ordered upstream->downstream, but that does not fix the
# trap/path interleaving once tributaries are present (two traps feeding the same
# path are unordered relative to each other), so we sort the combined graph.
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
# Reverse of the `merges` lists: `merge_target[p]` is the path that tributary
# path `p` feeds into, or 0 if `p` is not a tributary.  Each path is a tributary
# of at most one other (paths truncated by `_resolve_cell_overlaps!`), so this is
# well defined.
function _merge_targets(net::DynNetwork)
    merge_target = zeros(Int, length(net.flow_paths))
    for (a, path) in enumerate(net.flow_paths), (m, _) in path.merges
        merge_target[m] = a
    end
    return merge_target
end

# ----------------------------------------------------------------------------
# Culvert routing data and helpers.
#
# A culvert connects an inlet endpoint to an outlet endpoint, each owned by a flow
# path or a trap (resolved at construction; see the DynFlowPath/DynTrap culvert
# lists).  `CulvertPlan` precomputes, per culvert, which kind of element owns each
# endpoint and its network-local index, plus the barrel diameter and a handle on the
# topography for the invert-elevation lookups used by `culvert_rate`.

struct CulvertPlan
    topo::AbstractMatrix          # tstruct.topography (invert elevations)
    diam::Vector{Float64}         # barrel diameter (2r) per culvert
    inlet_is_path::BitVector
    inlet_owner::Vector{Int}      # local path or trap index owning the inlet
    outlet_is_path::BitVector
    outlet_owner::Vector{Int}     # local path or trap index owning the outlet
end

function _build_culvert_plan(net::DynNetwork, tstruct)
    nc   = length(net.culverts)
    plan = CulvertPlan(tstruct.topography, [2 * cv.r for cv in net.culverts],
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

# Per-path in-order event template for the culvert-aware router: tributary junctions
# and culvert inlet/outlet positions merged and sorted by cell position, as
# `(position, :trib|:cvin|:cvout, tributary-or-culvert index)`.  The list is static
# for a solve (only the flow values it drives are dynamic), so it is built once and
# reused every rate-function call instead of being rebuilt and re-sorted each time.
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
        ih   = max(trap_level[plan.inlet_owner[ci]] - plan.topo[cv.inlet], 0.0)
        isub = ih >= D
    end
    if plan.outlet_is_path[ci]
        oh, osub = D, false
    else
        oh   = max(trap_level[plan.outlet_owner[ci]] - plan.topo[cv.outlet], 0.0)
        osub = oh >= D
    end
    return culvert_rate(cv, (; topography = plan.topo);
                        inlet_submerged = isub,  inlet_head  = ih,
                        outlet_submerged = osub, outlet_head = oh,
                        allow_reverse = false)
end

# ----------------------------------------------------------------------------
"""
    _route_flow(net, external_inflow, spilling, footprint_infil, path_cell_infil;
                path_inflow=zeros(np)) -> Vector{Float64}

Total inflow arriving at each trap of `net` (in `net.traps` order): the caller-
supplied `external_inflow` (per trap) and `path_inflow` (per path, entering at
each path's inlet) plus everything routed in from upstream.

Each *spilling* trap emits `max(inflow - footprint_infil[i], 0)` into its spill
path.  Flow along a path is routed in segments between tributary junctions: the
head flow (from trap spills and `path_inflow`) travels through the pre-junction
cells losing infiltration, then each tributary's delivered flow is added at its
exact junction position, and the combined flow continues through the remaining
cells.  This gives each tributary the correct post-junction infiltration charge.

`path_cell_infil[p]` is a vector of per-cell infiltration rates for path `p`
(in cell order); both are passed in so this routine stays pure topology and
is testable on a hand-built network.  See [`_path_cell_infiltration`](@ref)
and [`_footprint_infiltration`](@ref) for computing them from a `TrapStructure`.
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
    sorted_tribs  = [sort([(j, m) for (m, j) in fp.merges]) for fp in net.flow_paths]
    path_events   = cvplan === nothing ? nothing : _path_event_templates(net)
    return _route_flow(net, external_inflow, spilling, footprint_infil,
                       prefix, path_inflow, sorted_tribs, order, _merge_targets(net),
                       cvplan, trap_level, path_events)
end

# Core router with all static routing data pre-supplied.  Tributaries are handled
# via segmented routing: flow is charged infiltration only over the path cells it
# actually travels through, so a tributary joining at junction j is not charged
# the pre-junction infiltration of the main path.  With no tributaries the loop
# over `sorted_tribs[p]` is empty and the result is identical to a simple lump
# deduction — one code path for both cases.
function _route_flow(net::DynNetwork,
                     external_inflow::AbstractVector{<:Real},
                     spilling::AbstractVector{Bool},
                     footprint_infil::AbstractVector{<:Real},
                     path_infil_prefix::AbstractVector{<:AbstractVector{<:Real}},
                     external_path_inflow::AbstractVector{<:Real},
                     sorted_trib_info::AbstractVector{<:AbstractVector},
                     order::AbstractVector{<:Integer},
                     merge_target::AbstractVector{<:Integer},
                     cvplan = nothing,
                     trap_level = nothing,
                     path_events = nothing)

    np = length(net.flow_paths)
    nt = length(net.traps)
    @assert length(external_inflow) == length(spilling) == length(footprint_infil) == nt
    @assert length(path_infil_prefix) == length(merge_target) ==
            length(external_path_inflow) == length(sorted_trib_info) == np

    trap_inflow = Float64.(external_inflow)        # accumulator, seeded with external trap inflow
    path_flow   = Float64.(external_path_inflow)   # head inflow per path (trap spills + path_inflow)
    trib_output = zeros(Float64, np)               # delivered output of each tributary path
    # actual flow each culvert carries (drawn at inlet == delivered at outlet)
    culvert_actual = cvplan === nothing ? Float64[] : zeros(Float64, length(net.culverts))

    for node in order
        if node <= np                               # a flow path
            p        = node
            prefix   = path_infil_prefix[p]
            fp       = net.flow_paths[p]
            current  = path_flow[p]
            prev_pfx = 0.0

            if cvplan === nothing
                # Segmented routing over tributary junctions only: each segment loses
                # its infiltration, then the tributary's delivered flow is added.
                for (junc, m) in sorted_trib_info[p]
                    current  = max(current - (prefix[junc] - prev_pfx), 0.0)
                    current += trib_output[m]
                    prev_pfx = prefix[junc]
                end
            else
                # Walk the precomputed in-order stream of tributary junctions and
                # culvert inlet/outlet positions along the path.  A culvert inlet
                # abstracts up to the flow passing its cell (mass cap); a culvert
                # outlet adds the amount its (earlier-routed) source actually drew.
                for (pos, kind, idx) in path_events[p]
                    current = max(current - (prefix[pos] - prev_pfx), 0.0)
                    if kind === :trib
                        current += trib_output[idx]
                    elseif kind === :cvout
                        current += culvert_actual[idx]            # deliver
                    else                                          # :cvin
                        a = min(_culvert_flow(cvplan, net, idx, trap_level), current)
                        culvert_actual[idx] = a                   # drawn == delivered
                        current -= a
                    end
                    prev_pfx = prefix[pos]
                end
            end

            delivered = max(current - (prefix[end] - prev_pfx), 0.0)
            tt = fp.target_trap
            if tt > 0
                trap_inflow[tt] += delivered        # into the downstream trap
            elseif merge_target[p] > 0
                trib_output[p] = delivered          # buffer for main path's segmented routing
            end                                      # else: exits the domain
        else                                        # a trap
            i = node - np
            if cvplan !== nothing
                trap = net.traps[i]
                # culvert outlets ending in this trap: deliver what their (already-
                # routed, earlier in topological order) source drew.
                for ci in trap.culvert_outlets
                    trap_inflow[i] += culvert_actual[ci]
                end
                # culvert inlets draining this trap (the matching delivery happens at
                # the outlet's owner -- trap or flow path -- later in topological order)
                for ci in trap.culvert_inlets
                    q = _culvert_flow(cvplan, net, ci, trap_level)
                    culvert_actual[ci] = q
                    trap_inflow[i]    -= q              # drawn out of this trap
                end
            end
            if spilling[i]
                spill = max(trap_inflow[i] - footprint_infil[i], 0.0)
                sp = net.traps[i].spill_path
                sp > 0 && (path_flow[sp] += spill)  # sp == 0: spills out of domain
            end
        end
    end
    return trap_inflow
end

# ----------------------------------------------------------------------------
"""
    _footprint_infiltration(tstruct, net, infiltration) -> Vector{Float64}

Whole-footprint infiltration rate of each trap in `net` (in `net.traps` order):
the infiltration grid summed over every cell of the trap's footprint.  This is the
loss a full, submerged trap incurs, used by [`_route_flow`](@ref).
"""
function _footprint_infiltration(tstruct::TrapStructure,
                                 net::DynNetwork,
                                 infiltration::AbstractMatrix{<:Real})
    return [sum(@view infiltration[tstruct.footprints[t.trap_ix]]) for t in net.traps]
end

# ----------------------------------------------------------------------------
"""
    _path_cell_infiltration(net, infiltration) -> Vector{Vector{Float64}}

Per-cell infiltration rate of each flow path in `net` (in `net.flow_paths` order):
the infiltration grid evaluated at each cell along the path, preserving cell order.
Used by [`_build_rate_params`](@ref) to compute prefix sums for [`_route_flow`](@ref).
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
# `dynNetworkRateFunction!` is the in-place ODE right-hand side over the trap-
# volume state, analogous to `NBSNetworkRateFunction!` in SUrbArea's NBS.jl.  Its
# parameters (`DynNetworkRateParams`) bundle everything that is static for a solve
# -- geometry, the constant external inflow, the infiltration sums, and the
# routing plan -- so the per-step work is just routing + a loop over traps.
# ============================================================================

"""
    DynNetworkRateParams

Static, precomputed inputs to [`dynNetworkRateFunction!`](@ref) for one solve.
Built once by [`_build_rate_params`](@ref) and passed as the `ODEProblem`
parameter.  Indexed network-locally (same order as `net.traps` / `net.flow_paths`).

# Fields
- `net`: the [`DynNetwork`](@ref) being solved
- `geom`: per-trap [`TrapGeometry`](@ref) (capacity + wetted-area infiltration)
- `external_inflow`: caller-supplied constant inflow per trap
- `external_path_inflow`: caller-supplied constant inflow directly onto each flow
  path (e.g. rainfall on path cells); added at the path head before infiltration
- `footprint_infil`: whole-footprint infiltration per trap
- `path_infil_prefix`: per-path prefix sums of per-cell infiltration
  (`prefix[k]` = sum of cell infiltrations 1..k-1, so `prefix[1]=0` and
  `prefix[end]` = total path infiltration)
- `sorted_trib_info`: per-path list of `(junction_cell_index, tributary_path_index)`
  sorted by junction position; empty for paths with no tributaries
- `order`, `merge_target`: the static routing plan (see [`_network_order`](@ref),
  [`_merge_targets`](@ref))
"""
struct DynNetworkRateParams
    net::DynNetwork
    geom::Vector{TrapGeometry}
    external_inflow::Vector{Float64}
    external_path_inflow::Vector{Float64}
    footprint_infil::Vector{Float64}
    path_infil_prefix::Vector{Vector{Float64}}
    sorted_trib_info::Vector{Vector{Tuple{Int,Int}}}
    order::Vector{Int}
    merge_target::Vector{Int}
    cvplan::Union{CulvertPlan,Nothing}   # culvert routing data, or nothing if none
    path_events::Union{Vector{Vector{Tuple{Int,Symbol,Int}}},Nothing}  # per-path event templates
end

# ----------------------------------------------------------------------------
"""
    _build_rate_params(tstruct, net, infiltration, external_inflow;
                       path_inflow=nothing, zvt=nothing) -> DynNetworkRateParams

Precompute the static [`DynNetworkRateParams`](@ref) for a solve: the per-trap
geometry, the footprint/path infiltration sums, and the routing plan.

`path_inflow` is the constant external inflow per flow path (e.g. rainfall on path
cells).  Defaults to zeros if `nothing`.

`zvt` is a pre-computed volume↔level table (see [`_build_trap_geometry`](@ref)).
Pass a cached value when solving many networks over the same [`TrapStructure`](@ref).
"""
function _build_rate_params(tstruct::TrapStructure,
                            net::DynNetwork,
                            infiltration::AbstractMatrix{<:Real},
                            external_inflow::AbstractVector{<:Real};
                            path_inflow = nothing,
                            zvt = nothing)
    nt = length(net.traps)
    np = length(net.flow_paths)
    @assert length(external_inflow) == nt
    path_inflow_vec = path_inflow === nothing ? zeros(Float64, np) : Float64.(path_inflow)
    @assert length(path_inflow_vec) == np
    order, _        = _network_order(net)
    cell_infil      = _path_cell_infiltration(net, infiltration)
    prefix          = [_infil_prefix(ci) for ci in cell_infil]
    sorted_tribs    = [sort([(j, m) for (m, j) in fp.merges]) for fp in net.flow_paths]
    cvplan = isempty(net.culverts) ? nothing : _build_culvert_plan(net, tstruct)
    events = cvplan === nothing ? nothing : _path_event_templates(net)
    return DynNetworkRateParams(net,
                                _build_trap_geometry(tstruct, net, infiltration; zvt=zvt),
                                Float64.(external_inflow),
                                path_inflow_vec,
                                _footprint_infiltration(tstruct, net, infiltration),
                                prefix,
                                sorted_tribs,
                                order,
                                _merge_targets(net),
                                cvplan,
                                events)
end

# ----------------------------------------------------------------------------
"""
    dynNetworkRateFunction!(dV, V, p::DynNetworkRateParams, t=0.0)

In-place ODE rate of the trap-volume state `V`, writing dV/dt into `dV`.

Each trap is either *full* (spilling) or *accumulating*.  A trap is taken to be
full when `V[i] >= capacity`; it then passes its excess inflow downstream and its
own volume is steady, so

    dV = inflow - footprint_infiltration - spill,  spill = max(inflow - footprint_infiltration, 0)

(which is 0 while the trap is adequately fed, and negative once its inflow drops
below its losses -- the onset of draining).  An accumulating trap spills nothing
and fills at its wetted-area rate

    dV = inflow - wetted_infiltration(V).

`inflow` is the caller's constant inflow plus everything routed in from upstream
full traps (see [`_route_flow`](@ref)).  Mirrors `NBSNetworkRateFunction!`.
"""
# The routed inflow into every trap at state `V`, plus which traps are full
# (spilling).  Shared by the rate function and the :unspill event condition, which
# needs the raw net (inflow - footprint loss) that the rate function clamps away.
function _routed_inflow(V, p::DynNetworkRateParams)
    geom = p.geom
    nt   = length(geom)
    spilling = Vector{Bool}(undef, nt)
    @inbounds for i in 1:nt
        spilling[i] = V[i] >= geom[i].capacity
    end
    # culverts need each trap's real water-surface elevation (not the Inf sentinel)
    trap_level = p.cvplan === nothing ? nothing :
                 [_surface_level(geom[i], V[i]) for i in 1:nt]
    inflow = _route_flow(p.net, p.external_inflow, spilling,
                         p.footprint_infil, p.path_infil_prefix, p.external_path_inflow,
                         p.sorted_trib_info, p.order, p.merge_target,
                         p.cvplan, trap_level, p.path_events)
    return inflow, spilling
end

function dynNetworkRateFunction!(dV, V, p::DynNetworkRateParams, t = 0.0)
    geom = p.geom
    nt   = length(geom)
    @assert length(dV) == length(V) == nt

    inflow, spilling = _routed_inflow(V, p)

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
    end
    return nothing
end

# ============================================================================
# Event detection.
#
# The solve integrates forward and stops at the first event.  Events are detected
# with a `VectorContinuousCallback` whose condition vector generalises the
# `[full, empty, stagnation]` triple of `fill_trap_until` to every evolving trap:
#
#   :fill        capacity - V        -> 0   trap fills, starts spilling (type-1)
#   :empty       V                   -> 0   trap empties, exposes its subtraps
#   :stagnation  dv0 * dV            -> 0   net rate changed sign: steady state.
#                                           Reachable in finite time because the
#                                           wetted-area infiltration is a step
#                                           function of V, so dV is piecewise
#                                           constant and crosses zero at a cell
#                                           elevation rather than asymptotically.
#
# Only the *evolving* traps (V < capacity at solve start) get the fill/empty/
# stagnation triple: a full trap already sits at `capacity - V == 0`, which would
# trigger :fill spuriously.  Each *full* trap instead gets an :unspill condition --
#
#   :unspill     inflow - loss       -> 0   a full trap's net inflow crosses zero
#                                           (downward) and it begins draining; e.g.
#                                           a feeding/draining culvert or (later)
#                                           an NBS inflow changing with the state.
#
# It never crosses while the net is constant, so it is inert for plain constant-
# inflow networks.
# ============================================================================

# One monitored event condition, tied to a network-local trap index.
struct EventCondition
    kind::Symbol
    trap::Int
end

"""
    DynNetworkEvent

The event a solve terminated on: its `kind` (`:fill`, `:empty`, `:stagnation`,
`:unspill`, or `:none`) and the network-local `trap` index it concerns (0 if
none).  The driver maps `trap` to a global trap index for its return value.
"""
mutable struct DynNetworkEvent
    kind::Symbol
    trap::Int
end
DynNetworkEvent() = DynNetworkEvent(:none, 0)

# ----------------------------------------------------------------------------
# Build the list of monitored conditions: the live triple for each evolving trap,
# plus dormant scaffolds that stay empty until drainage / culverts exist.
function _event_conditions(p::DynNetworkRateParams, evolving::AbstractVector{<:Integer})
    conds = EventCondition[]
    for e in evolving
        push!(conds, EventCondition(:fill,       e))
        push!(conds, EventCondition(:empty,      e))
        push!(conds, EventCondition(:stagnation, e))
    end

    # :unspill -- a full (non-evolving) trap stops spilling and begins draining when
    # its net inflow (inflow - footprint loss) crosses zero.  Its net can vary because
    # a feeding/draining culvert depends on the upstream state, and upstream spills
    # (and, later, NBS inflows) can drop too.  Harmless when the net is constant: the
    # condition simply never crosses.
    evolving_set = Set(evolving)
    for i in eachindex(p.net.traps)
        i in evolving_set || push!(conds, EventCondition(:unspill, i))
    end

    return conds
end

# ----------------------------------------------------------------------------
"""
    _build_event_callback(p, evolving, V0) -> (callback, event)

Construct the `VectorContinuousCallback` that halts the integration at the first
event and the [`DynNetworkEvent`](@ref) it will record into.  `evolving` lists the
network-local indices of the traps that evolve (those with `V < capacity` at the
start); `V0` is the initial state, used to fix the sign of each evolving trap's
rate for the stagnation test.
"""
function _build_event_callback(p::DynNetworkRateParams,
                               evolving::AbstractVector{<:Integer},
                               V0::AbstractVector{<:Real})

    conds = _event_conditions(p, evolving)
    event = DynNetworkEvent()

    # initial rate, to detect sign changes (stagnation) later
    dv0 = similar(V0, Float64)
    dynNetworkRateFunction!(dv0, V0, p, 0.0)
    dubuf = similar(V0, Float64)   # reused derivative scratch

    function condition(out, V, t, integrator)
        # one routing pass, reused for the stagnation rate (dubuf) and the :unspill
        # net (inflow - footprint loss), which the clamped rate would otherwise hide.
        inflow, spilling = _routed_inflow(V, p)
        @inbounds for i in eachindex(V)
            dubuf[i] = spilling[i] ?
                inflow[i] - p.footprint_infil[i] - max(inflow[i] - p.footprint_infil[i], 0.0) :
                inflow[i] - wetted_infiltration(p.geom[i], V[i])
        end
        @inbounds for (k, ec) in enumerate(conds)
            if ec.kind == :fill
                out[k] = p.geom[ec.trap].capacity - V[ec.trap]
            elseif ec.kind == :empty
                out[k] = V[ec.trap]
            elseif ec.kind == :stagnation
                out[k] = dv0[ec.trap] * dubuf[ec.trap]
            elseif ec.kind == :unspill
                # net inflow of a full trap; crosses zero (downward) as it starts draining
                out[k] = inflow[ec.trap] - p.footprint_infil[ec.trap]
            else
                out[k] = 1.0
            end
        end
    end

    function affect!(integrator, ix)
        # DiffEq v8+ may pass ix as a vector of flags; older versions a scalar index
        k = isa(ix, AbstractVector) ? findfirst(!iszero, ix) : ix
        event.kind = conds[k].kind
        event.trap = conds[k].trap
        terminate!(integrator)
    end

    return VectorContinuousCallback(condition, affect!, length(conds)), event
end

# ============================================================================
# Driver: solveDynNetwork.
#
# Wires the pieces together: build the static rate parameters, identify the
# evolving traps, integrate the trap-volume state forward as an `ODEProblem`, and
# stop at the first event.  Mirrors `solveNBSNetwork` in SUrbArea's NBS.jl, but
# the stop criterion is event-driven (a topology change) rather than a fixed time
# horizon, so no `tspan` is taken: integration runs to the first event, or to
# steady state (reported as an event time of `Inf`).
# ============================================================================

"""
    solveDynNetwork(tstruct, net, infiltration, inflow, state;
                    path_inflow=nothing, abstol=1e-8, reltol=1e-8, zvt=nothing)
        -> (; time, trap, kind, state)

Evolve the water content of the lakes (traps) of a [`DynNetwork`](@ref) forward in
time until the first *event* that changes the network topology, and return it.

# Arguments
- `tstruct`: the [`TrapStructure`](@ref) supplying the geometry `net` refers to
- `net`: the network to solve
- `infiltration`: grid (shape of `tstruct.topography`) of per-cell infiltration
  rates; 0 for impermeable / fully saturated cells
- `inflow`: constant inflow rate into each trap, indexed as `net.traps`
- `state`: current stored volume (net of subtraps) of each trap, indexed as
  `net.traps`.  Non-leaf traps must be full (at capacity); leaf traps may be
  partially filled.

# Keyword arguments
- `path_inflow`: constant inflow rate directly onto each flow path (e.g. rainfall
  on path cells), indexed as `net.flow_paths`.  Defaults to zeros if `nothing`.
  Added at the path head before infiltration losses are applied.
- `zvt`: pre-computed volume↔level tables from [`_compute_z_vol_tables`](@ref).
  Pass a cached value when solving many networks over the same [`TrapStructure`](@ref)
  to avoid recomputing these static tables on every call.
- `abstol`, `reltol`: ODE solver tolerances.

# Returns
A named tuple:
- `time`: time of the event, or `Inf` if the network reaches steady state with no
  further event
- `trap`: global trap index (into `tstruct`) where the event occurred, or `0`
  for the steady-state / no-event case
- `kind`: `:fill` (a trap filled and started spilling), `:empty` (a trap emptied),
  `:unspill` (a full trap's losses exceed its inflow and it immediately begins
  draining), or `:none` (steady state)
- `state`: the trap-volume state at `time`

The problem is integrated as an `ODEProblem` with [`dynNetworkRateFunction!`](@ref)
as the right-hand side and the first-event callback from
[`_build_event_callback`](@ref).
"""
function solveDynNetwork(tstruct::TrapStructure,
                         net::DynNetwork,
                         infiltration::AbstractMatrix{<:Real},
                         inflow::AbstractVector{<:Real},
                         state::AbstractVector{<:Real};
                         path_inflow = nothing,
                         abstol = 1e-8, reltol = 1e-8,
                         zvt = nothing)

    nt = length(net.traps)
    @assert length(state) == nt "state must have one entry per trap in net.traps"

    p  = _build_rate_params(tstruct, net, infiltration, inflow;
                            path_inflow=path_inflow, zvt=zvt)
    V0 = Float64.(state)

    # Compute initial rates once; used for both the :unspill check and the stagnation check.
    du0 = similar(V0, Float64)
    dynNetworkRateFunction!(du0, V0, p, 0.0)

    # A full trap whose losses exceed its inflow (du0 < 0) begins draining immediately.
    # Return an :unspill event at t=0 and clamp all such traps to just below capacity
    # so the next call treats them as evolving rather than re-triggering the same event.
    unspilling = [i for i in 1:nt if V0[i] >= p.geom[i].capacity && du0[i] < -abstol]
    if !isempty(unspilling)
        worst = unspilling[argmin(du0[unspilling])]
        final = copy(V0)
        for i in unspilling
            final[i] = prevfloat(p.geom[i].capacity)
        end
        return (time = 0.0, trap = net.traps[worst].trap_ix, kind = :unspill, state = final)
    end

    # Traps that evolve: those below their capacity.  Full traps with non-negative net
    # (handled above) are steady (dV ~ 0) and would trigger :fill spuriously.
    evolving = [i for i in 1:nt if V0[i] < p.geom[i].capacity]

    # Nothing evolves, or everything that could is already at rest: steady state.
    isempty(evolving) && return (time = Inf, trap = 0, kind = :none, state = V0)
    all(abs(du0[i]) <= abstol for i in evolving) &&
        return (time = Inf, trap = 0, kind = :none, state = V0)

    # Integrate to the first event.  Events occur in finite time (the wetted-area
    # infiltration is piecewise constant, so a leaf either fills or hits a steady
    # level), so an unbounded horizon with the terminating callback is safe.
    cb, event = _build_event_callback(p, evolving, V0)
    sol = solve(ODEProblem(dynNetworkRateFunction!, V0, (0.0, Inf), p);
                callback = cb, abstol = abstol, reltol = reltol)
    final = sol.u[end]

    # Stagnation (or no trigger) means steady state -> no further event.
    if event.kind == :stagnation || event.kind == :none
        return (time = Inf, trap = 0, kind = :none, state = final)
    end

    # Clamp the triggering trap to its exact threshold so the next call does not
    # re-trigger the same event due to floating-point residual at the crossing.
    final = copy(final)
    if event.kind == :fill
        final[event.trap] = p.geom[event.trap].capacity
    elseif event.kind == :empty
        final[event.trap] = 0.0
    end

    return (time  = sol.t[end],
            trap  = net.traps[event.trap].trap_ix,
            kind  = event.kind,
            state = final)
end
