import Interpolations
import Graphs
using DifferentialEquations   # brings DiscreteCallback, CallbackSet, SciMLBase into scope

export dynNetworkRateFunction!, solveDynNetwork!

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

        # Cells at/above the spillpoint never pond (water spills out instead), so they carry
        # no trap infiltration.  Excluding them keeps the wetted-area loss continuous at
        # V = capacity — otherwise a trap pinned at its spillpoint by a balancing through-flow
        # (e.g. a culvert outlet in drought) chatters full<->not-full.  Same rule as
        # `_ponding_infiltration` / `fill_trap_until`; test the actual terrain, not raised `bottom`.
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

A partially-filled trap infiltrates only through its wetted footprint, while a full
trap (`V >= capacity`, level `Inf`) infiltrates over its whole footprint.  Mirrors the
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
function _path_event_templates(net::DynNetwork, nbsplan = nothing)
    return [begin
                ev = Tuple{Int,Symbol,Int}[]
                for (m, j)   in fp.merges;          push!(ev, (j,   :trib,  m));  end
                for (ci, pos) in fp.culvert_inlets;  push!(ev, (pos, :cvin,  ci)); end
                for (ci, pos) in fp.culvert_outlets; push!(ev, (pos, :cvout, ci)); end
                if nbsplan !== nothing
                    for (pos, slot) in nbsplan.nbs_path_events[p]; push!(ev, (pos, :nbsout, slot)); end
                end
                sort!(ev; by = first)
                ev
            end
            for (p, fp) in enumerate(net.flow_paths)]
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

# Flow delivered at the end of a path: the head flow travels the path losing per-
# segment infiltration, tributaries add their delivered output at their junctions,
# and culvert inlets/outlets abstract/add at their cell positions.  `tribs` is the
# sorted `(junction, tributary)` list (culvert-free fast path); `events` is the
# precomputed `(pos, kind, idx)` stream (culvert-aware path) or `nothing`.  Mutates
# `culvert_actual` to record what each culvert inlet on this path drew.
function _path_delivered!(prefix, head_flow, tribs, events,
                          trib_output, culvert_actual, nbs_actual, cvplan, net, trap_level)
    current  = head_flow
    prev_pfx = 0.0
    if events === nothing
        for (junc, m) in tribs
            current  = max(current - (prefix[junc] - prev_pfx), 0.0)
            current += trib_output[m]
            prev_pfx = prefix[junc]
        end
    else
        for (pos, kind, idx) in events
            current = max(current - (prefix[pos] - prev_pfx), 0.0)
            if kind === :trib
                current += trib_output[idx]
            elseif kind === :cvout
                current += culvert_actual[idx]                # deliver
            elseif kind === :nbsout
                current += nbs_actual[idx]                    # NBS layer overflow at its exit/outlet
            else                                              # :cvin
                a = min(_culvert_flow(cvplan, net, idx, trap_level), current)
                culvert_actual[idx] = a                       # drawn == delivered
                current -= a
            end
            prev_pfx = prefix[pos]
        end
    end
    return max(current - (prefix[end] - prev_pfx), 0.0)
end

# Trap node: apply the culvert flows touching this trap (deliver what culverts
# routed in drew; drain what culverts draining it carry), then, if it is full,
# spill its surplus into its spill path.  Mutates `trap_inflow`, `path_flow`, and
# `culvert_actual`.
function _route_trap_node!(i, net, trap_inflow, path_flow, footprint_infil,
                           spilling, cvplan, trap_level, culvert_actual,
                           nbs_actual, nbs_trap_outlets)
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
    if nbs_trap_outlets !== nothing                   # deliver NBS overflow whose exit/outlet is this trap
        for slot in nbs_trap_outlets[i]
            trap_inflow[i] += nbs_actual[slot]
        end
    end
    if spilling[i]
        spill = max(trap_inflow[i] - footprint_infil[i], 0.0)
        sp = net.traps[i].spill_path
        sp > 0 && (path_flow[sp] += spill)            # sp == 0: spills out of domain
    end
    return nothing
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
                     path_events = nothing;
                     scratch = nothing,
                     nbs_actual = Float64[],
                     nbs_trap_outlets = nothing)

    np = length(net.flow_paths)
    nt = length(net.traps)
    @assert length(external_inflow) == length(spilling) == length(footprint_infil) == nt
    @assert length(path_infil_prefix) == length(merge_target) ==
            length(external_path_inflow) == length(sorted_trib_info) == np

    # Working accumulators.  With a caller-supplied `scratch` (the hot per-solve path)
    # they are reused in place; without one (the standalone/test path) they are freshly
    # allocated.  Either way: `trap_inflow` seeded with external trap inflow, `path_flow`
    # with per-path head inflow, `trib_output`/`culvert_actual` zeroed.
    if scratch === nothing
        trap_inflow = Float64.(external_inflow)
        path_flow   = Float64.(external_path_inflow)
        trib_output = zeros(Float64, np)
        culvert_actual = cvplan === nothing ? Float64[] : zeros(Float64, length(net.culverts))
    else
        trap_inflow    = copyto!(scratch.trap_inflow, external_inflow)
        path_flow      = copyto!(scratch.path_flow, external_path_inflow)
        trib_output    = fill!(scratch.trib_output, 0.0)
        culvert_actual = fill!(scratch.culvert_actual, 0.0)
    end

    for node in order
        if node <= np                               # a flow path
            p      = node
            events = path_events === nothing ? nothing : path_events[p]
            delivered = _path_delivered!(path_infil_prefix[p], path_flow[p],
                                         sorted_trib_info[p], events,
                                         trib_output, culvert_actual, nbs_actual,
                                         cvplan, net, trap_level)
            tt = net.flow_paths[p].target_trap
            if tt > 0
                trap_inflow[tt] += delivered        # into the downstream trap
            elseif merge_target[p] > 0
                trib_output[p] = delivered          # buffer for the main path's routing
            end                                      # else: exits the domain
        else                                        # a trap
            _route_trap_node!(node - np, net, trap_inflow, path_flow, footprint_infil,
                              spilling, cvplan, trap_level, culvert_actual,
                              nbs_actual, nbs_trap_outlets)
        end
    end
    return trap_inflow
end

# ----------------------------------------------------------------------------
"""
    _footprint_infiltration(geom) -> Vector{Float64}

Whole-footprint infiltration rate of each trap: the infiltration summed over every
cell of the trap's footprint (already net of the spillpoint-level cells excluded by
[`_build_trap_geometry`]).  This is the loss a full, submerged trap incurs, used by
[`_route_flow`](@ref).  Reads the per-cell `infil` from the precomputed geometry so it
stays consistent with `wetted_infiltration`.
"""
function _footprint_infiltration(geom::Vector{TrapGeometry})
    return [sum(g.infil) for g in geom]
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

# ----------------------------------------------------------------------------
# NBS terrain re-emit: the natural (no-NBS) exit split of a footprint.
#
# The top `n_terrain` layers of an NBS re-emit their overflow at the footprint's
# lower-edge exit boundary, spread the way the terrain would naturally spread it.
# This computes that static split: each footprint cell contributes a unit of
# drainage, followed downstream through
# the terrain flow tree until it first leaves the footprint; the first outside
# cell is that unit's exit.  A cell's weight is the fraction of footprint cells
# whose drainage exits there, so a funnel footprint (all flow leaving one cell)
# collapses to a single target and a fan spreads across several.
#
# `tstruct.flowgraph` is a flow tree — every cell has at most one downstream
# neighbour — so each unit's downstream walk is a simple chain with no branching
# and (the graph being a DAG) is guaranteed to terminate.  A footprint cell whose
# chain ends at a cell with no downstream neighbour is resolved the same way the
# `watercourses` sweep resolves such a cell: a trap bottom (`regions > 0`) is an
# internal pit — its water ponds inside the footprint and contributes no exit; a
# cell with `regions <= 0` drains off the domain edge, which is a valid terrain
# exit (the ordinary off-domain case) represented by the sentinel cell 0.
#
# Returns the exit cells (linear indices, all *outside* the footprint so the sink
# overlay in `watercourses` does not re-swallow the re-emit; the sentinel 0 marks
# the off-domain fraction) paired with weights summing to 1.  The off-domain
# weight stays in the sum so mass is conserved when it is later dropped at
# delivery.  Errors on a fully endorheic footprint (every cell ponds internally,
# no terrain or off-domain exit) — such a placement must route via piped outlets.
function _nbs_exit_weights(tstruct::TrapStructure, footprint::AbstractVector{<:Integer})
    g       = tstruct.flowgraph
    footset = Set{Int}(footprint)
    exits   = Dict{Int,Float64}()                # exit cell (0 = off-domain) -> unit count
    for c in footprint
        cur = Int(c)
        while cur in footset
            ds = Graphs.outneighbors(g, cur)
            isempty(ds) && break
            cur = ds[1]
        end
        if cur in footset                        # chain ended at a no-downstream footprint cell
            tstruct.regions[cur] > 0 && continue # trap bottom inside the footprint: ponds, no exit
            exits[0] = get(exits, 0, 0.0) + 1.0  # regions <= 0: drains off the domain edge
        else
            exits[cur] = get(exits, cur, 0.0) + 1.0   # first cell outside the footprint
        end
    end
    isempty(exits) &&
        error("_nbs_exit_weights: footprint has no terrain exit (fully endorheic); " *
              "its overflow cannot re-emit at terrain — use piped outlets instead")
    total = sum(values(exits))
    return sort!([(cell, w / total) for (cell, w) in exits]; by = first)
end

# ============================================================================
# Network rate function.
#
# `dynNetworkRateFunction!` is the in-place ODE right-hand side over the trap-
# volume state.  Its parameters (`DynNetworkRateParams`) bundle everything that is
# static for a solve
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
- `cvplan`, `path_events`: culvert routing data (or `nothing` when the net has no
  culverts)
- `scratch`: reusable per-solve working buffers (`RouteScratch`), so the hot rate
  function allocates nothing per call
"""
# Reusable per-solve scratch buffers for `_routed_inflow`/`_route_flow`.  The rate
# function is called very many times per solve (RK stages × many small steps), so
# allocating its working arrays afresh each call dominated GC; these are allocated once
# per solve (in `_build_rate_params`) and overwritten in place on every call.  Safe
# because every caller consumes the result immediately, the solver is single-threaded,
# and no call is re-entrant (see the callers of `_routed_inflow`).
struct RouteScratch
    spilling      ::Vector{Bool}
    trap_level    ::Vector{Float64}
    trap_inflow   ::Vector{Float64}
    path_flow     ::Vector{Float64}
    trib_output   ::Vector{Float64}
    culvert_actual::Vector{Float64}
    nbs_actual    ::Vector{Float64}    # per NBS overflow-delivery slot (like culvert_actual)
end
RouteScratch(nt::Int, np::Int, ncv::Int, nnbs::Int = 0) =
    RouteScratch(Vector{Bool}(undef, nt), Vector{Float64}(undef, nt),
                 Vector{Float64}(undef, nt), Vector{Float64}(undef, np),
                 Vector{Float64}(undef, np), Vector{Float64}(undef, ncv),
                 Vector{Float64}(undef, nnbs))

# ============================================================================
# NBS overlay routing plan.
#
# An NBS overlay element (a `DynNBSPlacement`) carries a vertical stack of storage layers
# whose state lives in the solver's ODE vector, appended after the `nt` trap states.
# Its topmost layer is fed by the static footprint capture (`nbs_inflow`, the sink-
# overlay tally from `watercourses`); each lower layer is fed by the infiltration of
# the layer above.  Each outflowing layer re-emits its overflow: the top `n_terrain`
# layers spread it over the footprint's natural terrain exit boundary (weights from
# `_nbs_exit_weights`), every outflowing layer below them through one explicit piped
# outlet.  Each re-emit landing point owned by an in-network path/trap gets a
# delivery slot in `nbs_actual` (filled from the layer storage each rate call, like
# `culvert_actual`); landings off the network / off the domain are dropped (that
# water leaves the domain).  Fluxes are power-law in the layer storage
# (`compute_outflow`), computed in the mm units the model is defined in with the
# mm<->m^3 conversion (`S_mm = V*1000/A`, `Q_m3 = Q_mm*1e-3`) done in the plan/rate
# function.
struct NBSLayerParams
    Kout::Float64; nout::Float64; Smax_mm::Float64
    Kinf::Float64; ninf::Float64; Smin_mm::Float64
    A::Float64                            # layer area (== footprint area, m^2)
end

struct NBSPlan
    placement_ix    ::Vector{Int}         # per NBS: its DynNBSPlacement `id` (== position in the
                                          # caller's vector); the key into `external_nbs_inflow`
    state_base      ::Vector{Int}         # per NBS: 0-based offset of its layer block, after the nt trap states
    layers          ::Vector{Vector{NBSLayerParams}}
    layer_slots     ::Vector{Vector{Vector{Tuple{Int,Float64}}}}  # per NBS, per layer: (slot, weight) deliveries
    nlayer_total    ::Int                 # total appended layer states
    n_slots         ::Int                 # number of on-network overflow-delivery slots
    nbs_path_events ::Vector{Vector{Tuple{Int,Int}}}  # per flow path: (cell position, slot)
    nbs_trap_outlets::Vector{Vector{Int}}             # per trap: slots delivered into it
    nbs_into        ::Vector{Vector{Int}}             # per NBS element: slots re-emitted onto its
                                                      # footprint by an upstream element (NBS→NBS
                                                      # capture, read as extra layer-1 inflow)
    # --- submergence ---
    n_terrain       ::Vector{Int}         # per NBS: size of the above-grade surface block (layers 1..n_terrain)
    containing_trap ::Vector{Int}         # per NBS: local trap covering the lowest footprint cell (0 = none)
    z_sub           ::Vector{Float64}     # per NBS: submergence threshold = lowest footprint elevation
    submerged       ::Vector{Bool}        # per NBS: current mode (surface block flooded / merged)
end

# Saturated-conduit intake of a (submerged) surface layer: once the surface is flooded,
# water is drawn down at the layer's infiltration
# *capacity*, evaluated at its own Smax (capacity-limited, constant, reuses the calibrated
# params in-range).  Returns m^3/time.  Isolated in one helper so switching to a
# flood-head-driven form later is a single-line change.
_nbs_saturated_draw(lp::NBSLayerParams) =
    compute_outflow(lp.Kinf, lp.ninf, lp.Smin_mm, lp.Smax_mm) * 1e-3

# Number of appended NBS layer states for a rate-params object (0 when no NBS).
_nbs_state_count(p) = p.nbsplan === nothing ? 0 : p.nbsplan.nlayer_total

# Build the NBS routing plan for `net` (nothing when it has no NBS elements).  Reads the
# `DynNBSPlacement`s straight off `net.nbs` — each supplies its own layer model (`nb.system`)
# and its stable inflow key (`nb.id`).  Each outflowing layer's overflow is delivered to
# resolved landing cells: the terrain exit boundary (top `n_terrain` layers, weighted by
# `_nbs_exit_weights`) or a piped outlet (lower layers, weight 1).  A landing owned by an
# in-network path gets a `(position, slot)` event; one owned by a trap gets a trap-outlet
# slot; a landing off the network or off the domain is dropped (its water exits).
function _build_nbs_plan(net::DynNetwork, tstruct,
                         submerged_of::Dict{Int,Bool} = Dict{Int,Bool}())
    isempty(net.nbs) && return nothing
    LI = LinearIndices(tstruct.topography)

    nbs_cell = Dict{Int,Int}()                 # linear footprint cell -> local NBS element index
    for (e, nb) in enumerate(net.nbs), k in nb.footprint
        nbs_cell[k] = e
    end
    trap_cell = Dict{Int,Int}()                # linear cell -> local trap index
    for (ti, t) in enumerate(net.traps), k in tstruct.footprints[t.trap_ix]
        trap_cell[k] = ti
    end
    path_cell = Dict{Int,Tuple{Int,Int}}()     # linear cell -> (local path index, position)
    for (pi, p) in enumerate(net.flow_paths), (pos, c) in enumerate(p.cells)
        path_cell[LI[c]] = (pi, pos)
    end

    np = length(net.flow_paths); ntr = length(net.traps); nnb = length(net.nbs)
    placement_ix = Int[]; state_base = Int[]
    layers       = Vector{NBSLayerParams}[]
    layer_slots  = Vector{Vector{Tuple{Int,Float64}}}[]
    nbs_path_events  = [Tuple{Int,Int}[] for _ in 1:np]
    nbs_trap_outlets = [Int[]            for _ in 1:ntr]
    nbs_into         = [Int[]            for _ in 1:nnb]
    n_terrain        = Int[]; containing_trap = Int[]; z_sub = Float64[]; submerged = Bool[]
    base = 0; slot = 0

    # Assign a delivery slot for fraction `w` of a layer's overflow landing at linear
    # cell `cell`; returns (slot, w) or nothing if the cell is off-network (discharge
    # then exits the domain and is simply not delivered).  A landing on another element's
    # footprint (NBS→NBS) is captured by that element: the slot is filled by this layer's
    # overflow and read as extra inflow by the receiver (see the rate function).  Footprint
    # capture wins over path/trap membership.
    function assign!(cell::Int, w::Float64)
        if haskey(nbs_cell, cell)
            slot += 1; push!(nbs_into[nbs_cell[cell]], slot); return (slot, w)
        elseif haskey(trap_cell, cell)
            slot += 1; push!(nbs_trap_outlets[trap_cell[cell]], slot); return (slot, w)
        elseif haskey(path_cell, cell)
            slot += 1; (pi, pos) = path_cell[cell]; push!(nbs_path_events[pi], (pos, slot)); return (slot, w)
        end
        return nothing
    end

    for nb in net.nbs
        lyrs   = nb.system.layers
        A_foot = Float64(length(nb.footprint))   # @@@ 1 m^2/cell; use real cell area when available
        push!(placement_ix, nb.id)               # nb.id is the key into external_nbs_inflow
        push!(state_base, base)
        push!(layers, NBSLayerParams[
            NBSLayerParams(float(L.Kout), float(L.nout), float(L.Smax),
                           float(L.Kinf), float(L.ninf), float(L.Smin), A_foot) for L in lyrs])

        exitw = _nbs_exit_weights(tstruct, nb.footprint)   # [(cell, weight)], cell 0 = off-domain
        per_layer = Vector{Tuple{Int,Float64}}[]
        piped = 0
        for (l, L) in enumerate(lyrs)
            s = Tuple{Int,Float64}[]
            if L.Kout > 0.0
                if l <= nb.n_terrain
                    for (cell, w) in exitw
                        cell == 0 && continue                 # off-domain fraction: exits
                        r = assign!(cell, w); r !== nothing && push!(s, r)
                    end
                else
                    piped += 1
                    r = assign!(LI[nb.outlets[piped]], 1.0); r !== nothing && push!(s, r)
                end
            end
            push!(per_layer, s)
        end
        push!(layer_slots, per_layer)

        # Submergence geometry: the surface floods when the containing trap's
        # water level reaches the lowest footprint cell.  The containing trap is the
        # network trap covering that lowest cell (0 if the footprint is not inside a
        # network trap — then submergence cannot arise).
        lowcell = nb.footprint[argmin(Float64[tstruct.topography[c] for c in nb.footprint])]
        push!(n_terrain,       nb.n_terrain)
        push!(containing_trap, get(trap_cell, lowcell, 0))
        push!(z_sub,           Float64(tstruct.topography[lowcell]))
        push!(submerged,       get(submerged_of, nb.id, false))

        base += length(lyrs)
    end
    return NBSPlan(placement_ix, state_base, layers, layer_slots, base, slot,
                   nbs_path_events, nbs_trap_outlets, nbs_into,
                   n_terrain, containing_trap, z_sub, submerged)
end

# Fill `nbs_actual[slot]` with each on-network delivery's share of its layer's overflow
# (m^3/time): `S_mm = V_layer*1000/A`, `Q_m3 = compute_outflow(...)*1e-3`, split across the
# layer's slots by their static weights.  Slots are dense (1..n_slots), so any unfilled
# (there are none: every slot belongs to exactly one (layer, landing)) would stay stale —
# the caller zeroes the buffer first for safety.
function _nbs_fill_actual!(nbs_actual, V, p, nt::Int)
    plan = p.nbsplan
    @inbounds for k in 1:length(plan.placement_ix)
        base = plan.state_base[k]
        for (l, lp) in enumerate(plan.layers[k])
            slots = plan.layer_slots[k][l]
            isempty(slots) && continue
            # A submerged surface layer (1..n_terrain) is underwater and does not re-emit
            # at terrain; its slots stay zero.  Subsurface (piped) layers keep discharging.
            (plan.submerged[k] && l <= plan.n_terrain[k]) && continue
            S_mm = V[nt + base + l] * 1000.0 / lp.A
            qo   = compute_outflow(lp.Kout, lp.nout, lp.Smax_mm, S_mm) * 1e-3
            for (slot, w) in slots
                nbs_actual[slot] = qo * w
            end
        end
    end
    return nothing
end

# ----------------------------------------------------------------------------
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
    nbsplan::Union{NBSPlan,Nothing}      # NBS overlay routing data, or nothing if none
    external_nbs_inflow::Vector{Float64} # per-NBS-placement static footprint capture (nbs_inflow)
    scratch::RouteScratch                # reusable per-solve working buffers
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
                            nbs_inflow::AbstractVector{<:Real} = Float64[],
                            nbs_submerged::Dict{Int,Bool} = Dict{Int,Bool}(),
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
    cvplan  = isempty(net.culverts) ? nothing : _build_culvert_plan(net, tstruct)
    nbsplan = isempty(net.nbs)      ? nothing : _build_nbs_plan(net, tstruct, nbs_submerged)
    # The event-driven router path is needed only for culverts (inlet draw / outlet
    # deliver at exact cell positions) or NBS overflow delivery at exit/outlet cells;
    # a net with neither needs no path events.
    events = (cvplan === nothing && nbsplan === nothing) ? nothing :
             _path_event_templates(net, nbsplan)
    tgeom = _build_trap_geometry(tstruct, net, infiltration; zvt=zvt)
    footprint_infil = _footprint_infiltration(tgeom)
    n_slots = nbsplan === nothing ? 0 : nbsplan.n_slots
    return DynNetworkRateParams(net,
                                tgeom,
                                Float64.(external_inflow),
                                path_inflow_vec,
                                footprint_infil,
                                prefix,
                                sorted_tribs,
                                order,
                                _merge_targets(net),
                                cvplan,
                                events,
                                nbsplan,
                                Float64.(nbs_inflow),
                                RouteScratch(nt, np, length(net.culverts), n_slots))
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
full traps (see [`_route_flow`](@ref)).
"""
# The routed inflow into every trap at state `V`, plus which traps are full
# (spilling).  Shared by the rate function and the :unspill event condition, which
# needs the raw net (inflow - footprint loss) that the rate function clamps away.
function _routed_inflow(V, p::DynNetworkRateParams)
    geom = p.geom
    nt   = length(geom)
    sc   = p.scratch
    spilling = sc.spilling
    @inbounds for i in 1:nt
        spilling[i] = V[i] >= geom[i].capacity
    end
    # culverts need each trap's real water-surface elevation (not the Inf sentinel)
    trap_level = nothing
    if p.cvplan !== nothing
        trap_level = sc.trap_level
        @inbounds for i in 1:nt
            trap_level[i] = _surface_level(geom[i], V[i])
        end
    end
    # NBS overflow delivered at exit/outlet cells: fill the slot buffer from the layer
    # storage (appended after the nt trap states), then route it in like culvert outlets.
    nbs_actual = sc.nbs_actual
    nbs_trap_outlets = nothing
    if p.nbsplan !== nothing
        fill!(nbs_actual, 0.0)
        _nbs_fill_actual!(nbs_actual, V, p, nt)
        nbs_trap_outlets = p.nbsplan.nbs_trap_outlets
    end
    inflow = _route_flow(p.net, p.external_inflow, spilling,
                         p.footprint_infil, p.path_infil_prefix, p.external_path_inflow,
                         p.sorted_trib_info, p.order, p.merge_target,
                         p.cvplan, trap_level, p.path_events; scratch = sc,
                         nbs_actual = nbs_actual, nbs_trap_outlets = nbs_trap_outlets)
    # Submerged NBS exchange with its containing trap: the captured footprint
    # runoff joins the flood and a saturated-conduit draw is taken back out.  Injected as
    # a trap-inflow adjustment (not a dV term) so the trap's own full/accumulating logic
    # handles it — a full flooding trap spills the surplus rather than over-filling.  The
    # subsurface then consumes `draw` in the rate function's layer loop.
    if p.nbsplan !== nothing
        plan = p.nbsplan
        @inbounds for k in 1:length(plan.placement_ix)
            (plan.submerged[k] && plan.containing_trap[k] > 0) || continue
            draw = plan.n_terrain[k] < length(plan.layers[k]) ?
                   _nbs_saturated_draw(plan.layers[k][plan.n_terrain[k]]) : 0.0
            inflow[plan.containing_trap[k]] += _nbs_captured(p, k) - draw
        end
    end
    return inflow, spilling
end

# Captured footprint inflow of NBS element `k`: the static footprint capture
# (`nbs_inflow`) plus any upstream element's overflow re-emitted onto this footprint
# (NBS→NBS, from the already-filled `nbs_actual` slots).
function _nbs_captured(p::DynNetworkRateParams, k::Int)
    plan = p.nbsplan
    c = p.external_nbs_inflow[plan.placement_ix[k]]
    @inbounds for s in plan.nbs_into[k]
        c += p.scratch.nbs_actual[s]
    end
    return c
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

    # NBS layer-state derivatives (m^3/time).  Layer 1 (top) is fed by the static
    # footprint capture (`nbs_inflow`); each lower layer by the layer above's
    # infiltration; the bottom layer's infiltration leaves the system to ground.  Each
    # layer loses its own overflow (delivered to its exit/outlet cells, already routed in
    # `_routed_inflow`) and its infiltration.  Fluxes are power-law in the layer storage,
    # in mm, converted to m^3 (`S_mm = V*1000/A`, `Q_m3 = Q_mm*1e-3`).
    if p.nbsplan !== nothing
        plan       = p.nbsplan
        nbs_actual = p.scratch.nbs_actual   # overflow slots, already filled in _routed_inflow
        @inbounds for k in 1:length(plan.placement_ix)
            base = plan.state_base[k]
            lyrs = plan.layers[k]
            nlyr = length(lyrs)

            if !plan.submerged[k]
                # --- dry: full cascade, top layer fed by the captured footprint inflow
                # (static capture + any upstream NBS→NBS re-emit, storage-driven). ---
                prev_qi = _nbs_captured(p, k)             # layer-1 inflow
                for (l, lp) in enumerate(lyrs)
                    S_mm = V[nt + base + l] * 1000.0 / lp.A
                    qo   = compute_outflow(lp.Kout, lp.nout, lp.Smax_mm, S_mm) * 1e-3
                    qi   = compute_outflow(lp.Kinf, lp.ninf, lp.Smin_mm, S_mm) * 1e-3
                    # @@@ evapotranspiration deferred: ET is an explicit 0.0 placeholder, so
                    #     wiring it in (from EVCoeff/EVS11) is a one-line change per layer, not
                    #     a restructuring.  Do not drop the term.
                    ET = 0.0
                    dV[nt + base + l] = prev_qi - qo - qi - ET
                    # Physical floor: a layer cannot drain below empty.
                    V[nt + base + l] <= 0.0 && dV[nt + base + l] < 0.0 && (dV[nt + base + l] = 0.0)
                    prev_qi = qi
                end
            else
                # --- submerged: the surface block (layers 1..n_terrain) is under
                # the flood and merged into the containing trap — frozen, no re-emit.  The
                # trap exchange (captured runoff in, saturated-conduit draw out) was applied
                # to the trap inflow in `_routed_inflow`; here the subsurface consumes `draw`
                # and cascades, its piped outlets still discharging. ---
                nte = plan.n_terrain[k]
                for l in 1:nte
                    dV[nt + base + l] = 0.0                # surface block frozen (merged)
                end
                if nte < nlyr
                    prev_qi = _nbs_saturated_draw(lyrs[nte])   # flood -> subsurface intake
                    for l in (nte + 1):nlyr
                        lp   = lyrs[l]
                        S_mm = V[nt + base + l] * 1000.0 / lp.A
                        qo   = compute_outflow(lp.Kout, lp.nout, lp.Smax_mm, S_mm) * 1e-3
                        qi   = compute_outflow(lp.Kinf, lp.ninf, lp.Smin_mm, S_mm) * 1e-3
                        ET   = 0.0
                        dV[nt + base + l] = prev_qi - qo - qi - ET
                        V[nt + base + l] <= 0.0 && dV[nt + base + l] < 0.0 && (dV[nt + base + l] = 0.0)
                        prev_qi = qi
                    end
                end
            end
        end
    end
    return nothing
end

# ============================================================================
# Event detection.
#
# Two callbacks stop the integration at the first event:
#
# 1. VectorContinuousCallback (`cb_topo`): per-trap topology-changing events.
#    Fires on downward zero-crossings (LeftRootFind) of:
#
#      :fill    capacity - V  -> 0   trap fills, starts spilling
#      :empty   V            -> 0   trap empties, exposes its subtraps
#      :unspill inflow - loss -> 0   a full trap's net inflow goes negative
#
#    Every *evolving* trap (V < capacity at solve start) gets :fill; parent
#    evolving traps also get :empty.  Every *full* trap gets :unspill.
#    LeftRootFind means conditions starting at exactly 0 (full traps with
#    zero initial net inflow) do not cause degenerate root-finding — they
#    wait silently for a genuine downward crossing (e.g. a culvert that later
#    makes the inflow positive and then negative mid-integration).
#
# 2. DiscreteCallback (`cb_ss`): global steady-state detector.
#    Fires at the first accepted ODE step where max(|dV/dt|) < abstol over
#    all evolving traps.  Using a DiscreteCallback (not ContinuousCallback)
#    is essential: wetted_infiltration is a step function of V, so evaluating
#    the condition at interpolated states mid-step would detect spurious
#    crossings.  Spilling-veto: if any evolving trap has V >= C at this step,
#    skip (the VCB :fill event fires first in the CallbackSet ordering).
# ============================================================================

# One monitored event condition.  For a topology event (:fill/:empty/:unspill) `trap`
# is the network-local trap and `nbs` is 0.  For an NBS submergence event (:submerge)
# `nbs` is the local NBS element index and `trap` is its containing trap; such an event
# does not terminate the solve — it flips the element's regime and merges the surface
# into the flood in place, then integration continues (see `_apply_submergence!`).
struct EventCondition
    kind::Symbol
    trap::Int
    nbs::Int
end
EventCondition(kind, trap) = EventCondition(kind, trap, 0)

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
# Build the list of monitored conditions.
#
# Each non-FULL (evolving) trap gets :fill.  Only evolving *parent* traps
# (trap_ix > nreg) get :empty — when a leaf trap reaches V=0 that is its physical
# floor with no topology change, so the event is not interesting.
# Every FULL (non-evolving) trap gets :unspill.  The VCB uses LeftRootFind so
# conditions starting at zero (full traps with zero initial net inflow) do not
# trigger degenerate root-finding; they wait silently for a genuine downward
# crossing (e.g. a culvert that later makes the inflow positive then negative).
# Steady-state termination is handled by a separate DiscreteCallback (`cb_ss`).
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
    # NBS submergence: one condition per element whose footprint sits in a network trap
    # and has an above-grade surface block.  It fires (non-terminating) each time the
    # containing trap's surface level crosses the footprint's flood threshold z_sub.
    if p.nbsplan !== nothing
        plan = p.nbsplan
        for k in 1:length(plan.placement_ix)
            (plan.containing_trap[k] > 0 && plan.n_terrain[k] >= 1) || continue
            push!(conds, EventCondition(:submerge, plan.containing_trap[k], k))
        end
    end
    return conds
end

# Set each NBS element's submerged flag from the current trap levels: submerged
# iff the containing trap's surface level has reached the footprint's flood threshold.
# Derived (not persisted) so it is always self-consistent with the committed state at the
# start of every solve; the in-solve :submerge event keeps it consistent across a window.
function _reconcile_submergence!(p::DynNetworkRateParams, V::AbstractVector{<:Real})
    p.nbsplan === nothing && return
    plan = p.nbsplan
    @inbounds for k in 1:length(plan.placement_ix)
        ct = plan.containing_trap[k]
        if ct > 0 && plan.n_terrain[k] >= 1
            plan.submerged[k] = _surface_level(p.geom[ct], V[ct]) >= plan.z_sub[k]
        else
            plan.submerged[k] = false
        end
    end
    return
end

# Apply an NBS submergence transition in place during the solve.  On submerge, the surface
# block merges into the containing trap — as much of its stored water as the trap's headroom
# allows moves into `u[ct]` (conserved; any residual stays frozen in the surface and is
# released on emergence), and the block is deactivated.  On emerge the surface simply
# re-activates from whatever storage it holds.  Non-terminating: integration continues.
function _apply_submergence!(u::AbstractVector, plan::NBSPlan, k::Int, nt::Int,
                             geom::Vector{TrapGeometry})
    if !plan.submerged[k]
        ct  = plan.containing_trap[k]
        base = plan.state_base[k]; nte = plan.n_terrain[k]
        St = 0.0
        @inbounds for l in 1:nte
            St += u[nt + base + l]
        end
        if St > 0.0
            moved = min(St, max(geom[ct].capacity - u[ct], 0.0))
            u[ct] += moved
            scale  = (St - moved) / St
            @inbounds for l in 1:nte
                u[nt + base + l] *= scale
            end
        end
        plan.submerged[k] = true
    else
        plan.submerged[k] = false
    end
    return
end

# ----------------------------------------------------------------------------
# Steady-state detector: DiscreteCallback that fires at the first accepted ODE
# step where the global maximum |dV/dt| across all evolving traps drops below
# abstol.
#
# A DiscreteCallback is the correct tool here because wetted_infiltration is a
# step function of V.  A ContinuousCallback would evaluate its condition at dense-
# output (polynomial) states that straddle cell-volume boundaries mid-step, and
# would detect spurious "crossings" of the discontinuous condition before the
# accepted step actually lands in the settled zone.  The DiscreteCallback only
# evaluates at accepted step endpoints, where the ODE state is genuine.
#
# Spilling-veto: if any evolving trap is at V >= capacity at this accepted step
# (the VCB :fill rootfinder is still narrowing down the crossing), fire only
# after the VCB has had a chance to terminate first.  In practice the VCB fires
# before this check, but the veto prevents a race at the exact crossing step.
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

# Steady-state callback for networks with NBS overlay elements.  The
# `_routed_inflow`-based check above cannot express the NBS layer dynamics, so this
# evaluates the full rate function and settles when |dV/dt| < abstol (or the rate has
# crossed its starting sign) across every index in `ss_indices` — the evolving traps
# plus all NBS layer states.  An evolving trap that is spilling still vetoes (its
# :fill wins).  `dbuf` is reused across calls.
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

Construct the `VectorContinuousCallback` (with `LeftRootFind`) that halts the
integration at the first topology-changing event and the
[`DynNetworkEvent`](@ref) it will record into.
`evolving` lists the network-local indices of non-FULL traps; `nreg` is
`numregions(tstruct)`, used to restrict `:empty` to parent traps.

Steady-state termination is handled by a separate `DiscreteCallback` built by
`_build_steadystate_callback`; the two are combined with `CallbackSet` in
`solveDynNetwork!`.
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
            elseif ec.kind == :submerge
                # crosses 0 downward at the transition, in whichever direction the current
                # regime makes it approach: dry -> submerge as the level rises past z_sub,
                # submerged -> emerge as it falls below.
                lvl = _surface_level(p.geom[ec.trap], V[ec.trap])
                out[k] = p.nbsplan.submerged[ec.nbs] ? (lvl - p.nbsplan.z_sub[ec.nbs]) :
                                                       (p.nbsplan.z_sub[ec.nbs] - lvl)
            else   # :unspill: fire when the full trap's net inflow drops below its losses
                out[k] = inflow[ec.trap] - p.footprint_infil[ec.trap]
            end
        end
    end

    function affect!(integrator, ix)
        k  = isa(ix, AbstractVector) ? findfirst(!iszero, ix) : ix
        ec = conds[k]
        if ec.kind == :submerge
            # regime switch: flip + merge surface into the flood in place, keep integrating.
            _apply_submergence!(integrator.u, p.nbsplan, ec.nbs, length(p.geom), p.geom)
            return
        end
        event.kind = ec.kind
        event.trap = ec.trap
        terminate!(integrator)
    end

    return VectorContinuousCallback(condition, affect!, length(conds);
                                    rootfind = SciMLBase.LeftRootFind), event
end

# ============================================================================
# Driver: solveDynNetwork.
#
# Wires the pieces together: build the static rate parameters, identify the
# evolving traps, integrate the trap-volume state forward as an `ODEProblem`, and
# stop at the first event.  The stop criterion is event-driven (a topology change)
# rather than a fixed time horizon, so no `tspan` is taken: integration runs to the
# first event, or to steady state (reported as an event time of `Inf`).
# ============================================================================

# ----------------------------------------------------------------------------
"""
    _validate_network(tstruct, net, state, geom)

Check the three-state caller contract before every [`solveDynNetwork!`](@ref) call.
Throws an informative error on the first violation found.

**Contract** (enforced here):
- FULL (V == C): `spill_path != 0` — either `> 0` (spills into that in-network
  path) or the sentinel `-1` (spills straight out of the domain).  If a parent
  trap, all its in-network children must also be FULL.
- TRANSITORY (0 < V < C): `spill_path == 0`.
- EMPTY (V == 0): `trap_ix <= numregions(tstruct)` (only leaf traps may be EMPTY).
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
            # This branch also covers a *parent* node sitting at its floor (V == 0):
            # a parent node subsumes its (full) children, so V == 0 is the parent's
            # own-volume floor with the children submerged — a valid TRANSITORY state,
            # not EMPTY.  It must, like any not-yet-full trap, have no downstream spill path.
            trap.spill_path == 0 ||
                error("solveDynNetwork!: TRANSITORY trap $(trap.trap_ix) has spill_path > 0. " *
                      "A trap that is not yet full should not have a downstream spill path. " *
                      "Rebuild the network without this trap in full_traps.")
        end
    end
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
                     tmax=Inf, path_inflow=nothing, abstol=1e-8, reltol=1e-8, zvt=nothing)
        -> (; time, trap, kind)

Evolve the water content of the lakes (traps) of a [`DynNetwork`](@ref) forward in
time until the first *event* that changes the network topology, or until the
elapsed time reaches `tmax`.  `state` is updated in-place with the trap volumes at
the event time (or at `tmax` if no event fires first).

# Arguments
- `state`: stored volume (net of subtraps) of each trap, indexed as `net.traps`.
  **Mutated in place.**  Must satisfy the three-state contract (see below) on entry.
- `tstruct`: the [`TrapStructure`](@ref) supplying the geometry `net` refers to
- `net`: the network to solve
- `infiltration`: grid (shape of `tstruct.topography`) of per-cell infiltration
  rates; 0 for impermeable / fully saturated cells
- `inflow`: constant inflow rate into each trap, indexed as `net.traps`

# Keyword arguments
- `tmax`: integrate at most this far (elapsed time from the solve start).  When no
  topology event fires before `tmax`, integration stops at `tmax` and the call
  returns `(time = tmax, trap = 0, kind = :none)` with `state` at `tmax`.  This
  lets a caller obtain the network volumes at an externally-determined time (e.g. a
  globally-committed event time, or a weather-period boundary) rather than only at
  the network's own next event.  Defaults to `Inf` (run to the first event or to
  steady state — the original behaviour).
- `path_inflow`: constant inflow rate directly onto each flow path (e.g. rainfall
  on path cells), indexed as `net.flow_paths`.  Defaults to zeros if `nothing`.
- `zvt`: pre-computed volume↔level tables from [`_compute_z_vol_tables`](@ref).
  Pass a cached value when solving many networks over the same [`TrapStructure`](@ref).
- `abstol`, `reltol`: ODE solver tolerances.  `abstol` also sets the |dV/dt|
  threshold of the steady-state cutoff.

# Returns
A named tuple (no `state` field — state is updated in place):
- `time`: time of the event; `Inf` if steady state was reached; `tmax` if the
  `tmax` window elapsed with no event.  `time` is elapsed from the solve start.
- `trap`: global trap index (into `tstruct`) where the event occurred, or `0`
- `kind`: `:fill`, `:empty`, `:unspill`, or `:none`

`kind == :none` covers both steady state (`time == Inf`) and a `tmax` cutoff
(`time == tmax`); they are distinguished by `time`.

# Three-state contract (validated at entry)

The caller must guarantee before every call:
- **FULL** (V == C exactly): `DynTrap.spill_path != 0` — `> 0` to spill into an
  in-network path, or the sentinel `-1` to spill out of the domain; if a parent
  trap, all its in-network children are also FULL.
- **TRANSITORY** (0 < V < C strictly): `DynTrap.spill_path == 0`.
- **EMPTY** (V == 0 exactly): the trap is a lowest-level (leaf) trap
  (`trap_ix <= numregions(tstruct)`).

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
                          nbs_inflow::AbstractVector{<:Real} = Float64[],
                          # Loosened from 1e-8: physical accuracy needs only ~mL (abstol, m^3)
                          # and ~ms.  ~Halves the ODE step count on the culvert worst case;
                          # fill-time drift vs the analytic path is ~3e-5 (tens of µs, far under
                          # 1 ms).  See PARITY_TOL in test/dynamics_test.jl.
                          abstol = 1e-6, reltol = 1e-4,
                          zvt = nothing)

    nt = length(net.traps)
    p  = _build_rate_params(tstruct, net, infiltration, inflow;
                            path_inflow=path_inflow, nbs_inflow=nbs_inflow, zvt=zvt)
    @assert length(state) == nt + _nbs_state_count(p) """
        state must have one entry per trap in net.traps plus one per NBS layer \
        ($(nt) traps + $(_nbs_state_count(p)) NBS layer states)"""
    V0 = copy(state)   # immutable snapshot; ODE evolves this, state is updated at the end

    _validate_network(tstruct, net, V0, p.geom)

    # Derive each NBS element's submerged regime from the current trap levels, so the
    # rate function starts consistent with the committed state regardless of the
    # dict passed in; the in-solve :submerge event keeps it consistent through the window.
    _reconcile_submergence!(p, V0)

    # Compute initial rates once: used for the t=0 fast-path checks.
    du0 = similar(V0, Float64)
    dynNetworkRateFunction!(du0, V0, p, 0.0)
    nreg = numregions(tstruct)

    # t=0 FULL→TRANSITORY fast path.  If any FULL trap already has du0 < 0 (net inflow
    # below its footprint losses), it begins draining right now.  Return :unspill
    # immediately without running the ODE.  `state` is left unchanged — the caller sets
    # state[trap] = prevfloat(C) before the next call to make it TRANSITORY.
    for i in 1:nt
        V0[i] == p.geom[i].capacity && du0[i] < 0.0 || continue
        return (time = 0.0, trap = net.traps[i].trap_ix, kind = :unspill)
    end

    # t=0 parent-EMPTY fast path.  A parent node sitting at its floor (V == 0, its
    # children submerged/full) whose net inflow is negative cannot keep that floor —
    # its children must be exposed now, so report :empty immediately.  The floor guard
    # zeroes `du0` at V == 0, so recompute the unclamped net rate here.
    inflow0, _ = _routed_inflow(V0, p)
    for i in 1:nt
        (V0[i] == 0.0 && net.traps[i].trap_ix > nreg) || continue
        inflow0[i] - wetted_infiltration(p.geom[i], 0.0) < 0.0 || continue
        return (time = 0.0, trap = net.traps[i].trap_ix, kind = :empty)
    end

    # Zero/negative window: nothing to advance.  `state` is left unchanged (it already
    # holds the volumes at the window's start) and no event is reported.
    tmax <= 0.0 && return (time = max(tmax, 0.0), trap = 0, kind = :none)

    # Traps that evolve: everything that is not FULL.  Empty traps with zero initial
    # inflow are included because inflow can start mid-integration (culverts, NBS
    # overflow, ...) without a separate event.  Excluding them would leave such traps
    # without a :fill callback and their volume could silently exceed capacity.
    # FULL traps (V == C) are handled only by the :unspill callback.
    evolving = [i for i in 1:nt if V0[i] != p.geom[i].capacity]

    # Steady-state test index set.  NBS layer states also evolve but never trigger a
    # topology event, so they enter the steady test only (appended after the nt traps).
    # Without NBS this is just `evolving`.
    ss_indices = p.nbsplan === nothing ? evolving :
        vcat(evolving, collect((nt + 1):(nt + p.nbsplan.nlayer_total)))

    # Nothing evolves (no non-full trap and no NBS layer state): steady state.
    isempty(ss_indices) && return (time = Inf, trap = 0, kind = :none)

    # All evolving states at or near zero rate: steady state already.
    # (abstol is a rate tolerance here, not a state classification guard.)
    all(abs(du0[i]) <= abstol for i in ss_indices) &&
        return (time = Inf, trap = 0, kind = :none)

    # Integrate to the first topology-changing event or steady state.
    # cb_topo (VectorContinuousCallback, LeftRootFind): topology events (:fill, :empty, :unspill).
    # cb_ss   (DiscreteCallback): fires at the first accepted step where max|dV/dt| < abstol
    #         across all evolving traps (and NBS layer states) — step-function infil means it
    #         must check at accepted steps, not interpolated points.
    # A trap-free (NBS-only) net has no topology events, so cb_topo is skipped to avoid a
    # zero-length VectorContinuousCallback.
    cb_ss = p.nbsplan === nothing ?
        _build_steadystate_callback(p, evolving, abstol, du0) :
        _build_steadystate_callback_nbs(p, ss_indices, abstol, du0)
    if nt == 0
        event    = DynNetworkEvent()
        callback = cb_ss
    else
        cb_topo, event = _build_event_callback(p, evolving, V0, nreg)
        callback = CallbackSet(cb_topo, cb_ss)
    end
    # Cap tmax at a large but finite value so DiffEq never receives Inf as a tspan
    # endpoint (which triggers an internal range(t, Inf, n) that Julia rejects).  1e12
    # is many orders of magnitude beyond any physical simulation horizon.
    tmax_ode = min(tmax, 1e12)
    # Explicit `Tsit5`: the rate function is non-smooth (culvert regime switches + routing
    # `min` are C0 kinks) but NOT stiff.  DiffEq's default auto-switcher misreads the kinks
    # as stiffness, goes implicit, and pays for dense finite-diff Jacobians on the whole
    # coupled system — ~4x slower on full-coverage+culvert, no accuracy gain.  Explicit just
    # steps across the kinks and wins.
    # @@@ If a genuinely stiff config appears (e.g. huge culvert draining a tiny trap),
    #     revisit with an auto-switching solver + a sparse/colored Jacobian, not dense FD.
    sol = solve(ODEProblem(dynNetworkRateFunction!, V0, (0.0, tmax_ode), p), Tsit5();
                callback = callback,
                abstol = abstol, reltol = reltol)

    # Write ODE result back into state in place (saves one nt-length allocation per call).
    state .= sol.u[end]

    # No topology event fired.  Two sub-cases, distinguished by the stop time:
    #   - reached the `tmax_ode` cutoff (sol.t[end] ≈ tmax_ode): the window elapsed with
    #     the network still evolving — report the caller's `tmax`, with `state` at tmax_ode.
    #   - cb_ss terminated earlier (sol.t[end] < tmax_ode): genuine steady state — `Inf`.
    if event.kind == :none
        return sol.t[end] >= tmax_ode ? (time = tmax, trap = 0, kind = :none) :
                                        (time = Inf,  trap = 0, kind = :none)
    end

    # Clamp the triggering trap to its exact threshold so the validator passes on the
    # next call without requiring the caller to handle floating-point residual at the
    # event crossing.
    if event.kind == :fill
        state[event.trap] = p.geom[event.trap].capacity
    elseif event.kind == :empty
        state[event.trap] = 0.0
    end

    return (time = sol.t[end],
            trap = net.traps[event.trap].trap_ix,
            kind = event.kind)
end
