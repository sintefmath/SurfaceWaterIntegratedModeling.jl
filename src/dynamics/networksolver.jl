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
        for m in path.merges
            Graphs.add_edge!(g, m, p)            # tributary m feeds path p
        end
    end
    for (t, trap) in enumerate(net.traps)
        trap.spill_path > 0 && Graphs.add_edge!(g, np + t, trap.spill_path)
    end

    return Graphs.topological_sort_by_dfs(g), np
end

# ----------------------------------------------------------------------------
# Reverse of the `merges` lists: `merge_target[p]` is the path that tributary
# path `p` feeds into, or 0 if `p` is not a tributary.  Each path is a tributary
# of at most one other (paths truncated by `_resolve_cell_overlaps!`), so this is
# well defined.
function _merge_targets(net::DynNetwork)
    merge_target = zeros(Int, length(net.flow_paths))
    for (a, path) in enumerate(net.flow_paths), m in path.merges
        merge_target[m] = a
    end
    return merge_target
end

# ----------------------------------------------------------------------------
"""
    _route_flow(net, external_inflow, spilling, footprint_infil, path_infil;
                path_inflow=zeros(np)) -> Vector{Float64}

Total inflow arriving at each trap of `net` (in `net.traps` order): the caller-
supplied `external_inflow` (per trap) and `path_inflow` (per path) plus everything
routed in from upstream.

Each *spilling* trap (`spilling[i] == true`) emits `max(inflow - footprint_infil[i],
0)` into its spill path; non-spilling traps pass nothing on.  Flow travelling down
a path is seeded with the path's `path_inflow`, loses `path_infil[p]` to
infiltration (floored at zero), then is delivered to the path's `target_trap`, or
into the path it merges into, or out of the domain.

`footprint_infil[i]` is trap `i`'s whole-footprint infiltration (a full trap is
submerged over its entire footprint) and `path_infil[p]` is path `p`'s total
infiltration capacity; both are passed in so this routine stays pure topology and
is testable on a hand-built network.  See [`_footprint_infiltration`](@ref) and
[`_path_infiltration`](@ref) for computing them from a `TrapStructure`.

!!! note "Merge approximation"
    A tributary's delivered flow is added at the head of the path it merges into,
    so it is charged that path's full infiltration capacity rather than only the
    portion below the junction (the junction cell is not retained when paths are
    truncated).  Exact without merges; a slight over-infiltration with them.
"""
function _route_flow(net::DynNetwork,
                     external_inflow::AbstractVector{<:Real},
                     spilling::AbstractVector{Bool},
                     footprint_infil::AbstractVector{<:Real},
                     path_infil::AbstractVector{<:Real};
                     path_inflow = zeros(length(net.flow_paths)))
    order, _ = _network_order(net)
    return _route_flow(net, external_inflow, spilling, footprint_infil, path_infil,
                       path_inflow, order, _merge_targets(net))
end

# Core router with the routing plan (`order` from [`_network_order`](@ref) and
# `merge_target` from [`_merge_targets`](@ref)) supplied.  Both are static for a
# given network, so the rate function precomputes them once rather than rebuilding
# the graph and topological sort on every ODE step.
function _route_flow(net::DynNetwork,
                     external_inflow::AbstractVector{<:Real},
                     spilling::AbstractVector{Bool},
                     footprint_infil::AbstractVector{<:Real},
                     path_infil::AbstractVector{<:Real},
                     external_path_inflow::AbstractVector{<:Real},
                     order::AbstractVector{<:Integer},
                     merge_target::AbstractVector{<:Integer})

    np = length(net.flow_paths)
    nt = length(net.traps)
    @assert length(external_inflow) == length(spilling) == length(footprint_infil) == nt
    @assert length(path_infil) == length(merge_target) == length(external_path_inflow) == np

    trap_inflow = Float64.(external_inflow)        # accumulator, seeded with external trap inflow
    path_flow   = Float64.(external_path_inflow)   # accumulator, seeded with external path inflow

    for node in order
        if node <= np                                   # a flow path
            p = node
            delivered = max(path_flow[p] - path_infil[p], 0.0)
            tt = net.flow_paths[p].target_trap
            if tt > 0
                trap_inflow[tt] += delivered            # into the downstream trap
            elseif merge_target[p] > 0
                path_flow[merge_target[p]] += delivered # into the path it merges into
            end                                          # else: leaves the domain
        else                                            # a trap
            i = node - np
            if spilling[i]
                spill = max(trap_inflow[i] - footprint_infil[i], 0.0)
                sp = net.traps[i].spill_path
                sp > 0 && (path_flow[sp] += spill)      # sp == 0: spills out of domain
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
    _path_infiltration(net, infiltration) -> Vector{Float64}

Total infiltration capacity of each flow path in `net` (in `net.flow_paths`
order): the infiltration grid summed over the path's cells.  Used by
[`_route_flow`](@ref) as the in-transit loss along the path.
"""
function _path_infiltration(net::DynNetwork, infiltration::AbstractMatrix{<:Real})
    return [isempty(p.cells) ? 0.0 : sum(infiltration[c] for c in p.cells)
            for p in net.flow_paths]
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
- `path_infil`: total infiltration capacity per flow path
- `order`, `merge_target`: the static routing plan (see [`_network_order`](@ref),
  [`_merge_targets`](@ref))
"""
struct DynNetworkRateParams
    net::DynNetwork
    geom::Vector{TrapGeometry}
    external_inflow::Vector{Float64}
    external_path_inflow::Vector{Float64}
    footprint_infil::Vector{Float64}
    path_infil::Vector{Float64}
    order::Vector{Int}
    merge_target::Vector{Int}
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
    order, _ = _network_order(net)
    return DynNetworkRateParams(net,
                                _build_trap_geometry(tstruct, net, infiltration; zvt=zvt),
                                Float64.(external_inflow),
                                path_inflow_vec,
                                _footprint_infiltration(tstruct, net, infiltration),
                                _path_infiltration(net, infiltration),
                                order,
                                _merge_targets(net))
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
function dynNetworkRateFunction!(dV, V, p::DynNetworkRateParams, t = 0.0)
    geom = p.geom
    nt   = length(geom)
    @assert length(dV) == length(V) == nt

    # which traps are currently full (and therefore spilling / passing water on)
    spilling = Vector{Bool}(undef, nt)
    @inbounds for i in 1:nt
        spilling[i] = V[i] >= geom[i].capacity
    end

    inflow = _route_flow(p.net, p.external_inflow, spilling,
                         p.footprint_infil, p.path_infil, p.external_path_inflow,
                         p.order, p.merge_target)

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
# Only *evolving* traps (V < capacity at solve start) are monitored: a full trap
# already sits at `capacity - V == 0`, which would trigger :fill spuriously.
#
# Dormant scaffolds (:unspill when a spilling trap's net crosses zero, and
# culvert-outlet threshold crossings) plug into the same vector but produce no
# entries until a state-dependent drain (culvert / NBS) exists.
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

    # --- dormant scaffolding ---------------------------------------------------
    # :unspill -- a full trap stops spilling when its net (inflow - losses) crosses
    # zero.  With constant inflow and infiltration this net is constant, so it can
    # only happen once a trap has a state-dependent drain (a culvert outlet).
    for (i, trap) in enumerate(p.net.traps)
        isempty(trap.culvert_outlets) || push!(conds, EventCondition(:unspill, i))
    end
    # Culvert-outlet threshold crossings would be appended here per culvert; there
    # are none in the current scope (lakes and rivers only), so nothing is added.

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
        dynNetworkRateFunction!(dubuf, V, p, t)
        @inbounds for (k, ec) in enumerate(conds)
            if ec.kind == :fill
                out[k] = p.geom[ec.trap].capacity - V[ec.trap]
            elseif ec.kind == :empty
                out[k] = V[ec.trap]
            elseif ec.kind == :stagnation
                out[k] = dv0[ec.trap] * dubuf[ec.trap]
            else
                out[k] = 1.0   # dormant scaffolds: never cross in the current scope
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
