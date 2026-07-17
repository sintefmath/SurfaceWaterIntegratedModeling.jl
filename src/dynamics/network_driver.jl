# ----------------------------------------------------------------------------
# The incremental network driver: builds the components once per weather period, then mutates
# the live `Vector{DynNetwork}` in place per event (`apply_fill!`/`apply_unfill!`/`apply_empty!`)
# instead of retracing.  `vol_by_trapix` (committed node volumes) and `nbs_state` (NBS layers)
# are the source of truth surviving structural change; contexts are rebuilt from them.
# ----------------------------------------------------------------------------

"""
    NetworkDriver

The live dynamic-network state for one weather period of a `fill_sequence` run: the component
set, one [`DynNetworkContext`](@ref) per component, and the committed state that survives a
structural change.  Built by [`build_network_driver`](@ref).

Fields:
- `comps`: live components, mutated in place by `apply_fill!` / `apply_unfill!` /
  `apply_empty!` — the structural source of truth (culverts and NBS travel on each component,
  not on the driver).
- `contexts`: one context per entry of `comps`, same order; rebuilt after every event.
- `vol_by_trapix`: node trap index -> committed volume at `last_time`.
- `nbs_state`: NBS placement id -> layer states, persistent across events and periods.
- `last_time`: absolute time the committed volumes hold at.
- `precip`: the weather period's rain rate, handed to each context for the NBS plan.  Constant
  for the driver's lifetime — a driver is built per weather period.
"""
mutable struct NetworkDriver
    comps          ::Vector{DynNetwork}          # live components, mutated by apply_* (the
                                                 # structural source of truth — culverts/NBS
                                                 # travel on each component, not the driver)
    contexts       ::Vector{DynNetworkContext}   # one per comp, same order (rebuilt each event)
    vol_by_trapix  ::Dict{Int,Float64}           # node trap_ix -> committed volume at last_time
    nbs_state      ::Dict{Int,Vector{Float64}}   # placement id -> layer states (persistent)
    last_time      ::Float64                      # absolute time the committed volumes hold at
    precip         ::Union{Matrix{Float64},Float64}  # the period's rain rate, for the NBS plan
end

# ----------------------------------------------------------------------------

"""
    build_network_driver(tstruct, dyn_coords, culverts, full_traps, cur_amounts, rateinfo,
                         infiltration, z_vol_tables, cur_time, endtime;
                         nbs_placements = DynNBSPlacement[], precipitation = 0.0,
                         nbs_state = Dict{Int,Vector{Float64}}()) -> NetworkDriver

Build the [`NetworkDriver`](@ref) for one weather period: trace the components from the seeds,
seed the committed volumes from `cur_amounts`, and build and predict one context per component
over `[cur_time, endtime]`.

# Arguments
- `tstruct`: the [`TrapStructure`](@ref) supplying terrain, traps and spillpoints.
- `dyn_coords`, `culverts`, `nbs_placements`: the dynamic seeds to trace from.
- `full_traps`: indices of traps currently at capacity.
- `cur_amounts`, `rateinfo`, `infiltration`, `z_vol_tables`: the fill-sequence state supplying
  each node's initial volume and inflow.
- `cur_time`, `endtime`: the period the driver runs over.
- `precipitation`: the period's rain rate; required (with `rateinfo.runoff`) when there are NBS
  placements, since each one's throughput is the rain on its footprint plus what flows in.
- `nbs_state`: existing NBS layer states, carried in from an earlier period.

# Returns
The driver, with its contexts built and predicted.

`cur_amounts` and `full_traps` are read only, but `nbs_placements` is mutated by
[`setup_network`](@ref), and the caller's `nbs_state` is retained by reference — the driver
writes each period's layer storage back through it.
"""
function build_network_driver(tstruct, dyn_coords, culverts, full_traps, cur_amounts,
                              rateinfo, infiltration, z_vol_tables, cur_time, endtime;
                              nbs_placements = DynNBSPlacement[],
                              precipitation  = 0.0,
                              nbs_state      = Dict{Int,Vector{Float64}}())
    comps = setup_network(tstruct, full_traps;
                          dyn_coords = dyn_coords, culverts = culverts, nbs = nbs_placements)

    # committed volume of every node from the caller's amounts (all current at build, no projection)
    vol = Dict{Int,Float64}()
    for net in comps, t in net.traps
        vol[t.trap_ix] = cur_amounts[t.trap_ix].amount
    end

    precip = precipitation isa Real ? Float64(precipitation) : Matrix{Float64}(precipitation)
    driver = NetworkDriver(comps, DynNetworkContext[], vol, nbs_state, cur_time, precip)
    _rebuild_contexts!(driver, tstruct, rateinfo, infiltration, z_vol_tables,
                       cur_amounts, Set(full_traps), cur_time, endtime)
    return driver
end

"""
    _driver_state0(driver, tstruct, rateinfo, cur_amounts, full_set, cur_time, z_vol_tables)
        -> Function

Build the `state0(g)` closure that seeds a context with trap `g`'s committed volume at
`cur_time`, from three sources in order: a full trap sits at its exact own capacity `C`; a node
already in `driver.vol_by_trapix` takes its committed volume; a newly-absorbed trap has a
`cur_amounts` entry generally computed some time in the past, so it is projected forward to
`cur_time` under its saved inflow.

A non-full trap is a transitory frontier (`spill_path == 0`) that must hold `V < capacity`, so
its seed is clamped to `prevfloat(C)`.  When several traps reach `C` in the same instant only
one fires per event; the rest stay transitory until their own `:fill` fires on the next predict
(mass preserved to a ULP).

Nothing is mutated — the closure only reads `driver.vol_by_trapix`.
"""
function _driver_state0(driver::NetworkDriver, tstruct, rateinfo, cur_amounts, full_set,
                        cur_time, z_vol_tables)
    project(g) = first(fill_trap_until(g, rateinfo, cur_amounts[g], cur_time,
                                       tstruct, z_vol_tables, use_saved = true))
    return function (g)
        g in full_set && return _own_capacity(tstruct, g)
        v = haskey(driver.vol_by_trapix, g) ? driver.vol_by_trapix[g] : project(g)
        return min(v, prevfloat(_own_capacity(tstruct, g)))
    end
end

"""
    _rebuild_contexts!(driver, tstruct, rateinfo, infiltration, z_vol_tables, cur_amounts,
                       full_set, cur_time, endtime) -> driver

Discard `driver.contexts` and rebuild one context per current component, each seeded from the
committed volumes and re-predicted out to `endtime`.  A context is a pure function of its
component plus those volumes, so a full rebuild is correct after any structural change.

Mutates `driver.contexts` (replaced wholesale); the committed volumes and NBS layers it reads
are left as they are.
"""
function _rebuild_contexts!(driver::NetworkDriver, tstruct, rateinfo, infiltration,
                            z_vol_tables, cur_amounts, full_set, cur_time, endtime)

    state0 = _driver_state0(driver, tstruct, rateinfo, cur_amounts, full_set,
                            cur_time, z_vol_tables)
    driver.contexts = DynNetworkContext[]
    for net in driver.comps
        ctx = _make_context(net, tstruct, rateinfo, state0, cur_time,
                            driver.precip, driver.nbs_state)
        _predict_network!(ctx, tstruct, infiltration, z_vol_tables, endtime)
        push!(driver.contexts, ctx)
    end
    return driver
end

"""
    _driver_next_event(driver) -> (ctx_index, next_event) or nothing

The earliest predicted event across all components, as the index of the owning context and its
`(; time, trap, kind)` prediction.  `nothing` when every component predicts `:none`.  Reads the
cached predictions; nothing is mutated or re-solved.
"""
function _driver_next_event(driver::NetworkDriver)
    best_i = 0
    best_t = Inf
    for (i, ctx) in enumerate(driver.contexts)
        ev = ctx.next_event
        ev.kind == :none && continue
        if ev.time < best_t
            best_t = ev.time
            best_i = i
        end
    end
    best_i == 0 && return nothing
    return (best_i, driver.contexts[best_i].next_event)
end

"""
    _harvest!(driver) -> driver

Read every committed context back into the driver's persistent stores — node volumes into
`vol_by_trapix`, NBS layers into `nbs_state` — so they survive the structural mutation that
follows and can reseed the rebuilt contexts.

Mutates `driver.vol_by_trapix` and `driver.nbs_state`; the contexts are read only.  Call only
once the contexts are committed to the intended time, as it takes their state at face value.
"""
function _harvest!(driver::NetworkDriver)
    for ctx in driver.contexts
        for (i, g) in enumerate(ctx.global_ix)
            driver.vol_by_trapix[g] = ctx.state[i]
        end
        _store_nbs_state!(driver.nbs_state, ctx)
    end
    return driver
end

# ----------------------------------------------------------------------------

"""
    step_network_driver!(driver, tstruct, rateinfo, infiltration, z_vol_tables,
                         filled_traps, cur_amounts, endtime)
        -> (; time, trap, kind, migrated) or nothing

Advance the networks to their earliest predicted event: commit every context to that time,
apply the event's structural transition to the live components, hand the committed state off,
and rebuild and re-predict the contexts.

Mutates `driver` (components, contexts, committed volumes, NBS layers, `last_time`) and
`filled_traps`, which is flipped to match the fired event — an `:empty` parent also clears its
exposed children.  `cur_amounts` is read only; it only seeds newly-absorbed traps.

Returns the fired event, where `migrated` lists the trap indices that entered or left the
networks, or `nothing` when no event falls before `endtime` (the period boundary is
[`finalize_network_driver!`](@ref)'s job).
"""
function step_network_driver!(driver::NetworkDriver, tstruct, rateinfo, infiltration,
                              z_vol_tables, filled_traps, cur_amounts, endtime)
    sel = _driver_next_event(driver)
    sel === nothing && return nothing
    _, ev = sel
    T = ev.time
    T >= endtime && return nothing        # beyond the period boundary — finalize handles it

    # commit EVERY context to T under its cached inflow before the structural change; T is the
    # earliest event, so the advance is event-free and order-independent for the others
    for ctx in driver.contexts
        _commit_network!(ctx, tstruct, infiltration, z_vol_tables, T)
    end
    driver.last_time = T
    _harvest!(driver)

    # apply the fired transition: update filled_traps, clamp committed volumes, mutate topology
    trap     = ev.trap
    migrated = Int[]
    if ev.kind == :fill
        filled_traps[trap] = true
        driver.vol_by_trapix[trap] = _own_capacity(tstruct, trap)         # pin to capacity
        migrated = apply_fill!(driver.comps, tstruct, findall(filled_traps), trap)
    elseif ev.kind == :unspill
        filled_traps[trap] = false
        driver.vol_by_trapix[trap] = prevfloat(_own_capacity(tstruct, trap))
        migrated = apply_unfill!(driver.comps, tstruct, trap)
    elseif ev.kind == :empty
        # de-subsumption: parent drains to its floor, exposing its children.  They MUST leave
        # full_traps before the regrow, or it re-subsumes them at V=0 and :empty re-fires forever.
        for c in subtrapsof(tstruct, trap)
            filled_traps[c] = false
            driver.vol_by_trapix[c] = prevfloat(_own_capacity(tstruct, c))
        end
        filled_traps[trap] = false
        driver.vol_by_trapix[trap] = 0.0
        migrated = apply_empty!(driver.comps, tstruct, findall(filled_traps), trap)
    end

    # rebuild the contexts off the mutated component set (coverage/inflow read fresh)
    _rebuild_contexts!(driver, tstruct, rateinfo, infiltration, z_vol_tables,
                       cur_amounts, Set(findall(filled_traps)), T, endtime)
    return (; time = T, trap = trap, kind = ev.kind, migrated = migrated)
end

# ----------------------------------------------------------------------------

"""
    finalize_network_driver!(driver, cur_amounts, tstruct, infiltration, z_vol_tables, endtime)
        -> cur_amounts

Period-boundary finalize: advance each context to `endtime` and settle the node volumes (and
subsumed capacities) into `cur_amounts`, carrying NBS layer storage into the next period.

Mutates `cur_amounts` (returned) and `driver` (contexts committed to `endtime`, `nbs_state`
updated).
"""
finalize_network_driver!(driver::NetworkDriver, cur_amounts, tstruct, infiltration,
                         z_vol_tables, endtime) =
    _finalize_networks!(cur_amounts, driver.contexts, tstruct, infiltration,
                        z_vol_tables, driver.last_time, endtime, driver.nbs_state)

# ----------------------------------------------------------------------------

"""
    _touch_networks_driver!(driver, changetimeest, sgraph, tstruct, filled_traps, cur_amounts,
                            rateinfo, z_vol_tables, infiltration, fill_updates, old_covered,
                            cur_time, endtime) -> (net_trap_set, net_covered_set)

The `fill_sequence` touch entry, called once per event that affects a network: commit every
context to `cur_time` under its cached inflow, apply the fired transitions via `apply_fill!` /
`apply_unfill!` / `apply_empty!`, reconcile the spillgraph to the new coverage, then rebuild
and re-predict the contexts.

`filled_traps` is ALREADY updated by the caller and, like `cur_amounts`, `fill_updates` and
`old_covered` (the covered set before this event), is read only here.  Mutates `driver`
(components, contexts, committed volumes, NBS layers, `last_time`), `changetimeest` (the
network traps get exact estimates), and — only when coverage actually changed — `sgraph` and
`rateinfo`.

Returns the refreshed `(net_trap_set, net_covered_set)`, both empty when the driver holds no
contexts; the committed volumes live in `driver.vol_by_trapix`.
"""
function _touch_networks_driver!(driver::NetworkDriver, changetimeest, sgraph, tstruct,
                                 filled_traps, cur_amounts, rateinfo, z_vol_tables,
                                 infiltration, fill_updates, old_covered, cur_time, endtime)
    isempty(driver.contexts) && return Set{Int}(), Set{Int}()

    # which node traps fired, and how (captured before the mutation below)
    fired = _capture_fired_kinds(driver.contexts, fill_updates)

    # commit every context to cur_time under cached inflow, harvest, then overlay the fired-trap
    # boundary values the solver expects
    for ctx in driver.contexts
        _commit_network!(ctx, tstruct, infiltration, z_vol_tables, cur_time)
    end
    driver.last_time = cur_time
    _harvest!(driver)
    _apply_fired_boundaries!(driver.vol_by_trapix, fired, tstruct)

    # apply each fired structural transition to the live components
    full_traps = findall(filled_traps)
    for (ft, k) in fired
        k == :fill    ? apply_fill!(driver.comps, tstruct, full_traps, ft) :
        k == :unspill ? apply_unfill!(driver.comps, tstruct, ft)           :
        k == :empty   ? apply_empty!(driver.comps, tstruct, full_traps, ft) : nothing
    end

    # reconcile the spillgraph/rateinfo to the new coverage before the rebuild reads inflow
    new_covered = _covered_of(driver.comps, tstruct)
    old_covered != new_covered &&
        _reconcile_spillgraph!(sgraph, rateinfo, filled_traps, new_covered, tstruct)

    # rebuild the contexts off the mutated component set and re-predict
    _rebuild_contexts!(driver, tstruct, rateinfo, infiltration, z_vol_tables,
                       cur_amounts, Set(full_traps), cur_time, endtime)
    nts = _net_trap_set(driver.contexts)
    ncs = _net_covered_set(driver.contexts, tstruct)
    _apply_network_changetimeest!(changetimeest, driver.contexts, ncs)
    return nts, ncs
end
