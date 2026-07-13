# ----------------------------------------------------------------------------
# Incremental network driver (gate Phase C).
#
# Replaces the structural core of `network_context.jl` (full retrace every event via
# `setup_network_cached` / `SubnetCache` / `_reuse_plan`) with the incremental
# membership layer: build once per weather period, then per event mutate the live
# `Vector{DynNetwork}` in place via `apply_fill!` / `apply_unfill!` / `apply_empty!`
# and re-sync only the touched components.
#
# The driver owns the committed per-node volumes (`vol_by_trapix`) and the persistent
# NBS layer store (`nbs_state`) as the source of truth that survives structural
# change; contexts (holding the solver `state`) are rebuilt from them after every
# structural mutation via `_make_context` (seeded by the `_driver_state0` per-trap
# rule).  See agent/GATE_INTEGRATION_PLAN.md §2 / §7.
# ----------------------------------------------------------------------------

mutable struct NetworkDriver
    comps          ::Vector{DynNetwork}          # live components, mutated by apply_* (the
                                                 # structural source of truth — culverts/NBS
                                                 # travel on each component, not the driver)
    contexts       ::Vector{DynNetworkContext}   # one per comp, same order (rebuilt each event)
    vol_by_trapix  ::Dict{Int,Float64}           # node trap_ix -> committed volume at last_time
    nbs_state      ::Dict{Int,Vector{Float64}}   # placement id -> layer states (persistent)
    last_time      ::Float64                      # absolute time the committed volumes hold at
end

# ----------------------------------------------------------------------------
# C2 — build entry.  Construct the driver for one weather period: build the components
# from the new `setup_network`, seed the committed volumes from `cur_amounts`, make one
# context per component, and predict each next event.
#   returns the NetworkDriver (its contexts carry the initial predictions).
function build_network_driver(tstruct, dyn_coords, culverts, full_traps, cur_amounts,
                              rateinfo, infiltration, z_vol_tables, cur_time, endtime;
                              nbs_placements = DynNBSPlacement[],
                              nbs_state      = Dict{Int,Vector{Float64}}())
    comps = setup_network(tstruct, full_traps;
                          dyn_coords = dyn_coords, culverts = culverts, nbs = nbs_placements)

    # committed volume of every node, read from the caller's per-trap amounts (all current
    # at build time, so no projection is needed here — see `_driver_state0`)
    vol = Dict{Int,Float64}()
    for net in comps, t in net.traps
        vol[t.trap_ix] = cur_amounts[t.trap_ix].amount
    end

    driver = NetworkDriver(comps, DynNetworkContext[], vol, nbs_state, cur_time)
    _rebuild_contexts!(driver, tstruct, rateinfo, infiltration, z_vol_tables,
                       cur_amounts, Set(full_traps), cur_time, endtime)
    return driver
end

# The committed-volume seed rule for a node `g` (C1's per-trap policy, shared by build
# and every post-event rebuild), following `_assemble_contexts`:
#   * a full trap sits at its exact own capacity `C`;
#   * a node already in `vol_by_trapix` takes its committed volume;
#   * a newly-absorbed trap (from the static path) has a STALE `cur_amounts` entry, so
#     project it forward to `cur_time` under the saved inflow ([[stale-state0-...]] §7.2).
# A NON-full trap is a transitory frontier (`spill_path == 0`), which the solver requires to
# hold `V < capacity`.  When many traps reach capacity in the same instant (heavy rain / full
# coverage) the b&b fires one and commits the others to `cur_time` at exactly `C`; seeding such
# a trap at `C` would read as full on a `spill_path == 0` node (three-state-contract error).
# Clamp a non-full seed to `prevfloat(C)`: it stays transitory and its `:fill` fires on the next
# predict, serialising the simultaneous fills (mass preserved to a ULP).
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

# Rebuild every context from the current `comps` and committed volumes, then predict each
# next event.  Called after build and after any structural mutation: a context is a pure
# function of its component + the committed volumes, so a full rebuild is correct (and
# cheap — one solve per component).  `full_set`/`cur_amounts` feed the `state0` seed rule.
function _rebuild_contexts!(driver::NetworkDriver, tstruct, rateinfo, infiltration,
                            z_vol_tables, cur_amounts, full_set, cur_time, endtime)

    # initial state rule for every trap in the current component set, used by `_make_context`
    state0 = _driver_state0(driver, tstruct, rateinfo, cur_amounts, full_set,
                            cur_time, z_vol_tables)
    driver.contexts = DynNetworkContext[]
    for net in driver.comps
        ctx = _make_context(net, tstruct, rateinfo, state0, cur_time, driver.nbs_state)

        # sets `ctx.next_event` to the earliest predicted event for this component, or :none
        _predict_network!(ctx, tstruct, infiltration, z_vol_tables, endtime)
        push!(driver.contexts, ctx)
    end
    return driver
end

# The earliest predicted event across all components, as `(ctx_index, event)` where
# `event` is the context's `next_event` NamedTuple (`time`, `trap`, `kind`).  Returns
# `nothing` when no component predicts a real event (all `:none`).
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

# Read every committed context back into the persistent stores: node volumes into
# `vol_by_trapix`, NBS layer states into `nbs_state`.  Called after committing all
# contexts to the event time, so the stores survive the structural mutation that follows.
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
# C3–C5 — per-event step.  Advance to the earliest predicted event, apply the matching
# structural transition to the live components, hand off state per the solver protocol,
# and rebuild the contexts on the mutated component set.  `filled_traps` is the caller's
# authoritative fill mask (mutated here to match the fired event); `cur_amounts` supplies
# the projection seed for any newly-absorbed trap.  Returns the fired event as
# `(; time, trap, kind, migrated)` (the migrated `trap_ix` set), or `nothing` when no
# event falls before `endtime`.
function step_network_driver!(driver::NetworkDriver, tstruct, rateinfo, infiltration,
                              z_vol_tables, filled_traps, cur_amounts, endtime)
    sel = _driver_next_event(driver)
    sel === nothing && return nothing
    _, ev = sel
    T = ev.time
    T >= endtime && return nothing        # beyond the period boundary — finalize handles it

    # (§7.3) commit EVERY context to T under its cached inflow BEFORE the structural
    # change; T is the earliest event, so this advance is event-free for the others and
    # order-independent.  Harvest the committed volumes/NBS states into the driver stores.
    for ctx in driver.contexts
        _commit_network!(ctx, tstruct, infiltration, z_vol_tables, T)
    end
    driver.last_time = T
    _harvest!(driver)

    # Apply the fired transition: update `filled_traps`, clamp the committed volumes per the
    # solver protocol, and mutate the component topology via the incremental membership layer.
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
        # De-subsumption: the parent drains to its floor, exposing its immediate children.
        # They MUST leave `full_traps` before the regrow (§7.2 hazard) or it re-subsumes
        # them at V=0 and :empty re-fires forever.
        for c in subtrapsof(tstruct, trap)
            filled_traps[c] = false
            driver.vol_by_trapix[c] = prevfloat(_own_capacity(tstruct, c))
        end
        filled_traps[trap] = false
        driver.vol_by_trapix[trap] = 0.0
        migrated = apply_empty!(driver.comps, tstruct, findall(filled_traps), trap)
    end

    # (§7.1) rebuild the contexts off the MUTATED component set (coverage/inflow read fresh,
    # not diffed); newly-grown traps are seeded via the projection rule in `_driver_state0`.
    _rebuild_contexts!(driver, tstruct, rateinfo, infiltration, z_vol_tables,
                       cur_amounts, Set(findall(filled_traps)), T, endtime)
    return (; time = T, trap = trap, kind = ev.kind, migrated = migrated)
end

# ----------------------------------------------------------------------------
# C6 — period-boundary finalize.  Advance each context to `endtime` and settle the node
# volumes (and subsumed-descendant capacities) back into `cur_amounts`, carrying NBS layer
# storage into the next period.  Reuses `_finalize_networks!` (its contract is exactly ours).
finalize_network_driver!(driver::NetworkDriver, cur_amounts, tstruct, infiltration,
                         z_vol_tables, endtime) =
    _finalize_networks!(cur_amounts, driver.contexts, tstruct, infiltration,
                        z_vol_tables, driver.last_time, endtime, driver.nbs_state)

# ----------------------------------------------------------------------------
# D1 — the `fill_sequence` touch entry: the driver-based replacement for `_touch_networks!`.
# Where the old touch RETRACES the structure every event (`setup_network_cached` + `_reuse_plan`
# + `_assemble_contexts`), this applies the fired transitions incrementally via `apply_*`.  The
# order mirrors the old touch so the state handoff matches the solver protocol: commit under the
# cached inflow, mutate topology, reconcile the spillgraph to the new coverage (C5), then rebuild
# the contexts reading the reconciled inflow.  `filled_traps` is ALREADY updated by the caller
# (fired traps flipped, `:empty` children exposed).  Returns `(contexts, net_trap_set,
# net_covered_set, committed)` — the same tuple the old `_touch_networks!` returns.
function _touch_networks_driver!(driver::NetworkDriver, changetimeest, sgraph, tstruct,
                                 filled_traps, cur_amounts, rateinfo, z_vol_tables,
                                 infiltration, fill_updates, old_covered, cur_time, endtime)
    isempty(driver.contexts) &&
        return driver.contexts, Set{Int}(), Set{Int}(), Dict{Int,Float64}()

    # which node traps fired this event, and how (captured from predictions before mutation)
    fired = _capture_fired_kinds(driver.contexts, fill_updates)

    # commit every context to cur_time under its cached inflow, harvest into the stores, then
    # overlay the fired-trap boundary values (:unspill→prevfloat(C); :empty→parent 0 / children
    # prevfloat(C); :fill→C is applied by the full-trap branch of the rebuild's seed rule)
    for ctx in driver.contexts
        _commit_network!(ctx, tstruct, infiltration, z_vol_tables, cur_time)
    end
    driver.last_time = cur_time
    _harvest!(driver)
    _apply_fired_boundaries!(driver.vol_by_trapix, fired, tstruct)
    committed = copy(driver.vol_by_trapix)

    # apply each fired structural transition to the live components
    full_traps = findall(filled_traps)
    for (ft, k) in fired
        k == :fill    ? apply_fill!(driver.comps, tstruct, full_traps, ft) :
        k == :unspill ? apply_unfill!(driver.comps, tstruct, ft)           :
        k == :empty   ? apply_empty!(driver.comps, tstruct, full_traps, ft) : nothing
    end

    # reconcile the spillgraph/rateinfo to the new coverage BEFORE the rebuild reads inflow (C5)
    new_covered = _covered_of(driver.comps, tstruct)
    old_covered != new_covered &&
        _reconcile_spillgraph!(sgraph, rateinfo, filled_traps, new_covered, tstruct)

    # rebuild the contexts off the mutated component set, predict fresh events, publish estimates
    _rebuild_contexts!(driver, tstruct, rateinfo, infiltration, z_vol_tables,
                       cur_amounts, Set(full_traps), cur_time, endtime)
    nts = _net_trap_set(driver.contexts)
    ncs = _net_covered_set(driver.contexts, tstruct)
    _apply_network_changetimeest!(changetimeest, driver.contexts, ncs)
    return driver.contexts, nts, ncs, committed
end
