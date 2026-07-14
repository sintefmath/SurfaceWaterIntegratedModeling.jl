# ----------------------------------------------------------------------------
# The incremental network driver: builds the components once per weather period, then mutates
# the live `Vector{DynNetwork}` in place per event (`apply_fill!`/`apply_unfill!`/`apply_empty!`)
# instead of retracing.  `vol_by_trapix` (committed node volumes) and `nbs_state` (NBS layers)
# are the source of truth surviving structural change; contexts are rebuilt from them.
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
# Build the driver for one weather period: trace the components, seed committed volumes from
# `cur_amounts`, and build + predict one context per component.
function build_network_driver(tstruct, dyn_coords, culverts, full_traps, cur_amounts,
                              rateinfo, infiltration, z_vol_tables, cur_time, endtime;
                              nbs_placements = DynNBSPlacement[],
                              nbs_state      = Dict{Int,Vector{Float64}}())
    comps = setup_network(tstruct, full_traps;
                          dyn_coords = dyn_coords, culverts = culverts, nbs = nbs_placements)

    # committed volume of every node from the caller's amounts (all current at build, no projection)
    vol = Dict{Int,Float64}()
    for net in comps, t in net.traps
        vol[t.trap_ix] = cur_amounts[t.trap_ix].amount
    end

    driver = NetworkDriver(comps, DynNetworkContext[], vol, nbs_state, cur_time)
    _rebuild_contexts!(driver, tstruct, rateinfo, infiltration, z_vol_tables,
                       cur_amounts, Set(full_traps), cur_time, endtime)
    return driver
end

# Returns a function giving the current water volume for a given trap `g`.
#   * a full trap sits at its exact own capacity `C`;
#   * a node already in `vol_by_trapix` takes its already computed volume;
#   * a newly-absorbed trap has a `cur_amounts` entry that generally was computed at some time
#     in the past, so project it forward to `cur_time` under the saved inflow.
# NB: A non-full trap is a transitory frontier (`spill_path == 0`) that must hold `V < capacity`,
# so its seed is clamped to `prevfloat(C)`: when several traps reach `C` in the same instant,
# only one fires per event and the rest stay transitory until their `:fill` fires on the next
# predict (mass preserved to a ULP).
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

# Rebuild every context from the current `comps` and committed volumes and re-predict.  A
# context is a pure function of its component + committed volumes, so a full rebuild is correct.
function _rebuild_contexts!(driver::NetworkDriver, tstruct, rateinfo, infiltration,
                            z_vol_tables, cur_amounts, full_set, cur_time, endtime)

    state0 = _driver_state0(driver, tstruct, rateinfo, cur_amounts, full_set,
                            cur_time, z_vol_tables)
    driver.contexts = DynNetworkContext[]
    for net in driver.comps
        ctx = _make_context(net, tstruct, rateinfo, state0, cur_time, driver.nbs_state)
        _predict_network!(ctx, tstruct, infiltration, z_vol_tables, endtime)
        push!(driver.contexts, ctx)
    end
    return driver
end

# The earliest predicted event across all components as `(ctx_index, next_event)`, or `nothing`
# when every component predicts `:none`.
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

# Read every committed context back into the persistent stores (node volumes -> `vol_by_trapix`,
# NBS layers -> `nbs_state`), so they survive the structural mutation that follows.
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
# Advance to the earliest predicted event, apply its structural transition to the live
# components, hand off state, and rebuild the contexts.  Mutates `filled_traps` to match the
# fired event; `cur_amounts` seeds any newly-absorbed trap.  Returns the fired event
# `(; time, trap, kind, migrated)`, or `nothing` if none falls before `endtime`.
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
# Period-boundary finalize: advance each context to `endtime` and settle node volumes (and
# subsumed capacities) into `cur_amounts`, carrying NBS layer storage into the next period.
finalize_network_driver!(driver::NetworkDriver, cur_amounts, tstruct, infiltration,
                         z_vol_tables, endtime) =
    _finalize_networks!(cur_amounts, driver.contexts, tstruct, infiltration,
                        z_vol_tables, driver.last_time, endtime, driver.nbs_state)

# ----------------------------------------------------------------------------
# The `fill_sequence` touch entry: commit under cached inflow, apply the fired transitions via
# `apply_*`, reconcile the spillgraph to the new coverage, then rebuild + re-predict the
# contexts.  `filled_traps` is ALREADY updated by the caller.  Returns the refreshed
# `(net_trap_set, net_covered_set)`; committed volumes live in `driver.vol_by_trapix`.
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
