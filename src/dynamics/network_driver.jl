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
# structural mutation via `_make_context` (whose `state0` closure is `_state_for`'s
# per-trap rule).  See agent/GATE_INTEGRATION_PLAN.md §2 / §7.
# ----------------------------------------------------------------------------

mutable struct NetworkDriver
    comps          ::Vector{DynNetwork}          # live components, mutated by apply_*
    contexts       ::Vector{DynNetworkContext}   # one per comp, same order (rebuilt each event)
    vol_by_trapix  ::Dict{Int,Float64}           # node trap_ix -> committed volume at last_time
    nbs_state      ::Dict{Int,Vector{Float64}}   # placement id -> layer states (persistent)
    nbs_placements ::Vector{DynNBSPlacement}
    culverts       ::Vector{DynCulvert}
    dyn_coords     ::Vector{CartesianIndex{2}}   # explicit dynamic seed cells (period-constant)
    last_time      ::Float64                      # absolute time the committed volumes hold at
end

# ----------------------------------------------------------------------------
# C1 — re-indexing helper.  Build the trap-volume portion of a component's solver
# `state`, in `net.traps` order, from the committed `vol_by_trapix` map.  A node not
# yet in the map (freshly grown / absorbed) is seeded by `seed0(trap_ix)`.  NBS layer
# states are appended separately by `_make_context` (keyed by the persistent store),
# so this returns only the leading `length(net.traps)` trap entries.
#   net           : the component whose state vector is being (re)built
#   vol_by_trapix : trap_ix -> committed volume
#   seed0         : trap_ix -> initial volume for a node absent from the map
# returns a Vector{Float64} of length length(net.traps).
_state_for(net::DynNetwork, vol_by_trapix, seed0) =
    Float64[get(() -> seed0(t.trap_ix), vol_by_trapix, t.trap_ix) for t in net.traps]

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

    driver = NetworkDriver(comps, DynNetworkContext[], vol, nbs_state,
                           nbs_placements, culverts, dyn_coords, cur_time)
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
function _driver_state0(driver::NetworkDriver, tstruct, rateinfo, cur_amounts, full_set,
                        cur_time, z_vol_tables)
    project(g) = first(fill_trap_until(g, rateinfo, cur_amounts[g], cur_time,
                                       tstruct, z_vol_tables, use_saved = true))
    return g -> g in full_set                    ? _own_capacity(tstruct, g) :
                haskey(driver.vol_by_trapix, g)  ? driver.vol_by_trapix[g]   :
                                                   project(g)
end

# Rebuild every context from the current `comps` and committed volumes, then predict each
# next event.  Called after build and after any structural mutation: a context is a pure
# function of its component + the committed volumes, so a full rebuild is correct (and
# cheap — one solve per component).  `full_set`/`cur_amounts` feed the `state0` seed rule.
function _rebuild_contexts!(driver::NetworkDriver, tstruct, rateinfo, infiltration,
                            z_vol_tables, cur_amounts, full_set, cur_time, endtime)
    state0 = _driver_state0(driver, tstruct, rateinfo, cur_amounts, full_set,
                            cur_time, z_vol_tables)
    driver.contexts = DynNetworkContext[]
    for net in driver.comps
        ctx = _make_context(net, tstruct, rateinfo, driver.dyn_coords, state0, cur_time,
                            driver.nbs_placements, driver.nbs_state)
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
