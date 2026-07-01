# ============================================================================
# DynNetworkContext: drives one DynNetwork between topology changes inside
# `fill_sequence`.  This is the integration layer between the analytical
# fill-sequence machinery and the multi-trap ODE solver `solveDynNetwork!`.
#
# See agent/reports/integrate_networks_plan.md (§2 state, §3 lifecycle, §4
# external inflow, §5 bounded integration).  This file is slice 2 of the
# implementation order: the context type and its build / predict / commit
# helpers.  Wiring into the fill_sequence event loop is a later slice.
# ============================================================================

"""
    DynNetworkContext

Mutable state for running and maintaining one [`DynNetwork`](@ref) between topology
changes during a `fill_sequence` run.

Fields:
- `net`: the current network (rebuilt when its topology changes).
- `state`: committed trap volumes (net of subtraps) at `last_solve_time`, indexed
  as `net.traps`.  Never advanced past the committed time.
- `global_ix`: the `TrapStructure` trap index of each `net.traps` entry.
- `inflow_sources`: global trap indices whose `trap_inflow` feeds this network's
  external inflow (the leaf descendants summed by `_external_inflow`); used by the
  touch test in the event loop.
- `last_solve_time`: absolute time the committed `state` refers to.
- `extern_inflow`: per-node external inflow the `state` is currently evolving
  under (cached at the last touch; the rate the commit/predict solves use).
- `next_event`: cached prediction `(; time, trap, kind)` of the network's next
  event, with `time` an ABSOLUTE timestamp and `trap` a global trap index.
"""
mutable struct DynNetworkContext
    net             ::DynNetwork
    state           ::Vector{Float64}
    global_ix       ::Vector{Int}
    inflow_sources  ::Set{Int}
    seeds           ::Vector{CartesianIndex{2}}  # the seed cells that produced this
                                                 # component (for per-context rebuild)
    last_solve_time ::Float64
    extern_inflow   ::Vector{Float64}
    next_event      ::NamedTuple
end

# ----------------------------------------------------------------------------
# Per-node gross composite external inflow (plan §4): for each node, the sum of
# `trap_inflow` over its leaf descendants (= `trap_inflow[node]` for a leaf node).
# `lowest_subtraps_for[t]` lists exactly those leaf descendants (== [t] for a
# leaf).  The leaf-sum is used rather than `getinflow(rateinfo, t.trap_ix)`
# directly: `_reconcile_spillgraph!` removes covered children from the spillgraph
# and withdraws their flow from `trap_inflow[parent]`, so that value is stale by
# the time `_external_inflow` is called.  Leaf `trap_inflow` values are unaffected
# and always current.
function _external_inflow(net::DynNetwork, rateinfo, tstruct)
    return Float64[sum(getinflow(rateinfo, leaf)
                       for leaf in tstruct.lowest_subtraps_for[t.trap_ix])
                   for t in net.traps]
end

# Global trap indices whose `trap_inflow` feeds this network's external inflow:
# the union over all nodes of their leaf descendants.  Used for the touch test
# (a network is touched when one of these appears in `getinflowupdates`).
function _inflow_sources(net::DynNetwork, tstruct)
    s = Set{Int}()
    for t in net.traps
        union!(s, tstruct.lowest_subtraps_for[t.trap_ix])
    end
    return s
end

# All (transitive) sub-traps of `t`, walking the agglomeration hierarchy downward.
function _descendants(tstruct, t::Int)
    out   = Int[]
    stack = collect(subtrapsof(tstruct, t))
    while !isempty(stack)
        x = pop!(stack)
        push!(out, x)
        append!(stack, subtrapsof(tstruct, x))
    end
    return out
end

# The trap indices that are *nodes* of some active network.
_net_trap_set(contexts) = Set{Int}(g for ctx in contexts for g in ctx.global_ix)

# `net_trap_set` plus every descendant subsumed under a network parent node
# (Design A): those subsumed children are full and static while their parent is a
# node, so they must be excluded from the standard changetime machinery.
function _net_covered_set(contexts, tstruct)
    s = Set{Int}()
    for ctx in contexts, g in ctx.global_ix
        push!(s, g)
        union!(s, _descendants(tstruct, g))
    end
    return s
end

# ----------------------------------------------------------------------------
# Predict the network's next event on a COPY of the committed state (the committed
# `state` is left untouched), bounded by the time remaining to `endtime`.  Solves
# under the cached `ctx.extern_inflow` (the rate the state is evolving under) and
# stores an ABSOLUTE-time `next_event`.  Returns that `next_event`.
function _predict_network!(ctx::DynNetworkContext, tstruct, infiltration,
                           z_vol_tables, endtime)
    res = solveDynNetwork!(copy(ctx.state), tstruct, ctx.net, infiltration,
                           ctx.extern_inflow;
                           tmax = endtime - ctx.last_solve_time, zvt = z_vol_tables)
    ctx.next_event = (; time = res.time + ctx.last_solve_time,
                        trap = res.trap, kind = res.kind)
    return ctx.next_event
end

# Advance the committed `state` IN PLACE to absolute time `T_commit`, under the
# cached `extern_inflow` (the rates valid since the last touch — so this is
# order-independent of any intervening spillgraph/flow update).  `state` is left at
# the exact `T_commit` volumes and `last_solve_time` is updated.  Returns the solver
# result (so the caller can act on an event that coincides with `T_commit`).
function _commit_network!(ctx::DynNetworkContext, tstruct, infiltration,
                          z_vol_tables, T_commit)
    res = solveDynNetwork!(ctx.state, tstruct, ctx.net, infiltration,
                           ctx.extern_inflow;
                           tmax = T_commit - ctx.last_solve_time, zvt = z_vol_tables)
    ctx.last_solve_time = T_commit
    return res
end

# ----------------------------------------------------------------------------
"""
    _build_dyn_networks(tstruct, dyn_traps, culverts, filled_traps, cur_amounts,
                        rateinfo, infiltration, z_vol_tables, cur_time, endtime)
        -> (contexts, net_trap_set, net_covered_set)

Construct the [`DynNetworkContext`](@ref)s for one weather period.  Networks are
seeded by the union of the explicit `dyn_traps` (trap indices, each converted to a
representative footprint cell) and every culvert's inlet/outlet cells; `setup_network`
then builds the disjoint network components, and each becomes a context whose
initial `state` is read from `cur_amounts`.  Each context's first event is predicted
(stored in `next_event`).

`dyn_traps` must be valid trap indices (the only seed-side validation needed, since
trap-indexed seeds are always inside a trap — no "hanging" paths).
"""
# Seed cells for the dynamic networks: a representative footprint cell of each
# named `dyn_trap` (validated) plus every culvert's inlet/outlet cell.
function _dyn_seeds(tstruct, dyn_traps, culverts)
    CI     = CartesianIndices(size(tstruct.topography))
    ntraps = numtraps(tstruct)
    for t in dyn_traps
        1 <= t <= ntraps ||
            error("_dyn_seeds: dyn_trap $t is not a valid trap index (1:$ntraps)")
    end
    seeds = CartesianIndex{2}[CI[tstruct.footprints[t][1]] for t in dyn_traps]
    for cv in culverts
        push!(seeds, cv.inlet)
        push!(seeds, cv.outlet)
    end
    return unique!(seeds)
end

function _build_dyn_networks(tstruct, dyn_traps, culverts, full_traps, cur_amounts,
                             rateinfo, infiltration, z_vol_tables, cur_time, endtime)
    seeds = _dyn_seeds(tstruct, dyn_traps, culverts)
    isempty(seeds) && return (DynNetworkContext[], Set{Int}(), Set{Int}())

    components = setup_network(tstruct, seeds, full_traps; culverts=culverts)

    contexts = DynNetworkContext[]
    for net in components
        ctx = _make_context(net, tstruct, rateinfo, seeds,
                            g -> cur_amounts[g].amount, cur_time)
        _predict_network!(ctx, tstruct, infiltration, z_vol_tables, endtime)
        push!(contexts, ctx)
    end

    return contexts, _net_trap_set(contexts), _net_covered_set(contexts, tstruct)
end

# Build one context from a single network component.  `seed_pool` is the full seed
# list; the context records the subset that falls inside this component (used to
# re-trace the same component on a rebuild).  `state0(g)` supplies the initial
# committed volume for the global trap index `g`.
function _make_context(net::DynNetwork, tstruct, rateinfo, seed_pool, state0, cur_time)
    global_ix = Int[t.trap_ix for t in net.traps]
    occ       = _occupied_cells(tstruct, [net])
    seeds     = CartesianIndex{2}[s for s in seed_pool if s in occ]
    state     = Float64[state0(g) for g in global_ix]
    return DynNetworkContext(net, state, global_ix,
                             _inflow_sources(net, tstruct), seeds,
                             cur_time,
                             _external_inflow(net, rateinfo, tstruct),
                             (; time = Inf, trap = 0, kind = :none))
end

# Own storage capacity of a trap (net of subtraps) — the same quantity the solver
# uses as `geom.capacity` and `fill_sequence` records for a full trap.
_own_capacity(tstruct, t::Int) = tstruct.trapvolumes[t] - tstruct.subvolumes[t]

# Pin every full trap of `ctx` to its exact capacity, removing the small ODE drift a
# commit leaves on FULL pass-through traps (the solver's exact `V == C` check would
# otherwise misclassify them as TRANSITORY).  `full_set` is the authoritative set of
# full global trap indices.
function _clamp_full_traps!(ctx::DynNetworkContext, tstruct, full_set)
    for (i, g) in enumerate(ctx.global_ix)
        g in full_set && (ctx.state[i] = _own_capacity(tstruct, g))
    end
end

# Write each context's predicted next event into `changetimeest` as an EXACT
# (min == max) estimate, and force every covered trap that is not a predicted
# trigger to `Inf`.  Used by both `_set_initial_changetime_estimates` and the touch
# step so the branch-and-bound resolves network traps from the network prediction.
function _apply_network_changetimeest!(changetimeest, net_contexts, net_covered_set)
    for t in net_covered_set
        changetimeest[t] = ChangeTimeEstimate(false, Inf, Inf)
    end
    for ctx in net_contexts
        ev = ctx.next_event
        ev.kind == :none && continue
        changetimeest[ev.trap] = ChangeTimeEstimate(ev.kind == :fill, ev.time, ev.time)
    end
    return changetimeest
end

# ----------------------------------------------------------------------------
# Reconcile the spillgraph (and the flow it drives in `rateinfo`) to the invariant
# "the spillgraph contains exactly the FULL ∧ ¬COVERED traps".  A trap that filled
# while non-network deposited its outflow downstream via the spillgraph; once it is
# absorbed into a network its outflow is routed by the ODE instead, so that deposit
# must be withdrawn (and vice-versa when a full trap LEAVES coverage, e.g. draining).
# Recomputes the target spillgraph from the effective filled state (cheap — no
# `compute_flow`), diffs it against the current one, and propagates only the changed
# edges' flow.  Design A's subsumption keeps the effective state valid (a full
# non-covered parent never has a covered child).
function _reconcile_spillgraph!(sgraph, rateinfo, filled_traps, covered, tstruct)
    n         = numtraps(tstruct)
    effective = Bool[filled_traps[t] && !(t in covered) for t in 1:n]
    new_sg    = compute_complete_spillgraph(tstruct, effective)
    graph_updates = IncrementalUpdate{Tuple{Int,Int}}[]
    for t in union(keys(sgraph.edges), keys(new_sg.edges))
        old_t = get(sgraph.edges, t, 0)
        new_t = get(new_sg.edges, t, 0)
        old_t != new_t && push!(graph_updates, IncrementalUpdate(t, (old_t, new_t)))
    end
    sgraph.edges = new_sg.edges
    isempty(graph_updates) || _update_flow!(rateinfo, graph_updates, tstruct, sgraph)
    return graph_updates
end

# ----------------------------------------------------------------------------
# Advance every network to `cur_time` and rebuild the whole set after a status
# change, per plan §3.  Each context is committed to `cur_time` under its cached
# external inflow (full pass-through traps pinned back to exact capacity), then the
# networks are rebuilt from the global seed pool — letting `setup_network` merge or
# split components correctly as traps fill — with committed volumes remapped by
# global trap index.  Must run AFTER `_update_flow!` (it reads the refreshed
# `rateinfo`).  Returns the new context vector and refreshed `net_trap_set` /
# `net_covered_set`.
#
# @@@ This commits and re-predicts ALL contexts on every event.  The touch-gating
# optimisation (skip networks whose inflow did not change, via `getinflowupdates`)
# is deferred — it needs context-merge handling to remain correct when growing
# networks overlap.  See agent/reports/integrate_networks_plan.md §3/§8.
function _touch_networks!(net_contexts, changetimeest, sgraph, tstruct, dyn_traps, culverts,
                          filled_traps, cur_amounts, rateinfo, z_vol_tables, infiltration,
                          fill_updates, old_covered, cur_time, endtime)
    isempty(net_contexts) &&
        return net_contexts, Set{Int}(), Set{Int}(), Dict{Int,Float64}()
    full_traps = findall(filled_traps)
    full_set   = Set(full_traps)

    # Capture the predicted kind of each fired trap before the contexts are mutated.
    fired_kind = Dict{Int,Symbol}()
    for ctx in net_contexts
        ev = ctx.next_event
        ev.kind != :none && any(u.index == ev.trap for u in fill_updates) &&
            (fired_kind[ev.trap] = ev.kind)
    end

    # 1. Commit every context to cur_time, pin full traps to exact capacity, and
    #    collect the committed volumes by global trap index.
    committed = Dict{Int,Float64}()
    for ctx in net_contexts
        _commit_network!(ctx, tstruct, infiltration, z_vol_tables, cur_time)
        _clamp_full_traps!(ctx, tstruct, full_set)
        for (i, g) in enumerate(ctx.global_ix)
            committed[g] = ctx.state[i]
        end
    end
    # Boundary values for fired traps that did not just become full.  An :empty
    # parent drops to 0 and EXPOSES its immediate children — they go from full to
    # transitory (just below their own capacity), so they leave `full_traps` (the
    # caller already flipped them — §`_expand_empty_fill_updates`) and start at
    # prevfloat(C_child).  An :unspill trap drops to prevfloat(C).
    for (ft, k) in fired_kind
        if k == :empty
            committed[ft] = 0.0
            for c in subtrapsof(tstruct, ft)
                committed[c] = prevfloat(_own_capacity(tstruct, c))
            end
        elseif k == :unspill
            committed[ft] = prevfloat(_own_capacity(tstruct, ft))
        end
    end

    # 2. Rebuild the network STRUCTURE from the global seeds and compute the new
    #    coverage.
    seeds      = _dyn_seeds(tstruct, dyn_traps, culverts)
    components = setup_network(tstruct, seeds, full_traps; culverts=culverts)
    new_covered = Set{Int}()
    for net in components, t in net.traps
        push!(new_covered, t.trap_ix)
        union!(new_covered, _descendants(tstruct, t.trap_ix))
    end

    # 3. Reconcile the spillgraph + rateinfo to the new coverage BEFORE the network
    #    external inflow is read — so absorbed traps' double-counted deposits are
    #    withdrawn (and traps leaving coverage regain their edges).
    old_covered != new_covered &&
        _reconcile_spillgraph!(sgraph, rateinfo, filled_traps, new_covered, tstruct)

    # 4. Build the contexts (external inflow now reads the reconciled rateinfo) and
    #    predict each.
    # For traps newly absorbed into the network (not in committed and not full), the
    # cur_amounts entry is stale — it holds the volume at the last committed event, not
    # at cur_time.  Project it forward to cur_time using the pre-reconcile (saved) inflow,
    # the same inflow that was active while the plain path was accumulating water there.
    function _state0_project(g)
        vol, _ = fill_trap_until(g, rateinfo, cur_amounts[g], cur_time,
                                 tstruct, z_vol_tables, use_saved=true)
        return vol
    end
    state0(g) = g in full_set        ? _own_capacity(tstruct, g) :
                haskey(committed, g) ? committed[g]              :
                                       _state0_project(g)
    new_contexts = DynNetworkContext[]
    for net in components
        c = _make_context(net, tstruct, rateinfo, seeds, state0, cur_time)
        _predict_network!(c, tstruct, infiltration, z_vol_tables, endtime)
        push!(new_contexts, c)
    end

    nts = _net_trap_set(new_contexts)
    ncs = _net_covered_set(new_contexts, tstruct)
    _apply_network_changetimeest!(changetimeest, new_contexts, ncs)
    return new_contexts, nts, ncs, committed
end

# ----------------------------------------------------------------------------
# Weather-period boundary finalization (plan §10).  A network trap follows the
# multi-trap ODE and cannot be projected to `endtime` with the constant-rate
# `fill_trap_until`, so each context is advanced to `endtime` under its cached
# external inflow and the network traps' amounts are read from the settled state.
# The advance is event-free by construction — any network event ≤ `endtime` was
# already processed before the event loop broke — so it settles cleanly.  These
# exact boundary volumes are what the NEXT weather period rebuilds its networks
# from, so they must not use the constant-rate projection.
#
# Runs for EVERY context, including ones whose state still sits at an earlier
# `last_solve_time` (the commit's `tmax` guard makes advancing an already-current
# context a no-op).  A node takes its settled ODE volume (net of subtraps); a
# subsumed full descendant sits at its own capacity — matching `_network_amount_updates`.
function _finalize_networks!(cur_amounts, net_contexts, tstruct, infiltration,
                             z_vol_tables, cur_time, endtime)
    isempty(net_contexts) && return cur_amounts
    stamp = min(cur_time, endtime)
    for ctx in net_contexts
        ctx.last_solve_time < endtime &&
            _commit_network!(ctx, tstruct, infiltration, z_vol_tables, endtime)
        for (i, g) in enumerate(ctx.global_ix)
            cur_amounts[g] = FilledAmount(ctx.state[i], stamp)
            for d in _descendants(tstruct, g)   # subsumed full descendants
                cur_amounts[d] = FilledAmount(_own_capacity(tstruct, d), stamp)
            end
        end
    end
    return cur_amounts
end

# ----------------------------------------------------------------------------
# A network *parent* that fired :empty exposes its immediate children as transitory
# (draining) traps.  Without this they stay in `filled_traps`, the rebuild re-subsumes
# them under the parent at V == 0, and the parent re-fires :empty at the same instant
# forever.  Returns `fill_updates` augmented with a `(child, false)` entry per exposed
# child (the caller then flips `filled_traps` for them like any other update).
function _expand_empty_fill_updates(fill_updates, net_contexts, tstruct)
    kind_of = Dict{Int,Symbol}()
    for ctx in net_contexts
        ev = ctx.next_event
        ev.kind != :none && (kind_of[ev.trap] = ev.kind)
    end
    extra = IncrementalUpdate{Bool}[]
    seen  = Set(u.index for u in fill_updates)
    for u in fill_updates
        (!u.value && get(kind_of, u.index, :none) == :empty) || continue
        for c in subtrapsof(tstruct, u.index)
            c in seen && continue
            push!(extra, IncrementalUpdate{Bool}(c, false))
            push!(seen, c)
        end
    end
    return isempty(extra) ? fill_updates : vcat(fill_updates, extra)
end

# Amount updates for network traps whose volume the network owns (§9): a node carries
# its committed ODE volume; a subsumed full descendant is at capacity; a trap that just
# LEFT the networks (e.g. an emptied parent) takes its committed boundary value.
# `affected` is the union of the pre- and post-touch covered sets.  Only traps whose
# amount actually changed (vs `cur_amounts`) are emitted, keeping the event incremental.
function _network_amount_updates(net_contexts, affected, committed, tstruct,
                                 cur_amounts, cur_time)
    node_amt = Dict{Int,Float64}()
    for ctx in net_contexts, (i, g) in enumerate(ctx.global_ix)
        node_amt[g] = ctx.state[i]
    end
    covered = _net_covered_set(net_contexts, tstruct)
    updates = IncrementalUpdate{FilledAmount}[]
    for t in affected
        amt = haskey(node_amt, t)  ? node_amt[t]                :   # a current node
              t in covered         ? _own_capacity(tstruct, t)  :   # subsumed (full)
              get(committed, t, _own_capacity(tstruct, t))           # just exited
        amt == cur_amounts[t].amount && continue
        push!(updates, IncrementalUpdate(t, FilledAmount(amt, cur_time)))
    end
    return updates
end
