# ============================================================================
# DynNetworkContext: drives one DynNetwork between topology changes inside
# `fill_sequence` — the integration layer between the analytical fill-sequence
# machinery and the multi-trap ODE solver `solveDynNetwork!`.
# Design notes: agent/reports/integrate_networks_plan.md (§2 state, §3 lifecycle,
# §4 external inflow, §5 bounded integration).
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
# Per-node external inflow (plan §4): each node's sum of `trap_inflow` over its leaf
# descendants (`lowest_subtraps_for[t]`, == [t] for a leaf).  Summed from the leaves,
# not `getinflow(node)` directly, because `_reconcile_spillgraph!` withdraws covered
# children's flow from the parent's `trap_inflow`, leaving it stale — leaf values stay current.
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
# The NBS network elements for this terrain: one `DynNBS` per placement, carrying
# the artificial trap index (looked up from the footprint's region) and the
# resolved per-layer outlets.  Empty when the structure has no NBS placements.
function _nbs_elements(tstruct)
    isempty(tstruct.nbs) && return DynNBS[]
    return DynNBS[DynNBS(tstruct.regions[p.footprint[1]], pi, p.outlets)
                  for (pi, p) in enumerate(tstruct.nbs)]
end

# Seed cells for the dynamic networks: a representative footprint cell of each
# named `dyn_trap` (validated) plus every culvert's inlet/outlet cell and every
# NBS trap's footprint + outlet cells (so each outlet's downstream is traced in).
function _dyn_seeds(tstruct, dyn_traps, culverts, nbs_objs::Vector{DynNBS}=DynNBS[])
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
    for nb in nbs_objs
        push!(seeds, CI[tstruct.footprints[nb.trap_ix][1]])
        for oc in nb.outlets
            push!(seeds, oc)
        end
    end
    return unique!(seeds)
end

function _build_dyn_networks(tstruct, dyn_traps, culverts, full_traps, cur_amounts,
                             rateinfo, infiltration, z_vol_tables, cur_time, endtime)
    nbs_objs = _nbs_elements(tstruct)
    seeds = _dyn_seeds(tstruct, dyn_traps, culverts, nbs_objs)
    isempty(seeds) && return (DynNetworkContext[], Set{Int}(), Set{Int}())

    # Returns a vector of DynNetwork
    components = setup_network(tstruct, seeds, full_traps; culverts=culverts, nbs=nbs_objs)

    contexts = DynNetworkContext[]
    for net in components
        ctx = _make_context(net, tstruct, rateinfo, seeds,
                            g -> cur_amounts[g].amount, cur_time)
        _predict_network!(ctx, tstruct, infiltration, z_vol_tables, endtime) # set ctx.next_event
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
    # Trap volumes, then one appended state per NBS layer (in `net.nbs` order, matching
    # `_build_nbs_plan`'s state_base offsets).  Layer states start empty; cross-period
    # persistence is added in 3c.
    nlayers   = sum(length(tstruct.nbs[nb.placement_ix].system.layers)
                    for nb in net.nbs; init = 0)
    state     = vcat(Float64[state0(g) for g in global_ix], zeros(Float64, nlayers))
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
# "spillgraph = exactly the FULL ∧ ¬COVERED traps": a trap absorbed into a network has
# its outflow routed by the ODE, not the spillgraph, so its deposit must be withdrawn
# (and restored when it leaves coverage).  Recomputes the target spillgraph from the
# effective filled state, diffs it, and propagates only the changed edges' flow.
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
# The old contexts that CHANGED this event: one whose member just fired, or whose
# external inflow was updated (an upstream plain trap spilled into it).  These must be
# re-committed and re-predicted; an unchanged context evolves under the same cached
# extern_inflow, so its committed state and absolute-time prediction stay exact and it
# can be carried over untouched.  Returns a set of indices into `net_contexts`.
function _affected_contexts(net_contexts, fill_updates, inflow_updated)
    return Set{Int}(ci for (ci, ctx) in enumerate(net_contexts)
                    if any(u.index ∈ ctx.global_ix for u in fill_updates) ||
                       !isdisjoint(ctx.inflow_sources, inflow_updated))
end

# The predicted kind (:fill/:empty/:unspill) of each trap that fired this event,
# captured from the contexts' predictions before anything is mutated.
function _capture_fired_kinds(net_contexts, fill_updates)
    fired_kind = Dict{Int,Symbol}()
    for ctx in net_contexts
        ev = ctx.next_event
        ev.kind != :none && any(u.index == ev.trap for u in fill_updates) &&
            (fired_kind[ev.trap] = ev.kind)
    end
    return fired_kind
end

# All traps covered by `components` (nodes plus their subsumed full descendants).
function _covered_of(components, tstruct)
    s = Set{Int}()
    for net in components, t in net.traps
        push!(s, t.trap_ix)
        union!(s, _descendants(tstruct, t.trap_ix))
    end
    return s
end

# Map each node's global trap index to the index of the old context that held it.
function _trap_owner_map(net_contexts)
    owner = Dict{Int,Int}()
    for (ci, ctx) in enumerate(net_contexts), g in ctx.global_ix
        owner[g] = ci
    end
    return owner
end

# Per-context gate: decide which freshly-rebuilt `components` can carry over an old
# context verbatim — skipping both its commit and its ODE re-prediction.  A component
# is reusable iff it is byte-identical to exactly one old context (same node set), that
# context is NOT affected, and its external inflow is unchanged after the reconcile
# (recomputed cheaply — no ODE).  A merge/grow/split yields a node set that matches no
# single old context, so it always rebuilds — the correctness guard against overlap.
# Returns a vector aligned with `components`: the reusable old-context index, or 0.
function _reuse_plan(components, net_contexts, affected, trap_owner, tstruct, rateinfo)
    plan = zeros(Int, length(components))
    for (k, net) in enumerate(components)
        S      = Set(t.trap_ix for t in net.traps)
        owners = Set(get(trap_owner, g, 0) for g in S)
        (length(owners) == 1 && 0 ∉ owners) || continue
        j = first(owners)
        (j ∉ affected && Set(net_contexts[j].global_ix) == S) || continue
        _external_inflow(net, rateinfo, tstruct) == net_contexts[j].extern_inflow &&
            (plan[k] = j)
    end
    return plan
end

# The old contexts that must be committed to `cur_time`: those feeding any rebuilt
# component (their nodes supply its initial volumes) plus every affected context (its
# boundary/just-exited volumes are read from the committed dict).  Reused contexts are
# deliberately excluded — leaving their state at the earlier `last_solve_time` keeps
# `_network_amount_updates` from emitting spurious updates for them (plan §9).
function _contexts_to_commit(components, reuse, trap_owner, affected)
    need = Set{Int}(affected)
    for (k, net) in enumerate(components)
        reuse[k] == 0 || continue
        for t in net.traps
            j = get(trap_owner, t.trap_ix, 0)
            j == 0 || push!(need, j)
        end
    end
    return need
end

# Commit the selected old contexts to `cur_time`, pin their full traps to exact
# capacity, and collect committed node volumes by global trap index.
function _commit_contexts!(net_contexts, which, tstruct, infiltration, z_vol_tables,
                           cur_time, full_set)
    committed = Dict{Int,Float64}()
    for ci in which
        ctx = net_contexts[ci]
        _commit_network!(ctx, tstruct, infiltration, z_vol_tables, cur_time)
        _clamp_full_traps!(ctx, tstruct, full_set)
        for (i, g) in enumerate(ctx.global_ix)
            committed[g] = ctx.state[i]
        end
    end
    return committed
end

# Overlay the boundary volumes of fired traps that did NOT just become full.  An :empty
# parent drops to 0 and EXPOSES its immediate children — they go from full to transitory
# (just below their own capacity), so they leave `full_traps` (the caller already flipped
# them — §`_expand_empty_fill_updates`) and start at prevfloat(C_child).  An :unspill trap
# drops to prevfloat(C).
function _apply_fired_boundaries!(committed, fired_kind, tstruct)
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
    return committed
end

# Assemble the post-touch context vector.  A reusable component carries its old context
# object over verbatim; a rebuilt one is constructed from `state0` — the committed node
# value, exact capacity for a full trap, or a constant-rate projection for a trap newly
# absorbed from the plain path (its `cur_amounts` entry is stale, holding the volume at
# its last committed event; project it forward under the pre-reconcile saved inflow) —
# and gets a fresh next-event prediction.
function _assemble_contexts(components, reuse, net_contexts, committed, full_set, seeds,
                            tstruct, rateinfo, infiltration, z_vol_tables,
                            cur_amounts, cur_time, endtime)
    project(g) = first(fill_trap_until(g, rateinfo, cur_amounts[g], cur_time,
                                       tstruct, z_vol_tables, use_saved=true))
    state0(g)  = g in full_set        ? _own_capacity(tstruct, g) :
                 haskey(committed, g) ? committed[g]              :
                                        project(g)
    out = DynNetworkContext[]
    for (k, net) in enumerate(components)
        if reuse[k] != 0
            push!(out, net_contexts[reuse[k]])
        else
            c = _make_context(net, tstruct, rateinfo, seeds, state0, cur_time)
            _predict_network!(c, tstruct, infiltration, z_vol_tables, endtime)
            push!(out, c)
        end
    end
    return out
end

# ----------------------------------------------------------------------------
# Advance the touched networks to `cur_time` and rebuild after a status change (plan §3):
# the STRUCTURE is retraced from the seed pool (via the incremental cache) so components
# merge/split correctly, but the expensive ODE solves (commit + predict) are GATED per
# context — only changed components are re-solved, unchanged ones carried over verbatim
# (`_reuse_plan`).  Two-level gating: the call site skips quiet events (plan D4/§8); this
# per-context gate then keeps the solve count near the number of changed components (usually
# one).  Must run AFTER `_update_flow!`.  Returns the new contexts, refreshed net_trap_set /
# net_covered_set, and the committed-volume dict (for `_network_amount_updates`).
function _touch_networks!(net_contexts, changetimeest, sgraph, tstruct, dyn_traps, culverts,
                          filled_traps, cur_amounts, rateinfo, z_vol_tables, infiltration,
                          fill_updates, old_covered, cur_time, endtime, subnet_cache)
    isempty(net_contexts) &&
        return net_contexts, Set{Int}(), Set{Int}(), Dict{Int,Float64}()
    full_traps = findall(filled_traps)
    full_set   = Set(full_traps)

    # Which old contexts changed this event, and what each fired trap did.
    inflow_updated = Set(u.index for u in getinflowupdates(rateinfo))
    affected   = _affected_contexts(net_contexts, fill_updates, inflow_updated)
    fired_kind = _capture_fired_kinds(net_contexts, fill_updates)

    # Retrace the structure from the global seeds, then reconcile the spillgraph to the
    # new coverage BEFORE any external inflow is read (absorbed traps' double-counted
    # deposits withdrawn; traps leaving coverage regain their edges).
    nbs_objs    = _nbs_elements(tstruct)
    seeds       = _dyn_seeds(tstruct, dyn_traps, culverts, nbs_objs)
    # Incremental retrace: the base per-seed subnet traces are cached (only chains a fired
    # trap touched are retraced); the culvert endpoint-expansion and culvert-aware merge run
    # on top, exactly as in `setup_network`.  NBS elements (like culverts) force the full merge.
    components  = setup_network_cached(tstruct, seeds, full_traps, subnet_cache;
                                       culverts=culverts, nbs=nbs_objs)
    new_covered = _covered_of(components, tstruct)
    old_covered != new_covered &&
        _reconcile_spillgraph!(sgraph, rateinfo, filled_traps, new_covered, tstruct)

    # Gate the solves: carry over unchanged components; commit only the contexts feeding
    # a rebuilt component (plus the affected ones), then overlay fired-trap boundaries.
    trap_owner = _trap_owner_map(net_contexts)
    reuse      = _reuse_plan(components, net_contexts, affected, trap_owner, tstruct, rateinfo)
    to_commit  = _contexts_to_commit(components, reuse, trap_owner, affected)
    committed  = _commit_contexts!(net_contexts, to_commit, tstruct, infiltration,
                                   z_vol_tables, cur_time, full_set)
    _apply_fired_boundaries!(committed, fired_kind, tstruct)

    # Build the new context vector (reuse or rebuild+predict) and publish the estimates.
    new_contexts = _assemble_contexts(components, reuse, net_contexts, committed, full_set,
                                      seeds, tstruct, rateinfo, infiltration, z_vol_tables,
                                      cur_amounts, cur_time, endtime)
    nts = _net_trap_set(new_contexts)
    ncs = _net_covered_set(new_contexts, tstruct)
    _apply_network_changetimeest!(changetimeest, new_contexts, ncs)
    return new_contexts, nts, ncs, committed
end

# ----------------------------------------------------------------------------
# Weather-period boundary finalization (plan §10): advance every context to `endtime`
# under its cached external inflow and read the settled volumes.  Network traps follow the
# multi-trap ODE, so they can't use the constant-rate `fill_trap_until` projection — and
# these exact volumes are what the NEXT period rebuilds from.  The advance is event-free
# (any event ≤ endtime was already processed), so it settles cleanly.  A node takes its
# settled ODE volume; a subsumed full descendant sits at capacity.
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
