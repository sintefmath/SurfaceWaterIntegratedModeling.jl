# ============================================================================
# DynNetworkContext: drives one DynNetwork between topology changes inside
# `fill_sequence` — the integration layer between the analytical fill-sequence
# machinery and the multi-trap ODE solver `solveDynNetwork!`.
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
- `inflow_regions`: the base spill regions whose `trap_inflow` feeds this network's
  external inflow — the leaf descendants of the nodes (a leaf trap and its spill region
  share the same id), summed by `_external_inflow`; used by the touch test in the event
  loop.
- `last_solve_time`: absolute time the committed `state` refers to.
- `extern_inflow`: per-node external inflow the `state` is currently evolving
  under (cached at the last touch; the rate the commit/predict solves use).
- `next_event`: cached prediction `(; time, trap, kind)` of the network's next
  event, with `time` an ABSOLUTE timestamp and `trap` a global trap index.
"""
mutable struct DynNetworkContext
    net             ::DynNetwork
    state           ::Vector{Float64}            # trap volumes, then appended NBS layer states
    global_ix       ::Vector{Int}
    inflow_regions  ::Set{Int}                   # base spill regions feeding external inflow
    last_solve_time ::Float64
    extern_inflow   ::Vector{Float64}
    runoff          ::Matrix{Float64}            # the oblivious runoff grid (rateinfo.runoff), read
                                                 # by the solver to build the NBS correction plan
    next_event      ::NamedTuple
end

# ----------------------------------------------------------------------------
# Each node's leaf spill regions (`lowest_subtraps_for[t]`, == [t] for a leaf), in
# `net.traps` order.  Both the per-node `extern_inflow` and the flat `inflow_regions`
# touch key are projections of this one relation, so the two stay in lock-step by
# construction.
_node_leaf_regions(net::DynNetwork, tstruct) =
    [tstruct.lowest_subtraps_for[t.trap_ix] for t in net.traps]

# Per-node external inflow: each node's sum of `trap_inflow` over its leaf regions.  Summed
# from the leaves, not `getinflow(node)` directly, because `_reconcile_spillgraph!` withdraws
# covered children's flow from the parent's `trap_inflow`, leaving it stale — leaves stay current.
_external_inflow(rateinfo, node_leaves) =
    Float64[sum(getinflow(rateinfo, leaf) for leaf in leaves) for leaves in node_leaves]

# The base spill regions feeding this network's external inflow: the nodes' leaf regions
# collected flat (a leaf trap and its spill region share the same id).  The nodes partition
# their leaves — none is shared — so this is a plain flatten, not a merge.  Used for the touch
# test (a network is touched when one of these appears in `getinflowupdates`).
_inflow_regions(node_leaves) = Set{Int}(Iterators.flatten(node_leaves))

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

# `net_trap_set` plus every descendant subsumed under a network parent node: those
# subsumed children are full and static while their parent is a node, so they must be
# excluded from the standard changetime machinery.
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
                           tmax = endtime - ctx.last_solve_time, zvt = z_vol_tables,
                           runoff = ctx.runoff)
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
                           tmax = T_commit - ctx.last_solve_time, zvt = z_vol_tables,
                           runoff = ctx.runoff)
    ctx.last_solve_time = T_commit
    return res
end

# ----------------------------------------------------------------------------
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

# The NBS layer-state block for `net`, read from the persistent `nbs_state` store in
# `net.nbs` order (matching `_make_context`'s state layout and `_build_nbs_plan`'s
# state_base offsets).  Each placement in `net.nbs` carries its own layer model and stable
# `id`; a placement not yet in the store starts empty (zeros).
function _nbs_layer_block(net::DynNetwork, nbs_state)
    block = Float64[]
    for nb in net.nbs
        L      = length(nb.system.layers)
        stored = get(nbs_state, nb.id, nothing)
        append!(block, stored === nothing ? zeros(Float64, L) : stored)
    end
    return block
end

# Write a context's current layer states back into the persistent `nbs_state` store
# (the single source of truth across events and weather periods).  Keyed by placement `id`,
# so a later rebuild that regroups NBS into new components restores each layer correctly.
function _store_nbs_state!(nbs_state, ctx::DynNetworkContext)
    nt   = length(ctx.global_ix)
    base = 0
    for nb in ctx.net.nbs
        L = length(nb.system.layers)
        nbs_state[nb.id] = ctx.state[(nt + base + 1):(nt + base + L)]
        base += L
    end
    return nbs_state
end

# Build one context from a single network component.  `state0(g)` supplies the initial
# committed volume for the global trap index `g`.  NBS placements travel on `net.nbs`.
function _make_context(net::DynNetwork, tstruct, rateinfo, state0, cur_time,
                       nbs_state = Dict{Int,Vector{Float64}}())
    global_ix = Int[t.trap_ix for t in net.traps]
    # Trap volumes, then one appended state per NBS layer (in `net.nbs` order, matching
    # `_build_nbs_plan`'s state_base offsets), restored from the persistent store so NBS
    # storage carries across events and weather periods.
    state     = vcat(Float64[state0(g) for g in global_ix],
                     _nbs_layer_block(net, nbs_state))
    # Both inflow fields are projections of the same per-node leaf-region relation (union of
    # the regions vs. their per-node summed inflow), so compute that relation once.
    node_leaves = _node_leaf_regions(net, tstruct)
    # Hold the oblivious runoff grid so the solver can build the NBS correction plan from it.
    return DynNetworkContext(net, state, global_ix,
                             _inflow_regions(node_leaves),
                             cur_time,
                             _external_inflow(rateinfo, node_leaves),
                             rateinfo.runoff,
                             (; time = Inf, trap = 0, kind = :none))
end

# Own storage capacity of a trap (net of subtraps) — the same quantity the solver
# uses as `geom.capacity` and `fill_sequence` records for a full trap.
_own_capacity(tstruct, t::Int) = tstruct.trapvolumes[t] - tstruct.subvolumes[t]

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

# Overlay the boundary volumes of fired traps that did NOT just become full.  An :empty
# parent drops to 0 and EXPOSES its immediate children — they go from full to transitory
# (just below their own capacity), so they leave `full_traps` (the caller already flipped
# them via `_expand_empty_fill_updates`) and start at prevfloat(C_child).  An :unspill trap
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

# ----------------------------------------------------------------------------
# Weather-period boundary finalization: advance every context to `endtime`
# under its cached external inflow and read the settled volumes.  Network traps follow the
# multi-trap ODE, so they can't use the constant-rate `fill_trap_until` projection — and
# these exact volumes are what the NEXT period rebuilds from.  The advance is event-free
# (any event ≤ endtime was already processed), so it settles cleanly.  A node takes its
# settled ODE volume; a subsumed full descendant sits at capacity.
function _finalize_networks!(cur_amounts, net_contexts, tstruct, infiltration,
                             z_vol_tables, cur_time, endtime,
                             nbs_state = Dict{Int,Vector{Float64}}())
    isempty(net_contexts) && return cur_amounts
    stamp = min(cur_time, endtime)
    for ctx in net_contexts
        ctx.last_solve_time < endtime &&
            _commit_network!(ctx, tstruct, infiltration, z_vol_tables, endtime)
        _store_nbs_state!(nbs_state, ctx)          # carry NBS layer storage into the next period
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

# Amount updates for network traps whose volume the network owns: a node carries
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
