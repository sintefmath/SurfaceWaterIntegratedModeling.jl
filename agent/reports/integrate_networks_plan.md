# Integration plan: `solveDynNetwork!` into `fill_sequence`

## Background

`fill_sequence` maintains constant rates between events and estimates change
times analytically (or via `fill_trap_until`'s single-trap ODE).  Once
rate-limited culverts are present, flows between traps in the same
hydraulically-connected sub-system are no longer constant between events, so
that analytical approach breaks.  `solveDynNetwork!` solves the multi-trap ODE
for such sub-systems and returns the first topology-changing event with an exact
time.

---

## 1. New `fill_sequence` arguments

```julia
function fill_sequence(tstruct, weather_events;
                       time_slack   = 0.0,
                       infiltration = nothing,
                       dyn_coords   ::Union{Nothing, Vector{CartesianIndex{2}}} = nothing,
                       culverts     ::Vector{DynCulvert} = DynCulvert[],
                       verbose      = false)
```

`dyn_coords` are the grid-cell entry points that seed dynamic networks (cells
upstream of or inside culvert-connected trap groups).  `culverts` is the list of
`DynCulvert` objects.  When both are empty/`nothing`, behaviour is unchanged.

---

## 2. New mutable state: `DynNetworkContext`

A small struct holds everything needed to run and maintain one `DynNetwork`
between topology changes:

```julia
mutable struct DynNetworkContext
    net             ::DynNetwork       # current network (rebuilt when topology changes)
    state           ::Vector{Float64}  # COMMITTED trap volumes at last_solve_time
                                       # (net of subtraps), indexed as net.traps
    global_ix       ::Vector{Int}      # tstruct trap index for each net.traps entry
    last_solve_time ::Float64          # absolute time the committed state refers to
    extern_inflow   ::Vector{Float64}  # external inflow per net.traps entry under which
                                       # `state` is currently evolving (cached at last touch)
    next_event      ::NamedTuple       # cached prediction of the next network event,
                                       # (; time, trap, kind), with time ABSOLUTE
end
```

Invariant: `state` always holds the network's volumes at the **committed**
absolute time `last_solve_time` — never a provisional/over-advanced state (see
§5).  `next_event.time` is an absolute timestamp (the relative elapsed time
returned by `solveDynNetwork!` plus the solve's start time).  `extern_inflow` is
the per-trap external inflow the network is currently evolving under; it is
constant between *touches* (§3), which is what lets an untouched network be
advanced later in a single bounded solve.

Two `Set{Int}`s of global trap indices are maintained in `fill_sequence`'s outer
scope for O(1) membership tests:

- `net_trap_set` — the trap indices that are *nodes* of some active `DynNetwork`.
- `net_covered_set` — `net_trap_set` plus all descendants subsumed under a network
  parent node (Design A, §3 Build).  These subsumed children are full and static
  while their parent is an active node, so they must be excluded from the standard
  changetime machinery (§7); their next change is governed by the parent's ODE.

---

## 3. Lifecycle of `DynNetworkContext` objects inside `fill_sequence`

A context is **built** once per `WeatherEvent`, and thereafter only **touched**
(advanced + refreshed) at events that actually affect it.  This gating is what
avoids re-solving every network on every event.

### Touch criterion

A network `N` is touched at a committed event time `T_commit` iff one of:

- a trap of `N` is in `fill_updates` — `N` fired its own event, or a member trap
  filled (topology change → rebuild); or
- `N`'s external inflow changed this event, i.e.
  `rateinfo.trap_inflow[global_ix] ≠ ctx.extern_inflow` on any entry.

Networks meeting neither condition are left **untouched**: their `state`,
`last_solve_time`, `extern_inflow`, cached `next_event`, and `changetimeest`
entries all remain valid, because their external inflow has been constant since
their last touch — so their predicted (absolute) event time is unchanged and
their state can be advanced later in a single bounded solve.

### Build

`setup_network(tstruct, dyn_coords, filled_traps; culverts=culverts)` — may return
several disjoint `DynNetwork` components, each becoming a `DynNetworkContext`.
`state` is initialised from `cur_amounts`, `extern_inflow` from the current
`rateinfo` (gross composite inflow, §4), and the first event is predicted (§5, §6).

The constructor enforces the **subsume-parents** invariant (Design A, bug B1): a
parent trap is always a single node subsuming all its descendants.  The
full-parent case is already handled by `flow_path_from` (it records the uppermost
full supertrap and deletes its footprint cells from the path); the constructor
additionally subsumes a *terminal-unfilled* parent's full descendants, dropping
the recorded child traps and the path segments internal to the parent's
footprint.  Consequently a network node is either a leaf trap or a parent that
subsumes its whole subtree — never a parent with some children also present as
nodes.

### Touch (`_touch_network!`)

1. **Commit** — advance `state` from `last_solve_time` to `T_commit` with a
   `tmax`-bounded solve (§5) under the *cached* `extern_inflow` (the rates valid
   over that interval).  Set `last_solve_time = T_commit`.  Because the commit
   uses the cached inflow, it does not depend on `rateinfo` still holding
   pre-event values, so it may run *after* `_update_flow!`.
2. **Boundary clamp** — if a trap of this context fired, set it to its post-event
   boundary value per `solveDynNetwork!`'s caller protocol (`:fill` → C,
   `:empty` → 0, `:unspill` → `prevfloat(C)`; kind taken from `ctx.next_event`).
3. **Rebuild** (only if a member trap filled) — `setup_network` with the updated
   `filled_traps`, then **remap** `state` onto the new `net.traps` by global trap
   index (see below).  Update `global_ix` and `net_trap_set`.
4. **Refresh** `extern_inflow = _external_inflow(net, rateinfo, tstruct)` (gross
   composite inflow per node, §4).
5. **Predict** the next event on a *copy* of `state` (§5/§6); store `next_event`
   (absolute) and the `changetimeest` entries.

### State remap across a rebuild

When `setup_network` returns a new `net.traps` ordering / membership, build the
new `state` by global trap index:

- a trap present in the *old* context carries its committed volume forward
  (looked up via the old `global_ix`);
- a **newly-introduced** trap is seeded from `cur_amounts[trap_ix]` (its
  fill_sequence-tracked amount), or from its boundary value if it just filled.

This is the only place `state` membership changes; non-rebuild touches keep
`net.traps` fixed.

---

## 4. The external-inflow problem

`solveDynNetwork!` expects `inflow[i]` = the inflow to trap `i` from **sources
outside** the `DynNetwork` (terrain and non-network traps).  Intra-network
routing is handled inside the ODE by `_route_flow`.  Using
`rateinfo.trap_inflow[t]` directly is wrong when another full network trap is
currently spilling to `t` — the ODE would double-count that contribution.

**Proposed approach**: treat network traps as always-empty in the spillgraph.
Network trap fill/unfill events are filtered out before being passed to
`update_spillgraph!`, so network traps never acquire out-edges in `sgraph`.
Since `compute_flow` only routes outflow from traps that have out-edges in the
spillgraph, it naturally skips network traps entirely — no modification to
`compute_flow` is needed.  Effect:

- `rateinfo.trap_inflow[t]` for network traps receives only terrain +
  non-network-trap contributions.
- No double-counting: the ODE handles intra-network flows; the spillgraph and
  `compute_flow` handle everything else.

**Why there is no outflow problem to non-network traps**: `setup_network` follows
the terrain flow graph all the way downstream from each entry point, including
every trap it encounters, and always terminates at the first unfilled trap.  By
construction, any trap that a full network trap would spill into is already part
of the network.  The only true exit from a `DynNetwork` is the domain boundary,
which `_route_flow` already handles correctly by dropping that flow.  There is
therefore no need to inject network outflow into `rateinfo` for downstream
non-network traps.

**Parent nodes — gross composite inflow (Design A, bug B1).**  A parent trap is
always represented as a *single* node that subsumes all its descendants (children
are never separate nodes; the constructor enforces this — §3 Build).  Such a node
is fed the **gross** composite inflow

```
external_inflow[node] = Σ_{leaf ∈ leaf-descendants(node)} trap_inflow[leaf]
                      = trap_inflow[node]   for a leaf node
```

This matches the solver's existing parent geometry, which charges the whole
composite footprint infiltration, so no solver change is needed.  Using
`trap_inflow[P]` directly would *double-count* infiltration — it is already net
of child infiltration (`signed_flow = inflow − Smax`), and the solver deducts the
composite footprint infiltration again (review bug **B1**).

---

## 5. Bounded integration and the state-commit protocol

**Problem.**  `solveDynNetwork!` advances `state` to *its own* next event T_X, but
the globally-committed next event in the branch-and-bound may be an earlier
non-network event at T_Y < T_X.  We must be able to obtain network volumes at the
committed time, not only at T_X.

**Design decision (API extension).**  Add an optional `tmax` to `solveDynNetwork!`:

```julia
solveDynNetwork!(state, tstruct, net, infiltration, inflow; tmax = Inf, ...)
```

Integration runs over `(0.0, tmax)` instead of `(0.0, Inf)`.  If a topology event
occurs at elapsed `t_e <= tmax`, it is returned as before (with `state` at
`t_e`).  Otherwise integration stops at `tmax` and returns
`(time = tmax, trap = 0, kind = :none)` with `state` at `tmax`.  As today, `time`
is **elapsed** time from the solve start; callers add the absolute start time.

**Two-phase use per context** (keeps `state` equal to committed volumes at
`last_solve_time`, never over-advanced):

- **Predict** — on a *copy* of the committed state, solve with
  `tmax = endtime - last_solve_time` to find the network's next event.  Store
  `next_event.time = elapsed + last_solve_time` (absolute) and the corresponding
  `changetimeest` entry.  The committed `state` is untouched.
- **Commit** — once the inner loop knows the committed time `T_commit` (the value
  `_identify_next_status_change!` returns), advance the committed `state` in place
  with `tmax = T_commit - last_solve_time`, then set `last_solve_time = T_commit`.
  `state` now holds exact volumes at `T_commit`; if `T_commit == T_X` this same
  call re-confirms the topology event.

The commit advance must use the rates valid over `[last_solve_time, T_commit]`.
Rather than ordering it before `_update_flow!`, the context caches those rates in
`ctx.extern_inflow` (refreshed only at a touch, §3), so the commit reads the
cached value and is order-independent of the spillgraph/flow update.

---

## 6. Changes to `_set_initial_changetime_estimates`

Current behaviour: loops over all traps, returns standard analytical estimates.

Proposed: after computing standard estimates for all traps, predict each
network's first event (on a copy of its state, per §5) and override its trap
estimates:

```
for each DynNetworkContext ctx:
    external_inflow = _external_inflow(ctx.net, rateinfo, tstruct)
        # external_inflow[i] = Σ trap_inflow[leaf descendants of global_ix[i]]
        # (= trap_inflow[global_ix[i]] for a leaf node; gross composite inflow, §4)
    result = solveDynNetwork!(copy(ctx.state), tstruct, ctx.net, infiltration,
                              external_inflow; tmax = endtime - cur_time, zvt = z_vol_tables)
    ctx.next_event = (; time = result.time + cur_time, result.trap, result.kind)
    if result.kind != :none:
        changetimeest[result.trap] = ChangeTimeEstimate(result.kind == :fill,
                                                        ctx.next_event.time,
                                                        ctx.next_event.time)
    for all other traps t in ctx.global_ix:
        changetimeest[t] = ChangeTimeEstimate(false, Inf, Inf)
```

The exact-estimate property (`min == max`) means network traps are resolved
immediately in `_identify_next_status_change!`'s branch-and-bound without
additional ODE calls — they never reach `_compute_exact_changetime`.

---

## 7. Changes to `_update_changetime_estimates!`

This function is called at the top of `_identify_next_status_change!` to refresh
estimates for traps whose inflow changed (via `getinflowupdates(rateinfo)`).

**Change**: skip any trap in `net_covered_set` (network nodes *and* their subsumed
descendants).  Network nodes' estimates are set by the most recent network
prediction; subsumed descendants are full and static under their parent node.
Calling `_compute_changetime_estimate` on either would produce incorrect results
(it assumes constant rates, and would let a subsumed full child spuriously become
a candidate to empty).

`net_covered_set` is threaded through as an additional argument.  The candidate
selection in `_identify_next_status_change!` must likewise ignore
`net_covered_set` traps.

---

## 8. Changes in the main event loop (`_fill_sequence_for_weather_event!`)

The existing sequence becomes (new steps marked ►):

```julia
cur_time, fill_updates = _identify_next_status_change!(...)   # cur_time == T_commit
# break as before if cur_time > endtime or fill_updates empty

for u in fill_updates; filled_traps[u.index] = u.value; end

# Only non-network traps enter the spillgraph machinery
non_net_fill_updates = filter(u -> u.index ∉ net_trap_set, fill_updates)
graph_updates = update_spillgraph!(sgraph, non_net_fill_updates, tstruct)
setsavepoint!(rateinfo)
_update_flow!(rateinfo, graph_updates, tstruct, sgraph)

► # Touch only the AFFECTED networks (commit → clamp → rebuild → refresh → predict, §3)
►  for ctx in net_contexts
►      fired          = any(u.index ∈ ctx.global_ix for u in fill_updates)
►      new_ext        = _external_inflow(ctx.net, rateinfo, tstruct)   # gross, §4
►      inflow_changed = new_ext != ctx.extern_inflow
►      (fired || inflow_changed) &&
►          _touch_network!(ctx, tstruct, dyn_coords, culverts, filled_traps,
►                          cur_amounts, rateinfo, z_vol_tables, cur_time, endtime)
►  end

# amount updates + SpillEvent as before (see §9)
```

`_touch_network!` performs the full §3 touch sequence on one context: commit to
`cur_time` under the cached `extern_inflow`, clamp a fired trap, rebuild + remap
on a member fill, refresh `extern_inflow`, and re-predict (updating
`changetimeest` as in §6).  Networks that neither fired nor saw an inflow change
are skipped entirely — no ODE solve.

---

## 9. `cur_amounts` and SpillEvent recording

**Network-trap volumes come only from the network solver.**  A network trap's
water amount is *never* computed with the constant-rate `fill_trap_until` — it is
read exclusively from the integrated `ctx.state`.  When a context is touched (§3),
its commit step leaves `ctx.state` holding exact network-trap volumes at
`T_commit`.  Therefore:

- `cur_amounts[t]` for a network trap `t` is read directly from `ctx.state`
  (mapped via `global_ix`), at the events where its context is touched.
- `amount_updates` must include entries for every trap of a touched network whose
  volume changed materially, not just the trigger trap.  (Untouched networks have
  unchanged dynamics; their amounts are interpolated between touch events, the
  same approximation the existing code already makes for ODE-governed traps.)

The two existing call sites that apply `fill_trap_until` to *all* affected/stale
traps must therefore exclude network traps (`net_trap_set`) and take their
amounts from `ctx.state` instead:

- `_update_affected_amounts` — skip any trap in `net_trap_set`; emit its amount
  from the corresponding `ctx.state` entry.
- the end-of-weather-period finalization loop in
  `_fill_sequence_for_weather_event!` — same exclusion; network-trap amounts at
  `endtime` come from advancing `ctx.state` to `endtime` via the `tmax`-bounded
  commit mechanism (§5; see open question 3), not from `fill_trap_until`.

---

## 10. `_compute_exact_changetime` — no change needed

Network traps always have `min == max` in their `ChangeTimeEstimate`, so
`_compute_exact_changetime` returns early on the existing
`changetimes[trap].min == changetimes[trap].max` fast path.

---

## Summary of touched locations

| Location | Change |
|---|---|
| `fill_sequence` signature | Add `dyn_coords`, `culverts` |
| `fill_sequence` body | Build `DynNetworkContext`s; maintain `net_trap_set` and `net_covered_set` |
| `setup_network` / `_subnetwork` | Subsume a terminal-unfilled parent's full descendants so a parent is always one node (Design A, B1) |
| `solveDynNetwork!` | Add optional `tmax`; return `(tmax, 0, :none)` with `state` at `tmax` when no event fires before `tmax` |
| `_fill_sequence_for_weather_event!` | Filter `fill_updates` from the spillgraph; after `_update_flow!`, touch only the affected networks (gated on fired-or-inflow-changed) |
| `_set_initial_changetime_estimates` | After standard estimates, predict each network's first event (on a copy) and override its trap estimates |
| `_update_changetime_estimates!` | Skip traps in `net_trap_set` |
| `_update_affected_amounts` | Skip `net_trap_set`; emit network-trap amounts from `ctx.state` |
| end-of-period finalization loop | Skip `net_trap_set`; advance `ctx.state` to `endtime` instead of `fill_trap_until` |
| NEW `_external_inflow` | Per-node gross composite inflow `Σ trap_inflow[leaf descendants]` (= `trap_inflow` for a leaf), §4 |
| NEW `_build_dyn_networks` | Construct `DynNetworkContext`s from `dyn_coords`, `culverts`, and `filled_traps` |
| NEW `_touch_network!` | Commit one context to `cur_time` (cached `extern_inflow`), clamp a fired trap, rebuild + remap on a member fill, refresh `extern_inflow`, re-predict |

---

## Open questions / deferred decisions

1. **`dyn_coords` vs culvert-derived coords**: The caller currently supplies
   `dyn_coords` explicitly.  An alternative is to derive them automatically from
   culvert inlet/outlet cells, removing the burden from the caller.
   `setup_network` handles merging, so both approaches are equivalent in effect.

2. **`path_inflow` population**: Flow-path cells that receive terrain runoff
   (from cells outside any trap footprint, upstream of the first network trap)
   should contribute to `path_inflow` in `solveDynNetwork!`.  The current plan
   leaves this at zero.  Populating it correctly requires summing
   `rateinfo.runoff` at each path's first cell; this can be added as a
   refinement without changing the overall architecture.

3. **Weather-event boundary**: At the end of a weather period, each network's
   committed `state` is advanced to `endtime` using the same `tmax`-bounded
   commit mechanism (§5), so no `fill_trap_until` projection is needed for network
   traps.

4. **`time_slack`**: Still unimplemented for the standard path; same deferral
   applies here.

5. **Uphill / reverse culverts**: `solveDynNetwork!` (via `_topological_order`)
   will throw on a cyclic network caused by an uphill culvert.  No new handling
   is needed here; the existing error is informative.
