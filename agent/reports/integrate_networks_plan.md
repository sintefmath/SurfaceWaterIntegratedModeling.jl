# Integration plan: `solveDynNetwork!` into `fill_sequence`

## Background

`fill_sequence` maintains constant rates between events and estimates change
times analytically (or via `fill_trap_until`'s single-trap ODE). Once
rate-limited culverts are present, flows between traps in the same
hydraulically-connected sub-system are no longer constant between events, so that
analytical approach breaks. `solveDynNetwork!` solves the multi-trap ODE for such
sub-systems and returns the first topology-changing event with an *exact* time.

The integration isolates the affected traps into one or more `DynNetwork`s, solves
each as an ODE only when something actually changes it, and feeds the resulting
exact change-times back into `fill_sequence`'s existing branch-and-bound. All
non-network traps keep the current analytical fast path.

### Already in place (landed and tested)

These constructor-side pieces exist on the branch and the full dynamics suite is
green; the rest of this document is not yet implemented.

- **Design A — subsume parents** (`_subnetwork` / `_subsume_terminal_parent`): a
  parent trap is always one network node subsuming its whole subtree.
- **Terminal-detection fix** in `_subnetwork`: handles the trap-terminated chain
  produced by the current `flow_path_from` (wrap-around among full siblings of an
  unfilled parent).
- **`spill_path == -1` sentinel**: a full trap spilling straight out of the domain
  carries `spill_path == -1`; `_validate_network` accepts FULL with `!= 0`
  (`> 0` or `-1`) and still rejects FULL with `0`. Routing/ordering already treat
  any `spill_path <= 0` as "no in-network successor".

---

## 1. New `fill_sequence` arguments

```julia
function fill_sequence(tstruct, weather_events;
                       time_slack   = 0.0,
                       infiltration = nothing,
                       dyn_traps    ::Union{Nothing, Vector{Int}} = nothing,
                       culverts     ::Vector{DynCulvert} = DynCulvert[],
                       verbose      = false)
```

`dyn_traps` are explicit caller-supplied **trap indices** to include as dynamic
networks; `culverts` is the list of `DynCulvert` objects. Networks are seeded by
the **union** of `dyn_traps` and the traps/paths touched by every culvert's
inlet/outlet cells (Decision D1), so a culvert-connected trap group is
auto-detected even with no `dyn_traps` naming it. Passing trap indices (rather than
grid cells restricted to footprints) makes "no hanging paths" true by construction
— `_build_dyn_networks` seeds `setup_network` from a representative footprint cell
of each named trap, so the network always starts *at* a trap (Decision D2). The
only validation is that each index is a valid trap. When both inputs are
empty/`nothing`, behaviour is identical to today.

---

## 2. New mutable state: `DynNetworkContext`

One context runs and maintains a single `DynNetwork` between topology changes:

```julia
mutable struct DynNetworkContext
    net             ::DynNetwork       # current network (rebuilt on a topology change)
    state           ::Vector{Float64}  # COMMITTED trap volumes at last_solve_time,
                                       # net of subtraps, indexed as net.traps
    global_ix       ::Vector{Int}      # tstruct trap index for each net.traps entry
    inflow_sources  ::Set{Int}         # global trap indices whose trap_inflow feeds
                                       # this network's external inflow (the leaf
                                       # descendants summed by _external_inflow);
                                       # used for the touch test (§3, §8)
    last_solve_time ::Float64          # absolute time the committed state refers to
    extern_inflow   ::Vector{Float64}  # per-node external inflow the state is evolving
                                       # under (cached at last touch, for the commit)
    next_event      ::NamedTuple       # (; time, trap, kind); time ABSOLUTE
end
```

Invariants and notes:

- `state` always holds the network's volumes at the **committed** absolute time
  `last_solve_time` — never over-advanced past it (§5).
- `next_event.time` is absolute (the elapsed time `solveDynNetwork!` returns plus
  the solve's start time).
- `extern_inflow` is the per-node inflow the state is currently evolving under. It
  is the *cached rate* the commit advance uses (§5); it is **not** used to detect
  inflow changes — that is the job of `inflow_sources` + `getinflowupdates`
  (Decision D4).

Two `Set{Int}`s of global trap indices live in `fill_sequence`'s outer scope:

- `net_trap_set` — trap indices that are *nodes* of some active `DynNetwork`.
- `net_covered_set` — `net_trap_set` plus all descendants subsumed under a network
  parent node (Design A). Subsumed children are full and static while their parent
  is an active node, so they are excluded from the standard changetime machinery
  (§7); their next change is governed by the parent's ODE.

---

## 3. Lifecycle of a `DynNetworkContext`

A context is **built** once per `WeatherEvent`, then only **touched** (advanced +
refreshed) at events that actually affect it. This gating is what avoids
re-solving every network on every event.

### Touch criterion (Decision D4)

A network `N` is touched at committed time `T_commit` iff one of:

- a trap of `N` is in `fill_updates` — `N` fired its own event, or a member trap
  filled (topology change → rebuild); or
- some trap in `ctx.inflow_sources` appears in `getinflowupdates(rateinfo)` — the
  network's external inflow changed this event.

`getinflowupdates(rateinfo)` is the savepoint-tracked set of traps whose
`trap_inflow` changed since the last `setsavepoint!`; the rest of the inner loop
already drives its staleness off the same signal. We therefore reuse it directly
rather than recomputing `_external_inflow` for every network and float-comparing —
cheaper and exactly consistent with how `_update_changetime_estimates!` decides
what is stale.

Networks meeting neither condition are left **untouched**: their `state`,
`last_solve_time`, `extern_inflow`, cached `next_event`, and `changetimeest`
entries all remain valid, because their external inflow has been constant since
their last touch — so the predicted (absolute) event time is unchanged and the
state can be advanced later in a single bounded solve.

### Build (`_build_dyn_networks`)

1. Validate each `dyn_trap` is a valid trap index (Decision D2), throwing otherwise.
2. Convert each `dyn_trap` to a representative footprint cell
   (e.g. `CI[tstruct.footprints[t][1]]`) and form
   `seeds = unique(footprint_cells(dyn_traps) ∪ {cv.inlet, cv.outlet for cv in culverts})`.
3. `setup_network(tstruct, seeds, filled_traps; culverts=culverts)` — may return
   several disjoint components, each becoming a `DynNetworkContext`.
4. Initialise `state` from `cur_amounts`, `inflow_sources` and `extern_inflow`
   from the current `rateinfo` (gross composite inflow, §4), and predict the first
   event (§5/§6).

`setup_network` already enforces Design A (parent = one node subsuming its
subtree) and assigns the `-1` sentinel to any full domain-exiting terminal.

### Touch (`_touch_network!`)

1. **Commit** — advance `state` from `last_solve_time` to `T_commit` with a
   `tmax`-bounded solve (§5) under the *cached* `extern_inflow`. Set
   `last_solve_time = T_commit`. Because it uses the cached inflow, the commit does
   not depend on `rateinfo` still holding pre-event values, so it may run *after*
   `_update_flow!`.
2. **Boundary clamp** — if a trap of this context fired, set its committed volume
   to the post-event boundary per `solveDynNetwork!`'s caller protocol
   (`:fill` → C, `:empty` → 0, `:unspill` → `prevfloat(C)`; kind from
   `ctx.next_event`).
3. **Rebuild** (only if a member trap filled) — `setup_network` with the updated
   `filled_traps`, then **remap** `state` onto the new `net.traps` by global trap
   index (below). Update `global_ix`, `inflow_sources`, `net_trap_set`,
   `net_covered_set`.
4. **Refresh** `extern_inflow = _external_inflow(net, rateinfo, tstruct)` (§4).
5. **Predict** the next event on a *copy* of `state` (§5/§6); store `next_event`
   (absolute) and the `changetimeest` entries.

### State remap across a rebuild

When `setup_network` returns a new `net.traps` ordering/membership, rebuild `state`
by global trap index:

- a trap present in the *old* context carries its committed volume forward (via the
  old `global_ix`);
- a **newly-introduced** trap is seeded from `cur_amounts[trap_ix]`, or from its
  boundary value if it just filled.

This is the only place `state` membership changes; non-rebuild touches keep
`net.traps` fixed.

---

## 4. The external-inflow problem

`solveDynNetwork!` expects `inflow[i]` = inflow to trap `i` from sources **outside**
the network (terrain and non-network traps); intra-network routing is handled
inside the ODE by `_route_flow`. Using `rateinfo.trap_inflow[t]` directly is wrong
when another full network trap is currently spilling to `t` — the ODE would
double-count that contribution.

**Approach: treat network traps as always-empty in the spillgraph.** Network-trap
fill/unfill events are filtered out before `update_spillgraph!`, so network traps
never acquire out-edges. `compute_flow` only routes outflow from traps that have
out-edges, so it skips network traps with no change to `compute_flow`. Effect:

- `rateinfo.trap_inflow[t]` for a network trap receives only terrain +
  non-network-trap contributions.
- No double-counting: the ODE handles intra-network flows; the spillgraph and
  `compute_flow` handle everything else.

**No outflow problem to non-network traps.** `setup_network` follows the terrain
flow graph downstream from each seed, including every trap it meets, and terminates
at the first unfilled trap. So any trap a full network trap would spill into is
already a member. The only true exit is the domain boundary (`spill_path == -1`),
which `_route_flow` drops. No network outflow needs injecting into `rateinfo`.

**Parent nodes — gross composite inflow.** A parent node subsumes all its
descendants, so it is fed the gross composite inflow

```
external_inflow[node] = Σ_{leaf ∈ leaf-descendants(node)} trap_inflow[leaf]
                      = trap_inflow[node]      for a leaf node
```

This matches the solver's parent geometry (it charges the whole composite footprint
infiltration), so no solver change is needed. Using `trap_inflow[P]` directly would
double-count infiltration (it is already net of child infiltration via
`signed_flow = inflow − Smax`, and the solver deducts the composite footprint
infiltration again — review bug B1). `inflow_sources` for a context is the union of
these leaf-descendant trap indices across its nodes.

---

## 5. Bounded integration and the state-commit protocol

**Problem.** `solveDynNetwork!` advances `state` to *its own* next event `T_X`, but
the globally-committed next event may be an earlier non-network event at
`T_Y < T_X`. We must be able to read network volumes at the committed time, not
only at `T_X`.

**API extension — optional `tmax`** (the first implementation slice, §
Implementation order):

```julia
solveDynNetwork!(state, tstruct, net, infiltration, inflow; tmax = Inf, ...)
```

Integration runs over `(0.0, tmax)`. If a topology event occurs at elapsed
`t_e <= tmax`, it is returned as today (with `state` at `t_e`). Otherwise
integration stops at `tmax` and returns `(time = tmax, trap = 0, kind = :none)`
with `state` at `tmax`. `time` is **elapsed** from the solve start; callers add the
absolute start time.

**Two-phase use per context** (keeps `state` at committed volumes, never
over-advanced):

- **Predict** — on a *copy* of the committed state, solve with
  `tmax = endtime − last_solve_time` to find the network's next event. Store
  `next_event.time = elapsed + last_solve_time` and the `changetimeest` entry. The
  committed `state` is untouched.
- **Commit** — once the inner loop knows `T_commit`, advance the committed `state`
  in place with `tmax = T_commit − last_solve_time`, then set
  `last_solve_time = T_commit`. If `T_commit == T_X`, the same call re-confirms the
  topology event.

The commit advance must use the rates valid over `[last_solve_time, T_commit]`. The
context caches those in `ctx.extern_inflow` (refreshed only at a touch), so the
commit is order-independent of the spillgraph/flow update.

---

## 6. `_set_initial_changetime_estimates`

After computing the standard analytical estimates for all traps, predict each
network's first event (on a copy of its state, §5) and override its trap estimates:

```
for each DynNetworkContext ctx:
    ext = _external_inflow(ctx.net, rateinfo, tstruct)         # gross composite, §4
    result = solveDynNetwork!(copy(ctx.state), tstruct, ctx.net, infiltration, ext;
                              tmax = endtime - cur_time, zvt = z_vol_tables)
    ctx.next_event = (; time = result.time + cur_time, result.trap, result.kind)
    if result.kind != :none:
        changetimeest[result.trap] =
            ChangeTimeEstimate(result.kind == :fill, ctx.next_event.time, ctx.next_event.time)
    for every other trap t in ctx.global_ix:
        changetimeest[t] = ChangeTimeEstimate(false, Inf, Inf)
```

Because network-trap estimates are *exact* (`min == max`), the branch-and-bound
resolves them immediately without extra ODE calls — they never reach
`_compute_exact_changetime`.

---

## 7. `_update_changetime_estimates!`

Called at the top of `_identify_next_status_change!` to refresh estimates for traps
whose inflow changed (`getinflowupdates(rateinfo)`).

**Change:** skip any trap in `net_covered_set` (network nodes *and* their subsumed
descendants). Node estimates come from the most recent network prediction; subsumed
descendants are full and static under their parent. `_compute_changetime_estimate`
assumes constant rates and would, e.g., let a subsumed full child spuriously become
a candidate to empty. `net_covered_set` is threaded in as an extra argument; the
candidate selection in `_identify_next_status_change!` must ignore those traps too.

---

## 8. Main event loop (`_fill_sequence_for_weather_event!`)

New steps marked ►:

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
►  inflow_updated = Set(u.index for u in getinflowupdates(rateinfo))
►  for ctx in net_contexts
►      fired          = any(u.index ∈ ctx.global_ix for u in fill_updates)
►      inflow_changed = !isdisjoint(ctx.inflow_sources, inflow_updated)
►      (fired || inflow_changed) &&
►          _touch_network!(ctx, tstruct, dyn_traps, culverts, filled_traps,
►                          cur_amounts, rateinfo, z_vol_tables, cur_time, endtime)
►  end

# amount updates + SpillEvent as before (§9)
```

`_touch_network!` runs the full §3 touch on one context. Networks that neither
fired nor saw an inflow change are skipped entirely — no ODE solve. The touch gate
reads `getinflowupdates` (Decision D4), the same staleness signal `_update_flow!`
just populated.

---

## 9. `cur_amounts` and SpillEvent recording

**Network-trap volumes come only from the solver.** A network trap's amount is
*never* computed with the constant-rate `fill_trap_until`; it is read from the
integrated `ctx.state`. When a context is touched, its commit step leaves
`ctx.state` holding exact volumes at `T_commit`. Therefore:

- `cur_amounts[t]` for a network trap is read directly from `ctx.state` (via
  `global_ix`), at the events where its context is touched.
- `amount_updates` includes every trap of a touched network whose volume changed
  materially, not just the trigger trap. (Untouched networks have unchanged
  dynamics; their amounts are interpolated between touch events — the same
  approximation the existing code already makes for ODE-governed traps.)

The two existing call sites that apply `fill_trap_until` to all affected/stale
traps must exclude `net_trap_set` and take amounts from `ctx.state`:

- `_update_affected_amounts` — skip `net_trap_set`; emit the amount from `ctx.state`.
- the end-of-period finalization loop — see §10.

---

## 10. Weather-event boundary finalization

An ordinary trap's volume at `endtime` is a constant-rate projection
(`fill_trap_until` short-circuits analytically). A network trap follows the
multi-trap ODE and cannot be projected that way, so at the period boundary each
network's committed `state` is advanced to `endtime` with the same bounded commit:

```
solveDynNetwork!(ctx.state, …, ctx.extern_inflow; tmax = endtime − ctx.last_solve_time)
```

The advance is guaranteed event-free (any network event ≤ `endtime` would already
have been processed before the loop broke), so it returns `:none` and leaves
`ctx.state` at exact `endtime` volumes, which become the network traps'
`cur_amounts`. These boundary volumes are the initial state the *next* weather
period rebuilds its networks from (§3), so they must be exact.

**Subtlety:** this final advance must run for **every** network, including ones
never touched during the period (whose `ctx.state` still sits at an earlier
`last_solve_time`) — not only the ones in `net_trap_set` that happened to be
touched.

---

## 11. `_compute_exact_changetime` — no change

Network traps always have `min == max`, so `_compute_exact_changetime` returns on
its existing `changetimes[trap].min == changetimes[trap].max` fast path.

---

## Summary of touched locations

| Location | Change |
|---|---|
| `fill_sequence` signature | Add `dyn_traps` (trap indices), `culverts` |
| `fill_sequence` body | Build `DynNetworkContext`s; maintain `net_trap_set`, `net_covered_set` |
| `solveDynNetwork!` | Add optional `tmax`; return `(tmax, 0, :none)` with `state` at `tmax` when no event fires before `tmax` |
| `_fill_sequence_for_weather_event!` | Filter `fill_updates` from the spillgraph; touch only affected networks (gated on fired-or-inflow-changed via `getinflowupdates`) |
| `_set_initial_changetime_estimates` | After standard estimates, predict each network's first event and override its trap estimates |
| `_update_changetime_estimates!` + candidate selection | Skip `net_covered_set` |
| `_update_affected_amounts` | Skip `net_trap_set`; emit network-trap amounts from `ctx.state` |
| end-of-period finalization | Advance **every** `ctx.state` to `endtime` (§10) instead of `fill_trap_until` |
| NEW `_external_inflow` | Per-node gross composite inflow `Σ trap_inflow[leaf descendants]` |
| NEW `_build_dyn_networks` | Validate `dyn_traps` are valid trap indices; seed `setup_network` with the **union** of each named trap's footprint cell and the culvert endpoints; build contexts |
| NEW `_touch_network!` | Commit one context to `cur_time`, clamp a fired trap, rebuild+remap on a member fill, refresh `extern_inflow`, re-predict |
| `setup_network` / `_subnetwork` (done) | Design A subsume; terminal-detection fix; `spill_path == -1` sentinel |

---

## Settled decisions

- **D1 — Seeding.** Networks are seeded by the union of the explicit `dyn_traps`
  (trap indices, converted to a footprint cell each) and every culvert's
  inlet/outlet cells, formed in `_build_dyn_networks`; `setup_network` is unchanged
  (its existing rule then includes all culverts). Explicit `dyn_traps` are kept
  (not auto-derived) so culvert-free networks are possible for parity testing and
  future non-culvert elements.
- **D2 — Traps, not cells (no hanging paths).** `dyn_traps` are trap indices rather
  than grid cells, so a seed is *always* inside a trap by construction — no
  bare-terrain seed, hence no "hanging" leading path; the only validation is index
  validity. (Culvert endpoints may still be on bare terrain and legitimately seed a
  leading path carrying the culvert's delivered flow.) This makes `path_inflow == 0`
  *exact*: all terrain runoff is already in `trap_inflow`, and the only path-head
  flows are trap spills (`spill_path`) and culvert deliveries (`culvert_outlets`).
- **D3 — Weather boundary.** Network-trap amounts at `endtime` come from advancing
  every `ctx.state` to `endtime` (§10), not from `fill_trap_until`.
- **D4 — Touch gate.** Driven off `getinflowupdates(rateinfo)` ∩ `inflow_sources`
  plus member fills, not recompute-and-compare of `_external_inflow`.

(Constructor decisions already landed: Design A subsume, terminal-detection fix,
`spill_path == -1` sentinel.)

## Deferred (out of scope for this pass)

- **`time_slack`** — unimplemented for the standard path; revisit only after a
  first working integration exists.
- **Uphill / reverse culverts** — `solveDynNetwork!` throws on the resulting cyclic
  network; the fail-loud error is accepted; not pursued for now.

---

## Validation strategy

No pre-made culvert scenario with a known reference trajectory exists, and none is
constructed for this pass. Validation has three legs:

1. **Parity (culvert-free), full *and* mixed coverage.** With no culverts present,
   the dynamic-network path must reproduce plain `fill_sequence`. Two variants:
   - **full coverage** — every trap placed in a network (`dyn_traps` = all traps)
     vs. plain `fill_sequence`;
   - **mixed coverage** — only a *subset* of traps placed in networks, the rest
     left to the analytical path, vs. the same scenario run as plain
     `fill_sequence`.

   The mixed variant is the important one: it exercises the *interplay* between the
   two subsystems (external-inflow handoff §4, spillgraph exclusion, touch gating,
   amount recording) that a full-vs-none comparison cannot reach. Each variant
   compares the **event sequence** (must match exactly: same `(trap, kind)` order)
   and the **timings** (equal to a tight tolerance). Timings are tolerance-based,
   not bit-exact, because the network path integrates ODEs numerically where the
   plain path is analytic; the event sequence/counts can still be checked exactly.
2. **Structural invariants (culvert cases).** Mass conservation, strict event
   ordering, `min == max` change-times for network traps, and the
   three-state / `spill_path` contract.
3. **Visual / animation (qualitative).** A non-CI `examples/` script that animates
   the evolving surface-water state on the 3D terrain and compares runs with vs.
   without culverts by eye. Built from existing helpers — `trap_states_at_timepoints`
   / `interpolate_timeseries` / `_fill_state_to_terrainmap` (`utils.jl`) feeding
   `plotgrid` + `drape_surface` (`IOandplot.jl`, Makie isolated there per AGENTS.md).

---

## Implementation order

1. **`tmax` on `solveDynNetwork!`** — self-contained and independently testable
   (run a network to a known event, then re-run with `tmax` short of it and confirm
   `(tmax, 0, :none)` with `state` at `tmax`). Everything below depends on it.
2. **`_external_inflow`, `_build_dyn_networks`** (trap-index → footprint-cell seed
   union), and `DynNetworkContext` with its prediction/commit helpers.
3. **Changetime overrides** (§6, §7) and the **main-loop touch wiring** (§8).
4. **`cur_amounts` / SpillEvent exclusions** (§9) and **boundary finalization**
   (§10).
5. **Validation**: culvert-free parity first, then structural invariants and the
   animation script.
