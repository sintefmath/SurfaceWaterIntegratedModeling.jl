# Gate: wiring the dynamic-membership layer into fill_sequence

Scoping of the "gate" — driving the incremental network layer
(`network_updating.jl`: `apply_fill!` / `apply_unfill!` / `apply_empty!` on the new
`build_network.jl` representation) from the live `fill_sequence` solve.

**Caveat (from the maintainer):** the *existing* dynamic integration
(`network_context.jl`, and the `_touch_networks!` / `setup_network_cached` path in
`fill_sequence.jl`) was never in active use, not properly reviewed, and is likely wrong
in several details. Treat it as a **sketch of the shape**, not a correct implementation —
salvage ideas, not code, and re-derive correctness.

## 1. What already exists (and how much is reusable)

- **The solver `solveDynNetwork!` (`networksolver.jl`) is in good shape and its event model
  is exactly ours.** It evolves per-component trap volumes and halts on one of three
  topology events, returning `(time, trap, kind)`:

  | solver `kind` | meaning | our entry point |
  |---------------|---------|-----------------|
  | `:fill`    | trap reaches capacity, starts spilling | `apply_fill!` (grows spill; regrows on subsumption/fusion) |
  | `:unspill` | a full trap's net inflow goes negative, stops spilling | `apply_unfill!` (detach) |
  | `:empty`   | a parent drains to its floor, exposes its subtraps | `apply_empty!` (de-subsume) |

  The solver's documented "caller protocol" already says, per event, to update
  `full_traps` and *"rebuild network via setup_network"* — that whole-rebuild step is what
  we **replace** with the incremental `apply_*` call. Everything else in the protocol
  (state clamping to `C` / `prevfloat(C)` / `0`, exposing children on `:empty`) is the
  state-handoff spec we must honour.

- **The rate / distributor layer (routing, culverts, NBS overlay, submergence) is largely
  built** in `networksolver.jl`. The "distributor re-impl" is smaller than feared — but it
  is currently built around the **`DynNBS` overlay element** (`_nbs_elements`, `NBSPlan`,
  `DynNBS.placement_ix`). Retiring `DynNBS` (deriving NBS coupling from `DynNBSPlacement`
  footprint geometry, per the routing redesign) is the real distributor work.

- **`network_context.jl` (the integration layer) is the main thing to rewrite.** Its
  *concepts* are reusable — `DynNetworkContext` (per-component `state` + `global_ix` +
  `last_solve_time` + `extern_inflow`), commit/predict, weather-boundary finalize, the
  `:empty` → expose-children handling (`_expand_empty_fill_updates`,
  `_apply_fired_boundaries!`), amount updates back to `cur_amounts`. Its *mechanism* — full
  structural retrace every event via `setup_network_cached` / `SubnetCache` +
  `_reuse_plan` gating + `_reconcile_spillgraph!` — is what the incremental `apply_*`
  layer supersedes.

## 2. Target architecture

`fill_sequence` keeps: `cur_amounts` (per-trap committed volumes/times), the spillgraph +
`rateinfo` flow, weather periods. A **network driver** owns `Vector{DynNetwork}` (the
components) with one `state`/context per component and:

1. **Build** (period start): new `setup_network(tstruct, full_traps; dyn_coords, culverts,
   nbs)`; one context per component; predict each component's next event.
2. **Per event** (from the earliest component prediction): dispatch on `kind` —
   `:fill → apply_fill!`, `:unspill → apply_unfill!`, `:empty → apply_empty!` — passing the
   event trap and the updated `full_traps`. Each returns the migrated `trap_ix`.
3. **State handoff** using the migrated sets and the solver's clamping rules:
   - grown/newly-dynamic traps: seed `state` from their static `cur_amounts` level;
   - detached traps: write their `state` back to `cur_amounts` (→ static handling);
   - `:empty` de-subsumption: distribute the parent's volume across the exposed children
     (parent → 0, each child → `prevfloat(C_child)` per the solver protocol);
   - re-index each touched component's `state` vector to the new `net.traps` order.
4. **Reconcile** the spillgraph/`rateinfo` for the changed coverage (a covered trap's
   outflow is ODE-routed, not spillgraph-routed) — the one piece of `_reconcile_spillgraph!`
   worth keeping, now driven by the `apply_*` migrated sets instead of a full recompute.
5. **Finalize** at the period boundary: advance each context to `endtime`, write settled
   volumes to `cur_amounts` (nodes) and capacity to subsumed descendants.

## 3. Work items (suggested order)

**A. Land `build_network.jl` in the module.** Add the include; resolve the `setup_network`
name collision with the old `elements.jl` one (they have different signatures today — retire
the old `setup_network` / `_subnetwork` / `setup_network_cached` / `SubnetCache` path once
nothing calls it). Remove the `@eval` injection from `test/dynamic_membership_test.jl`.

**B. Distributor: retire `DynNBS`.** Rework the solver's NBS plan (`_build_nbs_plan`,
`_nbs_elements`) to key off `DynNBSPlacement` + footprint geometry directly (the geometric
coupling we settled on), so `DynNetwork.nbs::Vector{DynNBSPlacement}` is the single NBS
representation end-to-end. This unblocks removing `DynNBS`.

**C. New network driver (replaces `network_context.jl`'s structural core).** Build via the
new `setup_network`; per-event `apply_*` dispatch; state handoff (§2.3); minimal
spillgraph reconcile (§2.4); finalize (§2.5). Keep/port the reusable context concepts;
drop `setup_network_cached` / `_reuse_plan` / full-retrace.

**D. `fill_sequence.jl` call sites.** Swap `_build_dyn_networks` / `_touch_networks!` /
`_finalize_networks!` for the new driver; delete the retired machinery.

**E. Tests.** Parity tests: the dynamic path's fill/drain event times vs the analytic
`fill_sequence` on cases with no true multi-trap coupling (should match to `PARITY_TOL`),
plus coupled cases (culvert / NBS / shared spill) exercised structurally. Extend
`test/dynamic_membership_test.jl` and `dynamics_test.jl`.

## 4. Open questions

- **How much of `network_context.jl` to port vs rewrite.** Given the "unreviewed, likely
  wrong" warning, default to rewrite of the structural core, cherry-picking only the
  state-handoff arithmetic after re-deriving it. Confirm the `DynNetworkContext` field set
  is still the right state to carry.
- **`apply_empty!`'s trigger.** The solver already emits `:empty` for a parent at its floor
  with negative net inflow (and a t=0 fast path) — that *is* the mid-level trigger we
  called gate-coupled. Confirm the exposed-children volumes and `full_traps` bookkeeping
  match `apply_empty!`'s post-regrow expectation.
- **State re-indexing on regrow.** After a fusion/subsumption regrow the component's
  `net.traps` order changes; the `state` vector must be rebuilt from per-`trap_ix`
  committed volumes. Define one helper (`state for the new net from a `trap_ix → volume`
  map + seeded new traps) and use it everywhere.
- **Coupling detection cost per event.** `apply_fill!` rebuilds `_index_components` each
  call (O(all live cells)); fine at fill-sequence event rates, revisit only if measured.
- **`DynNBS` removal blast radius.** Enumerate every `DynNBS` use (solver plan, context,
  tests) before B.

## 5. Non-goals / deferred

- Fission stays deferred (`_compact!` cleans indices; over-merged components remain correct).
- Reverse/uphill culverts remain out of scope (the solver already fails loud on a cycle).
- Perf tuning of regrow / index rebuild — correctness and wiring first.

## 6. Detailed action checklist

Ordered; run the relevant tests after each phase. `[inv]` = investigation, `[code]`,
`[test]`, `[doc]`. Key call sites found: `fill_sequence.jl` invokes the driver at L83
(`_build_dyn_networks`), L124 (`SubnetCache`), L140 (`_expand_empty_fill_updates`), L174
(`_touch_networks!`), L203 (`_network_amount_updates`), L232 (`_finalize_networks!`). Old
path to retire: `elements.jl` `setup_network` (L302), `setup_network_cached` (L443),
`_subnetwork` (L594), `_merge_networks` (L765+), `SubnetCache`. `DynNBS` lives in
`network_context.jl`, `networksolver.jl`, `build_network.jl` (17 sites).

### Phase A — land `build_network.jl` in the module (both `setup_network`s coexist) — DONE
- [x] A1 `[code]` added `include("dynamics/build_network.jl")` after `network_updating.jl`; the
  new `setup_network(tstruct, full_traps; …)` and old `setup_network(tstruct, dyn_coords,
  full_traps; …)` coexist by arity (no ambiguity).
- [x] A2 `[test]` dropped the `@eval … include(build_network.jl)` shim from
  `test/dynamic_membership_test.jl`; suite passes standalone (280/280).
- [~] A3 skipped — the `.jl~` / `#…#` / `.#…` files are the maintainer's live emacs
  lock/autosave (untracked, don't affect the build); leave them.
- [x] A4 `[test]` package precompiles clean; `dynamic_membership_test.jl` green (280/280).

### Phase B — distributor: retire `DynNBS` (geometric NBS coupling)
- [x] B1 `[inv]` enumerated `DynNBS` uses — solver (`_build_nbs_plan`, `NBSPlan`), the OLD
  `elements.jl` path (`_expand_with_nbs`/`_nbs_union_edges`/`_build_component`/`_merge_networks`),
  `network_context._nbs_elements`, tests. Found the layer half-migrated: `DynNetwork.nbs` is
  `Vector{DynNBSPlacement}` but the solver read `nb.placement_ix` — a `DynNBS` field it doesn't
  have (broken). See `agent/networksolver_callgraph.md`.
- [x] B2 `[code]` reworked `_build_nbs_plan` (+ `_build_rate_params` call) to read
  `nb.system.layers` and `nb.id` straight off `net.nbs`; `setup_network` now stamps
  `nb.id = position` (the `nbs_inflow` key). Verified on a hand-built net + solve: NBS layer
  fills from the id-keyed `nbs_inflow`, no `placement_ix` error.
- [x] B3 `[code]` deleted `DynNBS` (struct + docstring) and `_nbs_elements`; `DynNBSPlacement`
  is now the sole NBS type end-to-end. Gutted the OLD-path NBS wiring: `_expand_with_nbs`,
  `_nbs_union_edges`, and the `nbs` args on `setup_network` / `setup_network_cached` /
  `_merge_networks` / `_components` / `_build_component` all reverted to culvert-only forms;
  `DynNetwork`'s 3-arg compat ctor now defaults to `DynNBSPlacement[]`. Context call sites
  dropped `nbs=`. Removed the DynNBS-dependent testsets from `nbs_overlay_test.jl` (kept the
  two live ones — watercourses footprint-sink + terrain exit weights) and deleted
  `nbs_dynamic_test.jl`. Module loads; membership + overlay suites 303/303.
- [x] B4 `[test]` NBS end-to-end verified: over a footprint sweep, `setup_network(...; nbs=…)`
  builds NBS-carrying components and `solveDynNetwork!` runs (13/13 solved, NBS layer fills
  from the id-keyed `nbs_inflow`). No committed test yet (a proper one needs a curated
  coupling+draining footprint) — add with C's driver tests.

**Fixed during B2 — three pre-existing new-path NBS *build* bugs** (the NBS build had never run
with a real placement): (1) `build_network._downstream_cell` was shadowed by `utils.jl`'s
more-specific typed `_downstream_cell(::SimpleDiGraph,::Int)` (returns `(cell, terminus)`) —
deleted the duplicate, use the canonical one; (2) `Vector(::Set)` → `collect`; (3)
`inv_flow[cell]` KeyError → `get(…, Int[])`. `setup_network` now builds NBS components cleanly.

### Phase C — new network driver (replaces `network_context.jl`'s structural core)
All in the new `src/dynamics/network_driver.jl` (`NetworkDriver` owns `comps` +
`vol_by_trapix` + `nbs_state`; contexts rebuilt off the mutated component set each event).
- [x] C1 `[code]` the shared seed rule (`full→C`, committed, else project) is `_driver_state0`,
  a per-trap closure `_make_context` consumes directly (used by build and every rebuild). The
  earlier standalone `_state_for` vector form was redundant and has been removed.
- [x] C2 `[code]` `build_network_driver`: components from the new `setup_network`, one context
  each via `_make_context`, `_predict_network!` for the next event. `_driver_next_event` picks
  the earliest.
- [x] C3 `[code]` `step_network_driver!` dispatch: earliest `(trap, kind)` → `apply_fill!` /
  `apply_unfill!` / `apply_empty!`, then `_rebuild_contexts!` off the mutated `comps` (returns
  the migrated `trap_ix`).
- [x] C4 `[code]` state handoff: commit all contexts to the event time BEFORE the structural
  change (§7.3), `_harvest!` into the stores, clamp per protocol (`:fill`→C, `:unspill`→
  `prevfloat(C)`, `:empty`→parent 0 / children `prevfloat(C)` with children out of `full_traps`
  first, §7.2); grown traps seeded via the projection rule.
- [ ] C5 `[code]` spillgraph/`rateinfo` reconcile for the changed coverage — **deferred to D1**,
  where `filled_traps`/`sgraph`/`rateinfo` are threaded from `fill_sequence` (the driver mutates
  `filled_traps`; the reconcile belongs at the call site alongside `_reconcile_spillgraph!`).
- [x] C6 `[code]` `finalize_network_driver!` reuses `_finalize_networks!` (settles nodes +
  subsumed descendants into `cur_amounts`, carries NBS state forward).
- [~] C7 `[test]` `test/network_driver_test.jl` (14 assertions): build/predict,
  `_driver_next_event`, and a stepped single-basin evolution (time-ordered events + finalize at
  capacity) on synthetic terrain. Multi-trap absorption parity needs `fill_sequence`'s static
  loop to keep `cur_amounts` current (isolated, a grown trap's projection over-fills) → that
  parity lands in **E**. The driver's **NBS path is still untested** (no NBS in the C7 fixture)
  → first exercised in **E**.

**Cleanup after C (commit `340f976`).** Fixed a B3-induced defect: `_nbs_layer_block` /
`_store_nbs_state!` read `nb.placement_ix`, gone from `DynNBSPlacement` (now `.system` / `.id`);
the OLD path masked it (empty `net.nbs`) but the driver's NBS path would have thrown. That
removed the need to thread a placement list: `solveDynNetwork!`'s `nbs_placements` arg is unused
since B2, and `DynNetworkContext.seeds` was read nowhere — both context fields dropped,
`_make_context` slimmed to `(net, tstruct, rateinfo, state0, cur_time, nbs_state)`. `NetworkDriver`
dropped the redundant `nbs_placements` / `culverts` / `dyn_coords` fields (the components are the
structural source of truth). Suites 317/317.

### Phase D — swap `fill_sequence.jl` call sites, retire old path
- [x] D1 `[code]` driver wired into `fill_sequence` behind a `use_driver` kwarg, **alongside** the
  old path (reversible). Build → `build_network_driver`; touch → `_touch_networks_driver!` (the
  driver equivalent of `_touch_networks!`: commit → `apply_*` → **C5 spillgraph/`rateinfo`
  reconcile** → rebuild); finalize reuses `_finalize_networks!` on `driver.contexts`. Both paths
  expose the same `net_contexts` interface, so `_expand_empty_fill_updates` /
  `_network_amount_updates` / changetime plumbing are shared unchanged. `SubnetCache` (old-only)
  left until D2. Verified on mini.txt: `use_driver` off/on identical without a network (462/462,
  0.0 diff); a `dyn_traps=[233]` multi-trap chain runs end-to-end through the driver (462 events,
  monotone, all amounts finite & within capacity). `Sequencing` grid1/grid3 (default path) still
  pass. New test: `network_driver_test.jl` "through fill_sequence (D1)". Multi-trap absorption —
  the C7 caveat — now works (the static loop keeps `cur_amounts` current).
  - **Fixed a `build_network` tracer bug D1 surfaced**: when the trace stops at a FULL trap that
    spills straight out of the domain (no distinct downstream cell), `_grow_network_from_seed!`
    left its `spill_path` at 0 (reads as an unfilled frontier) instead of the `-1` out-of-domain
    sentinel → `solveDynNetwork!` three-state-contract error. Now sets `-1`. Reproduced by a fresh
    `setup_network` (not an `apply_*` divergence). Suites 317/317 after the fix.
- [x] D2 `[code]` retired the old path entirely. `fill_sequence` now builds the driver
  unconditionally (the `use_driver` flag is gone); removed the old build/touch branch,
  `SubnetCache`, and the unused `dyn_traps`/`culverts` params from the inner loop. Truncated
  `elements.jl` to the struct definitions only — deleted OLD `setup_network` / `_subnetwork` /
  `_subnet_deps` / `_merge_networks` / `_combine_subnets` / `_dedup_traps` / `_culvert_owners` /
  `_resolve_cell_overlaps!` / `_components` / `_build_component` / `_topological_order` /
  `_expand_with_culverts` / `_occupied_cells` / `setup_network_cached` / `SubnetCache` /
  `_build_network` / `_unfilled_*` / `_subsume_terminal_parent` / `_is_descendant` /
  `_spills_out_of_domain` (~780 lines). Deleted the superseded `network_context.jl` structural
  core (`_build_dyn_networks`, `_touch_networks!`, `_reuse_plan`, `_assemble_contexts`,
  `_contexts_to_commit`, `_commit_contexts!`, `_affected_contexts`, `_trap_owner_map`,
  `_clamp_full_traps!`). A `grep` confirms no live references to any retired symbol remain.
  Discovered the old path was already broken (`_subnet_deps` referenced an out-of-scope
  `ends_with_path`) — it never ran, confirming there was nothing to preserve.
- [x] D2-tests `[test]` ported `dynamics_test.jl` to the retirement: an `mk_network` shim forwards
  the old positional `setup_network(ts, coords, full)` to the new keyword form; deleted the
  old-internals testsets (`_merge_networks`/`_components`/`_build_component`/`_combine_subnets`/
  `_unfilled_trap_at`/`_topological_order`/…) and the slice-2 `_build_dyn_networks` set (now
  covered by `dynamic_membership_test` + `network_driver_test`); dropped 3 culvert-inclusion
  subtests obsolete under the new seed-all-outlets model; fixed 5-arg `DynFlowPath` calls to the
  6-arg form and added a `srcpath` helper for empty-cells source paths.
- [x] D3 `[test]` `Sequencing` grid1/grid3 pass (default path unchanged); the combined dynamics
  suite is 1046 pass / 0 fail / 3 broken. (`Trapping structure` grid2 is the pre-existing bay
  drift, unrelated.)

### Phase E — parity + regression tests
- [x] E1 `[test]` dynamic-vs-analytic parity. **Full-coverage bug FIXED** (was the last blocker).
  Root cause was NOT `grow_spill!` (the earlier guess): a full-coverage / heavy-rain instant
  fills several traps at once; the b&b fires one and commits the others to `cur_time` at exactly
  their capacity `C`. The seed then handed the solver a *transitory* frontier node (`spill_path ==
  0`) at `V == C` → three-state-contract error. Fix in `_driver_state0`: clamp a NON-full trap's
  seed to `prevfloat(C)`, so it stays transitory and its `:fill` fires on the next predict
  (serialising simultaneous fills; mass preserved to a ULP). Full coverage now matches plain on
  mini.txt: 462 events, same filled set, max fill-time drift 2.98e-5 < `PARITY_TOL`. Single /
  mixed / subtrap-seed parity already passed. Restored the `dynamics_test.jl` slice-3 assertions.
- [~] E2 `[test]` coupled cases (culvert / NBS / shared spill / subsumption / fusion) exercised
  end-to-end through `fill_sequence`, structural invariants asserted.
  - **Culvert coupling FIXED.** Root cause of the directional-effect gap: the new `setup_network`
    seeded culvert *outlets* only, so a culvert's *inlet trap* was not an evolving node — it
    filled statically and the culvert never drew from it (the old `_expand_with_culverts` seeded
    BOTH endpoints). Fix in `build_network.jl`: seed `culvert.inlet` as well. Now the culvert
    couples both traps — on mini.txt a large-bore culvert drains its inlet (233) below capacity
    so it never fills and delivers to its outlet (13), which fills far earlier. Verified
    mass-conserving (the ~0.1 stored-water difference is exactly the drained inlet's capacity
    leaving via 13's out-of-domain spill — nothing created/destroyed). Rewrote the
    `culverts through fill_sequence` directional + mass assertions to the correct coupled
    behaviour (no more `@test_skip`). Combined dynamics suite **1053 pass / 0 fail / 0 broken**.
  - **Still TODO:** NBS through the driver (no NBS fixture yet); shared-spill / subsumption /
    fusion end-to-end coupled cases.
- [ ] E3 `[doc]` update `AGENTS.md` (subsystem table), `DYNAMIC_MEMBERSHIP_PLAN.md` §8 (gate
  done), and the memory status once `build_network` is in the module.

### Open-question spikes (do before the phase that needs them)
- [x] S1 `[inv]` diff the old `network_context.jl` state-handoff arithmetic against the solver
  protocol — see §7.
- [x] S2 `[inv]` confirm `:empty` exposed-children volumes + `full_traps` deltas match
  `apply_empty!` — see §7.

## 7. Spike findings (S1, S2)

**S2 — `:empty` ↔ `apply_empty!` (confirmed match).** The solver's `:empty` condition is
`V[parent] → 0`; `affect!` reports `(kind=:empty, trap=parent_global_ix)`. Protocol after it:
clamp `state[parent]=0`, **remove the parent AND all its immediate children from `full_traps`**,
set each `state[child]=prevfloat(C_child)`, rebuild. That maps exactly onto `apply_empty!`:
- the driver removes the children from `full_traps` first, then calls
  `apply_empty!(comps, tstruct, full_traps, parent)`, which regrows the parent's component;
  with the children no longer full they de-subsume into separate nodes (parent disappears as a
  node — matching the tested cycle).
- **Critical ordering**: children must leave `full_traps` *before* the regrow, else it
  re-subsumes them under the parent at `V=0` and `:empty` re-fires forever. The old
  `_expand_empty_fill_updates` already encodes exactly this hazard — a reusable insight.
- Multi-level: only *immediate* children are exposed; deeper still-full descendants stay
  subsumed, which the regrow reproduces for free. `apply_empty!` needs no change.

**S1 — old state-handoff arithmetic is largely correct and reusable; structural core is not.**
Reusable, and matching the protocol (§1 table):
- `_apply_fired_boundaries!` — `:empty`→parent 0 / children `prevfloat(C)`; `:unspill`→`prevfloat(C)`. ✓
- `_clamp_full_traps!` — full trap pinned to exact `C` (kills ODE drift). ✓
- `_assemble_contexts.state0` — full→`C`, committed→committed, else `project` from the static
  path (`fill_trap_until` use_saved) — the seed for grown/newly-absorbed traps. ✓ (concept)
- `_finalize_networks!` — advance to `endtime`; node→settled volume, subsumed descendant→`C`. ✓
- `_external_inflow` / `_inflow_sources` — per-node inflow summed over **leaf** descendants
  (robust to `_reconcile_spillgraph!` withdrawing covered children's flow). ✓
Replaced by `apply_*` (do NOT port): `setup_network_cached` retrace, `_reuse_plan` gating, the
full-recompute `_reconcile_spillgraph!` diff, `_touch_networks!` orchestration.

Three details to **re-verify while building C** (the "likely wrong" risk lives here):
1. **Recompute `covered`/`inflow` from the post-`apply_*` components**, don't diff — the driver
   must read coverage/external-inflow off the mutated component set, not an incremental delta.
2. **Grown-trap state projection** (`project(g)`) is the stale-state absorption hazard
   ([[stale-state0-network-absorption-bug]]); seed a newly-absorbed trap from its value
   projected to `cur_time`, not its stale `cur_amounts` entry.
3. **Commit each affected component to `cur_time` under its cached `extern_inflow` BEFORE
   applying the structural change** (order-independence); mutate topology only on the
   committed state.
