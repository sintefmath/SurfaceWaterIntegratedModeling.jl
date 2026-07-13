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
- [ ] B1 `[inv]` enumerate every `DynNBS` use (solver `NBSPlan`/`_build_nbs_plan`/`_nbs_elements`,
  `network_context._nbs_elements`/`_nbs_layer_block`/`_store_nbs_state!`, `build_network`).
- [ ] B2 `[code]` rework `_build_nbs_plan` (+ helpers) to key off `DynNetwork.nbs ::
  Vector{DynNBSPlacement}` and footprint geometry directly — no `DynNBS.placement_ix` handle.
- [ ] B3 `[code]` delete `DynNBS` and `_nbs_elements`; `DynNBSPlacement` is the sole NBS type
  end-to-end. Fix `DynNetwork` constructors / call sites.
- [ ] B4 `[test]` NBS overlay + submergence solver tests still pass (port `nbs_dynamic_test.jl`
  cases as needed).

### Phase C — new network driver (replaces `network_context.jl`'s structural core)
- [ ] C1 `[code]` `_state_for(net, vol_by_trapix, seed0)`: build a component's `state` vector in
  `net.traps` order from a `trap_ix → volume` map, seeding new traps from `seed0` — the single
  re-indexing helper used after every `apply_*`/regrow.
- [ ] C2 `[code]` build entry: one context per component from the new `setup_network`; predict
  each next event (reuse `_predict_network!` / `solveDynNetwork!`).
- [ ] C3 `[code]` per-event dispatch: given the earliest component's `(trap, kind)`, update
  `full_traps`, call `apply_fill!` / `apply_unfill!` / `apply_empty!`, and re-sync the affected
  contexts' `state` via C1 using the returned migrated `trap_ix`.
- [ ] C4 `[code]` state handoff: seed grown/newly-dynamic traps from `cur_amounts`; write
  detached traps back to `cur_amounts`; on `:empty` distribute parent→children (parent 0,
  child `prevfloat(C)`), honouring the solver protocol (§1 table).
- [ ] C5 `[code]` minimal spillgraph/`rateinfo` reconcile for the changed coverage, driven by the
  migrated sets (keep the essence of `_reconcile_spillgraph!`, drop the full recompute).
- [ ] C6 `[code]` finalize: advance each context to `endtime`, write settled node volumes +
  subsumed-descendant capacities to `cur_amounts` (port `_finalize_networks!`).
- [ ] C7 `[test]` driver unit tests on real terrain: a scripted event sequence keeps
  `cur_amounts` consistent and matches an independent recompute.

### Phase D — swap `fill_sequence.jl` call sites, retire old path
- [ ] D1 `[code]` replace L83 `_build_dyn_networks` → new build; L174 `_touch_networks!` → new
  per-event dispatch (C3–C5); L232 `_finalize_networks!` → C6. Remove `SubnetCache` (L124) and
  the `_expand_empty_fill_updates` / `_network_amount_updates` plumbing where the driver subsumes it.
- [ ] D2 `[code]` delete the retired old path: `elements.jl` `setup_network` (L302),
  `setup_network_cached`, `_subnetwork`, `_merge_networks`, `SubnetCache`, and the old
  `network_context.jl` structural functions superseded by C.
- [ ] D3 `[test]` `Sequencing` + `Trapping structure` test sets pass (per AGENTS.md); no
  reference to retired symbols remains (`grep`).

### Phase E — parity + regression tests
- [ ] E1 `[test]` dynamic-vs-analytic parity: on cases with no true multi-trap coupling, the
  dynamic fill/drain event times match `fill_sequence` to `PARITY_TOL`.
- [ ] E2 `[test]` coupled cases (culvert / NBS / shared spill / subsumption / fusion) exercised
  end-to-end through `fill_sequence`, structural invariants asserted.
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
