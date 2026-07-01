# First-principles note: multi-level parent draining (the late-unspill discrepancy)

## Symptom

Two-weather drain (rain `t∈[0,0.5)`, drought `t≥0.5`), `dyn_traps=[414]`, no culverts.
Fill times match plain `fill_sequence` to ~4e-5, but **drain/expose times are off by
~0.37**: trap 414 is exposed at ~0.765 in dyn vs. unspilling at ~0.549 in plain.

414's ancestry is **414 ⊂ 443 ⊂ 455** (own volumes 0.109, 0.039, 0.430). While 455/443
are full, 414 is *subsumed* under the topmost full node 455 (Design A), so 414 is not
draining late — it is **exposed late**, because the composite 455 (then 443) drains too
slowly.

## The physics (what the correct rate is)

When 455 holds water in its own volume (above the children's spillpoints), that water
covers the whole footprint (455 has **no exclusive cells**: `I_455excl = 0`; the own
volume is a layer on top of the full children). So every footprint cell is wetted and
infiltrating. With external inflow `Q_in` into the composite and whole-footprint
infiltration `I_whole`, the own volume changes at

    dV_own/dt = Q_in − I_whole.

The water cascades: `Q_in` arrives at the leaves → leaves infiltrate `I_leaf` and overflow
the rest → 443 infiltrates `I_443excl` and overflows → 455 infiltrates `I_455excl` and
holds the remainder. Net into 455's own volume = `Q_in − I_leaf − I_443excl − I_455excl =
Q_in − I_whole`. So the **correct own-volume drain rate is `Q_in − I_whole`** (here
`I_whole = 2.3`, `I_455excl = 0`).

## Both formulations *should* agree (the algebra checks out)

- **dyn (solver, §4):** `external_inflow[455] = Σ_leaf trap_inflow[leaf]` (gross composite
  inflow), loss `= footprint_infil[455] = wetted_infiltration(cap) = I_whole = 2.3`
  (verified). So `dyn_net = Σ_leaf trap_inflow − I_whole`.
- **plain (`fill_trap_until`):** `trap_inflow[455] − (Smax[455] − Smin[455])`, where
  `trap_inflow[455]` is the *cascade* overflow `Σ_children (trap_inflow[child] −
  Smax[child])` and `Smax[455] − Smin[455]` is the exclusive-cell infiltration.

Writing `S = subtrap infiltration deducted on the cascade leaves→455` (= leaf + 443-excl
infiltration), one gets `external_inflow = trap_inflow[455] + S` and, for a clean
hierarchy, `Smin[455] = S` so `I_whole − (Smax−Smin) = S`. Substituting:

    dyn_net = (trap_inflow[455] + S) − I_whole
            =  trap_inflow[455] − (I_whole − S)
            =  trap_inflow[455] − (Smax[455] − Smin[455]) = plain_net.

**So this is NOT a fundamental formula error** — for a static full composite the gross
(dyn) and cascade (plain) accountings are algebraically identical. The own-volume drain
rate is the same `Q_in − I_whole` on both sides.

## Where the discrepancy actually is

Since the *rate* agrees, the ~0.37 gap is in the **drain TRANSIENT / unspill timing** —
specifically *when* the composite's feed `Q_in` decays below `I_whole` during the drought.
That feed comes from upstream (non-network) traps that were full at drought start and
drain over `[0.5, ~0.8]`. The two paths differ in HOW that feed reaches the composite:

- **plain:** the upstream spill cascades through the chain of *full children*
  (leaf → 443 → 455); each link must be full to pass it on, and each deducts its
  infiltration as it goes.
- **dyn:** the upstream spill is deposited into `trap_inflow[leaf]` (the leaves are
  *covered*, so in the reconciled spillgraph they have no out-edge and the flow stops
  there), and `external_inflow[455] = Σ_leaf trap_inflow[leaf]`.

For a *static* configuration these match (shown above). The suspicion is that they
diverge **dynamically** during the drain because:
1. the covered leaves short-circuit the cascade (dyn sums at the leaves; plain routes
   through the full-children chain), so the transient feed seen by 455 differs; and/or
2. the spillgraph **reconciliation** during the shrinking-coverage drain produces
   `trap_inflow[leaf]` values that lag/differ from plain's cascade.

Net effect observed: dyn's `external_inflow[455]` stays above `I_whole` longer than
plain's effective feed → 455 unspills late → descendants exposed late.

## Recommended next step (one clean measurement, then fix the lagging side)

The algebra says the steady-state rate is right, so do NOT change the §4 gross-inflow
formula or `footprint_infil` blindly. Instead, instrument a single drought event and dump
side by side:

- dyn: `external_inflow[455]` and `footprint_infil[455]`;
- plain: `trap_inflow[455]` and its `fill_trap_until` infiltration term;
- the per-leaf `trap_inflow[leaf]` in both runs at the same wall-time.

If the per-leaf `trap_inflow[leaf]` already differ → the bug is in the **reconciliation /
flow accounting for covered leaves during the drain** (coverage-exit side), and the fix
belongs there (mirror of the filling-side double-count fix). If the per-leaf values agree
but `external_inflow[455] − I_whole ≠ trap_inflow[455] − (Smax−Smin)` → the bug is in the
gross-vs-cascade composition for a *multi-level* node (e.g. `lowest_subtraps_for` vs.
immediate-child Smax bookkeeping), and the fix is in `_external_inflow`.

Either way it is a localized accounting fix in the network/spillgraph layer, not a change
to the physical model: **the target rate is `Q_in − I_whole`, and both formulations encode
it correctly for a static composite — the task is to make the draining transient honor it.**
