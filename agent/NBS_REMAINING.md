# NBS overlay — remaining work & open questions

Status as of the `nbs-dynamic-core` branch (Stages A–C of
`agent/NBS_OPTION1_OVERLAY_PLAN.md` complete). The overlay is wired end-to-end
through `fill_sequence`: static footprint capture → layered storage ODE →
weighted terrain re-emit + piped outlets → NBS→NBS coupling → context
persistence → submergence. dynamics `993/993`; overlay suite `66` tests in
`test/nbs_overlay_test.jl`. Mass conservation is checked at every seam.

This file lists what is **not** done and the deliberate simplifications, so a
reviewer can decide what (if anything) to change before merge.

---

## 1. `nbs_inflow` refresh granularity — **NEEDS REVIEW (flagged by user)**

**Current behaviour.** `nbs_inflow` (the per-placement static footprint capture)
is computed by `watercourses`' sink overlay inside `compute_flow`, which runs
**once per weather period**. Within a period it is held constant — it is a
"per-event constant" on `RateInfo`, deliberately *outside* the incremental
update machinery (`_update_flow!` / `setrunoff!` / the `stored_*` diffs) that
keeps the rest of `RateInfo` (`runoff`, `trap_inflow`) current as traps fire.

**The gap.** A footprint's capture depends on the trap-filled configuration. When
an upstream trap **fills mid-period** and begins spilling, more water reaches a
downstream footprint's catchment — but `nbs_inflow` does not rise to reflect it,
and the incremental `runoff` update (which is NBS-agnostic) routes that new spill
*past* the footprint uncaptured. So the NBS **under-captures** for the rest of
the period.

**Severity.** Accuracy, **not** conservation — no water is created or lost; the
would-be capture is simply routed downstream instead of into the NBS. It
self-corrects at the next weather-period boundary (fresh `watercourses` with the
current `full_traps`). So capture is exact at *period* granularity, stale at
*within-period `:fill`-event* granularity. It bites only when a footprint is
downstream of a trap that (a) fills *during* a period and (b) meaningfully feeds
that footprint's catchment.

**Note.** Plan §7a describes `nbs_inflow` as "piecewise-constant … refreshed at
`:fill`/`:empty` events". The current implementation refreshes it only at
*weather-period* boundaries, not at those events — so it does not yet match the
plan's stated intent. **This is the item the user wants to review before
deciding.**

**Fix options (not implemented):**
1. *Coarse / localized*: on each `:fill` event whose trap is in a footprint's
   catchment, recompute that footprint's capture (a footprint-local `watercourses`
   pass), update `nbs_inflow`, and re-apply the sink to `runoff`. Simple, but a
   partial recompute per relevant event.
2. *Precise*: make the incremental hot path (`_update_flow!` / `_update_runoff!`)
   NBS-aware — terminate-and-tally at footprint cells whenever spill is re-routed.
   Exact, but threads the NBS overlay into the performance-sensitive incremental
   machinery that is currently deliberately NBS-free.

---

## 2. Evapotranspiration (plan §Deferred 4)

The layer model carries `EVCoeff` / `EVS11`, but ET is **not** applied. The rate
function keeps an explicit `ET = 0.0` placeholder in every per-layer balance
(`dS/dt = inflow − infiltration − overflow − ET`), so wiring it in is a one-line
change per layer computed from `EVCoeff`/`EVS11` — not a restructuring. The mass
ledger already reserves the term (do not drop it).

---

## 3. Submerged-feed model — one isolated formula

The submerged subsurface intake uses the **saturated-conduit** model (user
decision): `_nbs_saturated_draw(lp) = compute_outflow(Kinf, ninf, Smin, Smax)·1e-3`
(capacity-limited, drawn from the containing trap). It is isolated in that single
helper so switching to the alternative **flood-head-driven** form
(`Kinf·(flood_depth_mm − Smin)^ninf`) later is a one-line change. See the C-stage
discussion for the trade-off (capacity-limited/bounded vs depth-responsive/
extrapolates the calibrated fit).

---

## 4. Submergence — binary threshold only (plan §5 upgrade path)

Submergence is a **binary** transition at the exit-boundary threshold
(`z_sub = min footprint elevation`): the whole surface block flips dry↔submerged
at once. The plan's deferred upgrade — scaling only the surface layer by a
dry-fraction `(1−f)` from the real per-cell footprint elevations, for a footprint
spanning a large elevation range within one trap — is **not** implemented. Would
be smoother but adds per-cell bookkeeping and tests.

Also within submergence:
- If a solve **starts** already submerged (e.g. the trap was flooded before the
  first solve), no crossing event fires, so the surface block is **not** merged
  into the trap — its water stays frozen in the surface layer (conserved, released
  on emergence) rather than moved into `cur_amounts`. Only the mid-solve *crossing*
  performs the merge. Minor; the frozen water is mm-scale and conserved.

---

## 5. NBS→NBS coupling — chains only, mutual is untested

NBS→NBS coupling is implemented and mass-conserving (element A's re-emit landing
on B's footprint is captured as extra layer-1 inflow to B; both land in one
component). The dynamics are storage-driven, so even a **mutual** A↔B cycle is a
well-posed ODE (no algebraic loop, no ordering needed) — but only **chains** are
tested (`test/nbs_overlay_test.jl` "NBS→NBS coupling"). A mutual/cyclic
configuration is believed fine but has no dedicated test.

---

## 6. Footprint area assumes 1 m²/cell

`_build_nbs_plan` sets each layer's area `A_foot = length(footprint)` (cell
count), i.e. **1 m²/cell** for the mm↔m³ conversion (`S_mm = V·1000/A`). Marked
`@@@` at the point of use; replace with the real cell size once grid resolution is
available (the culvert code has the same assumption).

---

## 7. Not addressed / out of scope here

- **Animation / visualisation** of NBS state over a run (plan pt 5 leg 3) — not
  started.
- **Old carve-based tests** (`test/nbs_test.jl`, `test/nbs_dynamic_test.jl`) remain
  disabled in `runtests.jl`; superseded by `test/nbs_overlay_test.jl`. They could
  be deleted once the overlay is accepted.
- Pre-existing unrelated items: `Base.copy(::RateInfo)` is broken (wrong arity) and
  dead (no callers) — left untouched; reverse/uphill culverts remain deferred.

---

## Test map (`test/nbs_overlay_test.jl`, 66 tests)

| Testset | Covers |
|---|---|
| NBS footprint-as-sink overlay (watercourses) | §1/§2 capture, no-double-count, infil-ignored, multi, overlap-reject, RateInfo plumbing |
| NBS terrain exit weights | natural exit split, funnel/fan, off-domain sentinel, endorheic error |
| NBS overlay element in setup_network | element placed in the right component, downstream pulled in |
| NBS layer ODE + weighted re-emit (rate function) | single/two-layer cascade, weighted delivery, per-layer mass |
| NBS overlay element solved to steady state | NBS-only net settles where overflow = captured inflow |
| NBS→NBS coupling | A→B chain, one component, capture, mass (no leak) |
| NBS through fill_sequence | end-to-end delay (upstream non-submerging NBS), empty-NBS byte-identical |
| NBS submerged rate regime | surface frozen, slots zeroed, mass with recirculation |
| NBS submergence (detect + merge + drain) | derivation, closed-loop mass, mid-solve crossing → merge → :fill |
