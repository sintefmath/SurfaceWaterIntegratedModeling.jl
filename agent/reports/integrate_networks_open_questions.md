# Integrate-networks: open questions before implementation

Status snapshot (updated 2026-06-29):

- **B1 constructor fix is done and committed**: `_subnetwork` →
  `_subsume_terminal_parent` (+ `_is_descendant`) in `src/dynamics/elements.jl`,
  with the regression test in `test/dynamics_test.jl`. Design A ("always
  subsume") is settled.
- **`_subnetwork` terminal-detection fix + `spill_path == -1` sentinel done**
  (this session): the changed `flow_path_from` shifted the terminal-unfilled-parent
  case to a trap-terminated chain; `_subnetwork` now handles it, and a full trap
  spilling straight out of the domain carries `spill_path == -1`. Full dynamics
  suite green. See `integrate_networks_plan_review.md` (B2).
- **Plan deferred decisions all settled** (this session): the five items at the end
  of `integrate_networks_plan.md` are now resolved/deferred. Together with the Q4
  decision below, this resolves Q3, Q4, and Q5.
- **No `fill_sequence` integration code exists yet.** `solveDynNetwork!` has no
  `tmax`; `fill_sequence` has no `dyn_traps`/`culverts` args; `_external_inflow`,
  `_build_dyn_networks`, `_touch_network!`, and `DynNetworkContext` don't exist.

**Remaining to address before coding: Q1 (sequencing) and Q2 (touch-gating).**
Q3, Q4, and Q5 are resolved.

The architecture in `integrate_networks_plan.md` reads as coherent and consistent
with the code. Questions below are scoping + one efficiency/consistency point —
not objections.

---

## 1. First slice / sequencing — STILL OPEN

Proposal: land `tmax` on `solveDynNetwork!` first — self-contained,
independently testable (run a network to a known event, then re-run with `tmax`
short of it and confirm `(tmax, 0, :none)` with `state` at `tmax`), and
everything else depends on it.

**Decision needed:** start with `tmax`, or develop the whole integration together
on a branch?

---

## 2. The touch-gating mechanism (main technical question) — STILL OPEN

Plan §8 detects "inflow changed" by recomputing `_external_inflow` for *every*
network on *every* event and comparing `new_ext != ctx.extern_inflow` (exact
float `!=`). That works — because `_update_flow!` updates `trap_inflow`
incrementally, unaffected entries stay bit-identical — but it's
O(networks × nodes) of recompute per event purely to detect change.

`RateInfo` already tracks exactly which traps had inflow changes via its
savepoint mechanism (`getinflowupdates`). Reusing that — *touch a network iff a
leaf-descendant appears in `getinflowupdates(rateinfo)`, or a member fired* — is
both cheaper and exactly consistent with how the rest of the loop decides what's
stale.

**Decision needed:** any reason to prefer recompute-and-compare over driving the
gate off `getinflowupdates`?

---

## 3. `dyn_coords` source (plan open-question 1) — RESOLVED

**Decision:** the argument is a list of **trap indices** (`dyn_traps`), explicit
and caller-supplied (not auto-derived). Each network is seeded from the **union**
of the named traps (converted to a representative footprint cell) and every
culvert's inlet/outlet cells (formed in `_build_dyn_networks`; `setup_network`
unchanged).

Rationale: explicit `dyn_traps` lets us build/run networks with **no culverts**
for parity testing against plain `fill_sequence`, and generalizes to future
non-culvert elements (NBS); culvert endpoints auto-seed so every
culvert-connected trap group is covered without the caller naming it.

Refinement: the original `dyn_coords` (grid cells restricted to footprints) was
replaced by trap indices, since a seed is required to be in a trap anyway. This
makes "no hanging paths" true by construction and reduces validation to a trivial
index-range check. Culvert endpoints may still be on bare terrain. See
`integrate_networks_plan.md` Decisions D1/D2.

---

## 4. Validation fixture — RESOLVED

Tests here are deterministic checksums; the rule is "fix logic, never adjust
expected values."  There is **no pre-made culvert scenario** with a known reference
trajectory to validate timings against, and we will not construct one for this
pass.  Validation has three legs instead:

1. **Parity (culvert-free), full *and* mixed coverage.**  With no culverts, the
   dynamic-network path must reproduce plain `fill_sequence`.  Two variants:
   *full* (every trap in a network) and *mixed* (only a subset networked, the rest
   analytic), each vs. the same scenario run as plain `fill_sequence`.  The mixed
   variant is the important one — it exercises the **interplay** between the two
   subsystems (external-inflow handoff, spillgraph exclusion, touch gating, amount
   recording) that a full-vs-none comparison cannot reach.  Compare the event
   sequence exactly (same `(trap, kind)` order) and timings to a tight tolerance
   (tolerance-based, not bit-exact: the network path integrates ODEs numerically
   where the plain path is analytic).  (Enabled by the Q3 decision to keep
   `dyn_traps` independent of culverts.)

2. **Structural invariants (culvert cases).**  Mass conservation, strict event
   ordering, `min == max` (exact) change-times for network traps, and the
   three-state / `spill_path` contract.  These hold regardless of a reference
   trajectory.

3. **Visual / animation validation (qualitative).**  A non-CI `examples/` script
   that animates the evolving surface-water state on the 3D terrain and lets us
   compare runs **with vs. without culverts** by eye (does the culvert visibly
   relieve/redirect filling as expected).  Building blocks already exist — no new
   core plotting:
   - `utils.jl`: `trap_states_at_timepoints` / `interpolate_timeseries` to
     reconstruct fill state at sampled times from the `fill_sequence` output, and
     `_fill_state_to_terrainmap` to turn a fill state into a water-surface map.
   - `IOandplot.jl`: `plotgrid` to draw the draped 3D surface and `drape_surface`
     to update its texture per frame (assembled with Makie's `record`);
     `set_camerapos` / `grid3Dgeometry` for the view.

   The script stays in `examples/` and uses `IOandplot.jl`'s functions, honoring the
   AGENTS.md rule that Makie-dependent code is isolated to `IOandplot.jl`.

---

## 5. Confirm deferrals — RESOLVED

Outcome of the deferral review (plan open-questions):

- **`path_inflow`** — *not deferred; RESOLVED.* Stays zero and is now **exact**, not
  an approximation: "hanging" upstream paths cannot arise (`dyn_traps` are trap
  indices, so every seed is inside a trap), so all terrain runoff is already in
  `trap_inflow` and the only path-head flows are trap spills and culvert
  deliveries. (plan open-question 2)
- **`time_slack`** — *DEFERRED by decision.* Still unimplemented for the standard
  path; if revisited, only **after** a first working integration exists.
- **uphill / reverse culverts** — *DEFERRED by decision.* `solveDynNetwork!` throws
  on the resulting cyclic network; the fail-loud error is the accepted behaviour;
  not pursued for now.
