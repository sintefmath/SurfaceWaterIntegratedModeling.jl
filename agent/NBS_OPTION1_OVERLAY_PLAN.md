# NBS as an infiltration/interception overlay (Option 1) — SELECTED

**Status: chosen approach.** The interior-region model
(`NBS_INTERIOR_REGION_PLAN.md`) is the rejected alternative, kept for reference.
Motivation and option analysis: `agent/prompts/nbs_integration_reconsideration.org`.

## Core idea

The NBS is **not** a trap and does **not** touch `spillanalysis`. Terrain traps
are determined by real geometry, independently. An NBS is a
capture-and-re-emit **overlay on the runoff field**: it intercepts all flow
entering its footprint, runs the lumped layer model, and re-emits only the
top-layer overflow onto the terrain (plus bottom-layer overflow to a piped
outlet). Structurally identical to how a *full trap* already rewrites the runoff
field (`_compute_trap_output` + `_update_runoff!`, `watercourses.jl:86-113`) —
so it slots into existing machinery, not new routing.

## Why this over the interior-region model

The interior model makes the NBS "a trap that doesn't fill by volume" — a hybrid
that fights the trap machinery at the merge/submerge seam:
- forces unphysical **ponding** on slope footprints (permeable pavement) via the
  `-4` flag;
- needs a **submerge/emerge state handoff** on the absorption seam that already
  produced two stale-state bugs;
- needs a **scheduler fill-time reconciliation** (layer storage vs terrain
  volume);
- concentrates **runoff** at one terrain-imposed spillpoint.

The overlay avoids all four: real traps stay real; the NBS is a pure overlay
that reads local water depth. No forced ponding, no hierarchy handoff, natural
terrain-spread runoff.

Cost it takes on instead: computing footprint inflow/outflow as a runoff-field
edit (below), and determining which downstream traps must join the dynamic
network. Both use existing patterns.

## Key semantics: `runoff` is accumulated

`watercourses` propagates `runoff[ds] += runoff[cur]` down the flowgraph
(`watercourses.jl:72`). So `runoff[c]` is the **accumulated** intensity at `c` —
downstream cells already contain upstream flow. Summing a cell and its
downstream neighbor double-counts. The inflow/outflow formulas below are written
to avoid this.

## 1. Inflow to the footprint `F` — footprint-as-sink overlay (static)

Realized as a **sink overlay in the `watercourses` accumulation sweep**, so it is
static per trap-config and independent of NBS layer state (it does *not* make
rateinfo depend on anything dynamic). Every footprint cell is treated as a sink:
when the topological sweep reaches it, its accumulated `runoff` **terminates**
(no downstream propagation) and is **tallied** into a per-placement accumulator.

- an external streamline entering at cell `v∈F` terminates at `v`, tallying
  `runoff[v]` (its upstream inflow + its own precip);
- interior footprint cells tally their own precip;
- `zero terrain infiltration over F` first — the NBS layer model *is* the
  infiltration, feed it gross water.

Summed over `F`: NBS inflow `= inflow_external + direct_precip`. This equals the
edge-sum `Σ runoff[u]` over boundary edges `u→v, u∉F, v∈F` plus footprint precip,
but is obtained for free from the existing accumulation machinery — no separate
edge bookkeeping. `watercourses` returns one extra `nbs_inflow` vector (per
placement); `rateinfo` carries it; the NBS element reads it as its external
inflow.

## 2. Capture (same overlay) — removes the downstream double-count

Terminating at footprint cells means the footprint throughput **does not
continue downstream**, so a downstream trap `T` no longer double-counts the
footprint's catchment (the reason the sink overlay is needed at all, vs. a pure
read). On a non-convex footprint, terminate-at-entry also makes exit-and-re-enter
impossible, so each external streamline is tallied once.

## 3. Outflow — re-emit only overflow

The overflow is layer-storage-driven (`compute_outflow`), hence **dynamic** — it
lives in the network solver, NOT in `watercourses`, reusing the existing
outlet-delivery path (`nbs_actual` → `:nbsout` / trap outlets). Routing rule
(decision 1): **the `n` topmost layers re-emit at terrain; every outflowing layer
below them requires an explicit piped outlet.** Each layer's outflow is delivered
to a resolved outlet cell:
- **terrain re-emit** (top `n` layers): **distributed** across the footprint's
  lower-edge exit cells (`w` in each edge `v→w, v∈F, w∉F`), weighted by the flow
  fraction that would naturally exit there *without* the NBS. Magnitude is the
  layer's dynamic overflow; the **weights are static**, precomputed once per
  placement by a footprint-local accumulation pass (the fraction of footprint
  throughput exiting at each edge cell). Delivery: `overflow × weight_i` to each
  exit cell `i`, each outside `F` (so not re-swallowed by the sink overlay §1),
  then propagating downhill through normal routing. Natural terrain spread — no
  single imposed spillpoint; a funnel footprint whose flow all exits one cell
  collapses to a single target. (If `F` touches the domain edge, an exit may go
  out of domain — ordinary off-domain case.)
- **piped** (each outflowing layer below the top `n`): explicit outlet cell,
  culvert-style delivery. Required — no default for lower layers.
- **bottom-layer infiltration** → lost to ground.

So `watercourses` provides only the static captured inflow (§1); the dynamic
overflow is delivered in the solver, exactly like a culvert outlet today.

## 4. Dynamic network inclusion

Causality moves with the flow. The NBS's outputs are **storage-driven**
(`compute_outflow` on layer storage, §3) with no dependence on downstream
level — so in the routine case coupling is one-way *both* upstream and
downstream:

- **Downstream of exit + piped outlets — included, one-way forcing.** The
  re-emitted overflow drives time-varying inflow into every trap it reaches, so
  those traps' fill timing depends on NBS output and they must be in the dynamic
  solve. The relationship is a forcing (NBS → downstream); downstream level does
  *not* feed back into the NBS rate (the outflow is storage-driven, not
  head-driven — a submerged outlet keeps discharging into the pool, it does not
  stop). Trace the affected set from `F`'s exit cells + piped outlet cells with
  `flow_path_from` (`watercourses.jl:142`), exactly as culvert outlets are
  included today. Saturation-dependent, so recompute on network touch, not once.

- **Upstream feeders — one-way input.** A feeder trap fills and spills on its own
  inflow; the NBS capturing that spill downstream does not change the feeder. The
  feeder is solved upstream-first and delivers its spill to the NBS as a one-way
  forcing input — ordinary topological-order routing. That `F` is not a region
  only means the crossing flow is diverted into the NBS at `F`; it does not
  couple the feeder.

- **The one genuine back-coupling: footprint backwater flooding.** If a
  downstream pool grows until it floods the *footprint itself*, the surface layer
  merges with the flood (§5) — downstream now affects the NBS. This is the
  existing supertrap/merge (backwater) mechanism, not NBS-specific, and not a
  routine outlet effect.

Needs a read-only footprint→cell map; `spillanalysis` stays untouched.

## 5. Submergence — surface phenomenon, binary threshold, NBS-as-drain

The footprint sits in whatever real trap the geometry gives it. That trap floods
via the ordinary machinery; the NBS is not in the hierarchy, so there is no state
handoff. Three facts shape the treatment:

**(a) Submergence is a *surface-layer* phenomenon.** Subsurface layers (soil,
drainage) are below grade — whether the water above is NBS-captured ponding or a
flood pool doesn't change how they infiltrate/percolate/drain (governed by their
own storage/saturation). So partial submergence reduces to "what happens to the
*surface* layer."

**(b) A submerged NBS is a drain, not frozen.** When flooded, the subsurface keeps
infiltrating, now fed by the flood → it **draws water out of the containing
trap** (green infrastructure at the bottom of a flood keeps soaking it up). The
piped outlet likewise keeps draining (storage-driven, independent of tailwater).
So a submerged NBS acts as a distributed infiltration + pipe drain on the trap it
sits in.

**(c) Lumped ⇒ no honest partial state ⇒ binary threshold.** A lumped store has no
spatial profile, so don't invent a fractional storage state. Threshold at
**exit-boundary submergence** (lowest footprint cells flooded) — exactly when
top-layer terrain re-emission stops working. At that event:
- surface layer merges into the containing trap's flood (stored water added to
  `cur_amounts`, conserved);
- subsurface layers keep running as the drain in (b);
- outlets keep discharging regardless of whether their cell is submerged — the
  outflow is storage-driven, not head-driven, so tailwater does not stop it;
- on recession below the threshold, surface layer re-activates.

**Same-trap outlet recirculation (explicitly handled).** When the footprint is
submerged, subsurface infiltration draws water *out* of the containing trap; if
an outlet's cell lies in that same flooded trap, the bottom-layer overflow flows
back *into* it. Net trap change = −(subsurface storage growth) −(ground-loss
infiltration); the recirculated overflow cancels exactly. No algebraic loop: the
outlet flow is storage-driven (a function of the NBS bottom-layer storage), the
abstraction is a function of trap level — so the ODE stays well-posed. No
special-casing is needed beyond delivering each outlet's flow to whatever trap
contains its cell (here, the same trap), which the position-exact delivery
already does. Pinned by a dedicated test (§6a).

Real per-cell footprint elevations are available (no carving), so the threshold
is well-defined and scheduled as a transition *event* (not per-step drift —
avoids the stale-state bug class).

**Upgrade path (deferred).** If the binary jump proves too crude for a footprint
spanning a large elevation range within one trap, scale only the *surface* layer
by dry-fraction `(1−f)`, `f` = submerged-cell fraction from the real elevations;
subsurface unaffected. Smooth but more bookkeeping/tests. Not in the first cut.

## 6. Mass ledger

```
inflow_external + direct_precip
   = Σ_layers dS/dt·A/1000 + top_overflow_reemit + bottom_overflow_piped
     + bottom_infil_to_ground + evap
```
Capture removes `inflow_external` from the downstream runoff exactly once;
re-emission adds the overflow back exactly once. Net terrain effect = the NBS
retained/infiltrated the difference. Same conservation shape as a full trap's
capture-and-spill.

## 6a. Mass-conservation testing (mandatory — bug-prone seam)

Capture-and-re-emit is the highest-risk path; tests must pin every route.
Minimum set:
- single-layer puddle, no submergence: `inflow = ΔS + overflow`, residual ≈ 0;
- multi-layer cascade: `total_in = ΔΣS + Σ overflow + ground_infil + evap`;
- **capture correctness:** downstream trap inflow drops by exactly
  `inflow_external − reemitted_overflow`;
- **non-convex footprint** (exit-and-re-enter): external inflow counted once;
- **submergence transition, both directions:** no water created/lost at the
  threshold event when the surface store folds into the trap and back;
- **submerged NBS as drain:** trap volume drops by exactly the infiltration drawn
  while submerged;
- **same-trap outlet recirculation:** submerged NBS whose outlet cell is in the
  same flooded trap — abstracted water returns to that trap, net change =
  subsurface storage growth + ground loss (recirculated part cancels), residual
  ≈ 0;
- **distributed re-emit:** weighted lower-edge re-emission conserves mass
  (Σ weights = 1, Σ delivered = overflow) and matches the natural (no-NBS) exit
  split;
- **NBS → NBS (variable-rate coupling):** upstream NBS re-emit reaching a
  downstream NBS's footprint is captured dynamically within the solve (§7a), and
  the chain conserves mass under a continuously-varying upstream rate;
- NBS-into-NBS ordering (see resolved decision 3);
- windowed-solve vs full-run residual (mirror existing culvert conservation
  tests).

## 7. Data model

```julia
struct NBSPlacement <: DynObject
    system::NBSSystem
    footprint::Vector{Int}                # linear indices
    n_terrain::Int                        # the n topmost layers re-emit at F's exit boundary
    outlets::Vector{CartesianIndex{2}}    # explicit piped outlet, one per outflowing
                                          # layer below the top n (required)
end
```
No spillpoint, no region, no `-4` flag, no carve. Surface layer is a **distinct
store** (decision 2), following the calibrated NBS model when dry and merging into
the containing trap on submergence (§5). Construction validates that every
outflowing layer (`Kout>0`) below the top `n` has an explicit outlet.

The terrain layers' re-emit targets are **not** a struct field: the weighted
lower-edge exit-cell list is *derived at setup* from the footprint geometry (§3),
stored in the NBS routing plan (`deliver_slot` gains a weighted slot-list per
terrain layer; piped layers keep a single slot).

## 7a. Dynamic representation: NBS is a network *element*, not a trap node

Orientation finding (reading `networksolver.jl` + `network_context.jl`): the
current PR architects the NBS as a **trap node** in `net.traps` — its top-layer
inflow is `inflow[plan.trap_local]`, the flow routed to that node, which exists
only because the dug trap forms a region whose catchment feeds `trap_inflow`;
`_nbs_elements`, `_dyn_seeds`, `_external_inflow`, `_make_context` all key off
`t.trap_ix` / `tstruct.footprints[nb.trap_ix]`. Removing the dug trap removes
this inflow path, so the NBS must be re-architected into a standalone network
**element** (it is already `DynNBS <: DynObject`), like a culvert: external
captured inflow + outlet delivery, not a trap.

**Inflow has two parts — the footprint is a flow terminus in two places:**

1. **Static terrain inflow** (per-event constant): the §1 sink-overlay tally in
   `watercourses`, fed as external inflow to the layer ODE. Piecewise-constant
   because ordinary traps only change output at `:fill`/`:empty` events, which
   terminate the solve and refresh rateinfo. This is the footprint as a *static*
   terminus (flow terminates there in the accumulation sweep).
2. **Dynamic internal inflow**: any flow *internal to the network* that reaches
   the footprint — most importantly an **upstream NBS's re-emit** (whose rate
   varies continuously, storage-driven, with no bounding event), and culvert
   deliveries. This must be captured **within the solve** (`_route_flow` runs
   every ODE-step). This is the footprint as a *dynamic* terminus: flow paths
   entering the footprint must terminate at the NBS element, so path flow is
   drawn in instead of continuing downstream.

So the NBS element is *both* fed by a static external constant *and* a dynamic
routing target — not a pure external-inflow sink. (Earlier drafts wrongly said
"external constant, not a routing sink"; the NBS→NBS variable-rate case proves it
must be both.) The current dug-trap design gets the dynamic terminus for free (the
dug region is a path target); the overlay must make the footprint a path terminus
explicitly at network-build time (`setup_network`).

**Reusable as-is:** layer ODE cascade (`compute_outflow`, mm↔m³), outlet-delivery
slots (`nbs_actual` / `:nbsout` / `nbs_trap_outlets`), placement-keyed persistence
(`nbs_state`), and the routing-node nature (delivery to footprint cells → the NBS).
**Must change:** the *terrain* inflow source (routed-dug-region-inflow → sink-tally
`nbs_inflow` from `watercourses`); `_nbs_elements` / `_dyn_seeds` / `_make_context`
/ `_external_inflow` to key off the real footprint (not a dug trap); and
`setup_network` to make the real footprint a path terminus.

## 8. Delta vs current PR (nbs-dynamic-core)

**Keep (reusable as-is):** `NBSLayer` / `NBSSystem` / `compute_outflow`,
factories, mm↔m³ conversion, the layer ODE cascade, the outlet-delivery slots
(`nbs_actual` / `:nbsout` / `nbs_trap_outlets`), placement-keyed persistence
(`nbs_state`). Keep the footprint→layer-area computation (currently in
`_prepare_nbs!`): area `A` still feeds the mm↔m³ conversion — move it into the
`fill_sequence` setup.

**Remove:** the static terrain modelling — `_dig_nbs_traps!` (dig/carve),
`_resolve_nbs!`'s outlet-on-spillpoint, the `-4`/region concept.
`spillanalysis` loses its `nbs` kwarg entirely (NBS becomes purely a
`fill_sequence` concern).

**Re-architect (the crux, §7a):** the DynNetwork NBS from a *dug-trap node* to an
*overlay routing element* — rewire `_nbs_elements` (no region lookup),
`_dyn_seeds`, `_make_context`, `_external_inflow`/`_inflow_sources`, and
`_routed_inflow`'s top-layer feed (static `nbs_inflow` + dynamic footprint
capture). Make the real footprint a **path terminus** in `setup_network` so
dynamic internal flow (upstream-NBS re-emit, culvert delivery) reaching it is
captured within the solve — the dug region did this implicitly.

**Add:** in `watercourses` — the static footprint-as-sink overlay (capture +
per-placement `nbs_inflow` tally, §1/§2), plus that field on `rateinfo`. In the
solver — feed `nbs_inflow` as the element's external inflow; deliver dynamic
overflow to resolved outlet cells (terrain exit-boundary / piped). Plus
footprint→containing-trap + exit-boundary network inclusion, and the submergence
event (drain + recirculation).

## 8a. Implementation staging (each a reviewable commit)

- **A** — strip static NBS (dig / `_resolve_nbs!` / `-4` / `spillanalysis`
  kwarg); keep footprint→area.
- **B1** — footprint-as-sink overlay in `watercourses` (§1 tally + §2 capture) +
  the `nbs_inflow` field on `rateinfo` + mass tests for capture/no-double-count.
- **B2** — re-architect the DynNetwork NBS trap-node → overlay element fed by
  B1's `nbs_inflow` (static) + dynamic footprint capture; rewire `_nbs_elements`
  / seeds / context / `_routed_inflow`; make the footprint a path terminus in
  `setup_network` (dynamic internal capture, incl. NBS→NBS); precompute the
  weighted lower-edge exit-cell lists; extend `deliver_slot` to a weighted
  slot-list per terrain layer; keep layer ODE + persistence; dynamic distributed
  re-emit (§3) delivered here. + mass tests (dry cases, distributed re-emit,
  NBS→NBS variable-rate).
- **C** — submergence event + drain + recirculation + tests (§6a, submerged
  cases).

## Resolved decisions

1. **n-topmost-terrain routing.** The `n` topmost layers re-emit at `F`'s exit
   boundary; every outflowing layer below them requires an explicit piped outlet
   (§3). No arbitrary per-layer mix — terrain-drain is a contiguous top block of
   size `n`.
2. **Surface layer is a distinct store.** Follows the calibrated NBS model when
   dry; merges into the containing trap on submergence (§5). Preserves calibration
   fidelity; submergence is surface-only (§5(a)) so the merge is contained.
3. **NBS-into-NBS resolved by topological ordering.** When one NBS's re-emission
   lands on another's footprint, order the NBS elements topologically (feeder
   before receiver), like culvert inlet-before-outlet, so capture-then-re-emit is
   mass-conservative and order-independent. Only fall back to "unsupported + TODO"
   if a genuine cycle (mutual NBS-into-NBS) appears — deferred like reverse
   culverts.

## Deferred

4. **Evapotranspiration.** The `EVCoeff`/`EVS11` layer params are not applied in
   the solver yet. Keep the per-layer rate function structured so an ET sink is a
   clean addition — each layer's balance is written as
   `dS/dt = inflow − infiltration_out − overflow_out − ET`, with the `ET` term
   present as an explicit `0.0` placeholder computed from `EVCoeff`/`EVS11` (not
   omitted), so wiring it in later is a one-line change per layer, not a
   restructuring. Do **not** bake in an assumption that ET is absent (e.g. don't
   drop the term from the mass ledger).
