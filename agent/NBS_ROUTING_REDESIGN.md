# NBS routing redesign — design notes

Status: **implemented**. `watercourses` is NBS-oblivious; the correction lives in the
network solver (`networksolver.jl`: `NBSElement`/`NBSPlan`, `_build_nbs_plan`,
`_propagate_correction`, `_apply_nbs_corrections!`).  The old footprint-as-sink +
positive re-emit path (sink overlay, `nbs_actual` slots, `_nbs_exit_weights`, `nbs_into`,
submergence, `RouteScratch`) is deleted.  Deferred: submergence and NBS→NBS series (§11).

Scope: how Nature-Based Solution (NBS) installations couple to the dynamic
surface-flow network. Supersedes the overlay-element approach in
`NBS_OPTION1_OVERLAY_PLAN.md` for the routing/coupling layer. The NBS *layer
model* primitives — `NBSSystem` / `NBSLayer` / `compute_outflow` in
`nbs_elements.jl` — are retained unchanged.

Network build entry point: `src/dynamics/build_network.jl` (`setup_network`).
The correction consumer lives in the solver / rateinfo layer
(`networksolver.jl`, `network_context.jl`).

---

## 1. Core idea

An NBS is **invisible to everything outside `dynamics/`**. `watercourses` /
`_compute_initial_rateinfo` route the fixed overland flow across the NBS
footprint exactly as if no NBS were present — the baseline that sets regular
trap fill times is NBS-oblivious. The NBS then acts *only* as a **correction**
to that baseline, applied entirely inside the dynamics layer:

- it captures water on its footprint,
- stores it in its layered `NBSSystem` ODE,
- releases it — as terrain re-emit at the footprint's exit boundary and as
  piped discharge at explicit outlets — from where ordinary flow paths carry it
  downstream.

The correction is a **signed flow** propagated along the existing network
routing: negative where the NBS is currently holding water back, positive where
its store is draining downstream faster than water is currently arriving. It is
never a rewrite of the surface-network topology.

This mirrors how the dynamic network already corrects the fixed flow elsewhere:
a static baseline, plus a dynamic correction where the network is active.

## 2. Water balance

For one NBS at any instant, the layer ODE keeps its internal balance:

```
capture         I         =  precip on footprint  +  inflow across the footprint boundary
storage         dS        =  d/dt of total layer storage
terrain re-emit O_terrain  =  overflow of the top n_terrain layers, re-emitted at the
                              footprint's natural drainage targets (see below)
piped discharge E_piped    =  discharge of the outflowing layers below the top n_terrain

I = dS + O_terrain + E_piped
```

Two **independent** surface effects, applied at different locations — they are
NOT lumped into a single emission:

- **Terrain re-emit** — a signed correction applied at the footprint's
  **natural drainage endpoints**: its **boundary exits** *and* its **internal
  depressions** (footprint trap-bottom cells). Each endpoint `c` carries an
  oblivious NBS flow `Q_c` (what the capture drained there in the no-NBS run;
  `I = Σ_c Q_c`), and the correction there is simply

  ```
  c_in(c) = −X · Q_c,   with the single scalar   X = 1 − O_terrain/I   (§4.3)
  ```

  No weight vector — `Q_c` *is* the split. Boundary corrections propagate
  downstream along the seeded path; internal ones deduct from the **containing
  trap** (§4.2). Summed over endpoints the total is `−X·I = O_terrain − I`, so
  everything the top layer stored or passed to the piped layers is deducted here
  — the only place the balance shows up.
- **Outlets** — the piped discharge `E_piped`, injected **as-is** at the outlet
  cells. `E_piped` is a direct output of the ODE (computed from the piped
  layers' storage); it is *not* a fraction of anything and *not* affected by the
  terrain balance. The water it carries was already deducted at the terrain
  targets (inside the `−I`), so re-adding it at the pipe is correct double-entry,
  not double-count.

Net water removed from the surface = `Σ_c (−X·Q_c) + E_piped =
(O_terrain − I) + E_piped = −dS`: exactly what the NBS is holding (negative while
filling, positive while draining). At steady state `dS = 0` and both effects
cancel the baseline back to correct. Mass closes at all times because
`O_terrain`, `E_piped`, `I` all come from the ODE.

The submergence coupling (the containing trap flooding the NBS when its water
level rises above the footprint's lowest cell) is a **separate**, geometric
interaction, not part of this balance.

**`P = I` via the zero-infiltration contract.** The correction is applied
against the oblivious pass-through `P` — the capture disposed of across **all**
the footprint's drainage endpoints (boundary exits *and* internal depressions),
`P = Σ_c Q_c`. For `P` to equal the capture `I`, the footprint must carry **zero
terrain infiltration** — else the oblivious run infiltrates under the footprint
*and* the layer model infiltrates, double-counting. This is the standing
contract (currently a `# @@@` note in `build_network.jl`); it must be
**enforced** at `setup_network`, not merely assumed. With it, oblivious `P = I`
and the correction is well defined.

## 3. The four cell-sets (computed in `build_network.jl`)

`_compute_nbs_inflow_outflow_cells!` classifies each footprint's coupling cells:

| Set | Definition | Role |
|-----|-----------|------|
| `footprint_inflow_cells`  | cells feeding *into* the footprint from outside (inverse-flow neighbours of footprint cells, minus the footprint) | where capture `I` is sourced |
| `footprint_outflow_cells` | footprint cells whose downstream neighbour lies *outside* the footprint | where the terrain re-emit correction is injected (at the external `ds` cell, which is seeded) |
| `internal_accumulation_cells` | footprint ∩ trap-bottoms | internal terrain re-emit target — its correction share lands on the containing trap |
| `outlets` (user-supplied) | piped-discharge cells for the outflowing layers below the top `n_terrain` | where piped emission re-enters the terrain (seeded) |

## 4. The signed correction — propagation

The NBS internals are not spatially modelled. The correction is a single signed
scalar per exit, injected at that exit cell and propagated downstream through
the **existing network routing** (paths, traps, spill paths) — not a bespoke
walk. This replaces the earlier equal-split-plus-clamp heuristic: the routing
below *is* the clamp, applied exactly per cell along the real flow route.

### 4.1 Per-cell propagation rule

`watercourses` already stores, per cell, the signed value we need: `runoff[c]`
starts at `−infiltration_capacity` and accumulates upstream flow, so the final
grid holds **positive = runoff actually emitted**, **negative = remaining
infiltration capacity**. Call it `V_c`; the cell's emitted runoff is
`max(V_c, 0)`.

A signed correction `c_in` arriving at cell `c` propagates to `c`'s downstream
neighbour as:

```
c_out = max(V_c + c_in, 0) − max(V_c, 0)
```

(and, if `V_c` is reused downstream, update `V_c ← max(V_c + c_in, −infil_c)`).

One expression, both signs:

- **Retaining** (`c_in < 0`): while `|c_in| ≤ V_c` the full correction passes
  (`c_out = c_in`); once it would push the cell's runoff below zero it
  attenuates to `−max(V_c, 0)` (all the runoff there was) and dies — you cannot
  remove water that was already infiltrating.
- **Draining** (`c_in > 0`): the extra first refills any spare capacity
  (`V_c < 0`); only the surplus above zero continues as new runoff. So a release
  onto a segment that had run dry in the oblivious field is re-absorbed before
  it reaches a trap.

A naïve `c_out = c_in` (correction passes unchanged) is wrong: it over-subtracts
where the cell's runoff can't absorb the full correction, and over-delivers a
release where the terrain has spare capacity. The `max`/`max` form is the exact
attenuation.

### 4.2 Where the correction enters

- **Terrain re-emit** (top `n_terrain` layers): one scalar
  `X = (dS + E_piped)/I` per NBS per step (§4.3), applied cell-by-cell as
  `c_in(c) = −X · Q_c` at each
  drainage endpoint `c`. `Q_c` is the oblivious NBS flow there, read straight
  from the field (the runoff at the footprint-side boundary cell — pure NBS,
  since the footprint holds the capture; `I = Σ_c Q_c`). No weight vector — `Q_c`
  *is* the split. A **boundary** endpoint's `c_in` is injected at its external
  `ds` cell and propagated along the seeded path by §4.1; an
  **internal-depression** endpoint's `c_in` is deducted from the **containing
  trap** (the `DynTrap` whose footprint contains the cell — resolved via the net,
  not by re-deriving `regions[cell]`; build a cell→local-trap map once per solve,
  this runs in the ODE hot path). This subsumes the old separate "internal
  accumulation deduction" — it is just an internal drainage endpoint.

  `−X·Q_c` needs no precomputed exit weights at all: it never over-deducts (a
  fraction of the flow actually present, so no clamp/redistribute) and reflects
  real inflow concentration. `build_network` still supplies *which* cells are the
  endpoints (`footprint_outflow_cells` / `internal_accumulation_cells`), just not
  their weights.

  **`I ≈ 0` fallback:** during pure drainage (dry period, store emptying) every
  `Q_c = 0` and `X` is singular, so `−X·Q_c` is indeterminate yet
  `O_terrain > 0` must still be placed. Split it **equally** among all the
  endpoint cells — boundary exits and internal depressions alike. No weights, no
  graph traversal.
- **Piped outlets** (outflowing layers below the top `n_terrain`): `E_piped`,
  a direct positive injection at each outlet cell — the ODE's piped-layer
  discharge, injected as-is (§2), still propagated downstream by the §4.1 rule
  (it just always starts positive). Not tied to the terrain balance.

Terrain correction + piped injection = `(O_terrain − I) + E_piped = −dS`. Mass
conserved (§2).

### 4.3 `X` is used unclamped

`X` is a single scalar per NBS per step, used directly as the multiplier
`c_in(c) = −X·Q_c`:

```
X = 1 − O_terrain/I
```

the fraction of capture *not* re-emitted over the terrain. `O_terrain` is the
top-`n_terrain` layers' overflow, read directly from the layer state. **Use this
form**, not `(dS + E_piped)/I` — the two are equal only when nothing infiltrates
to ground. The layer model's bottom layer *does* infiltrate to ground (a fourth
sink `G`), so the real balance is `I = dS + O_terrain + E_piped + G` and
`(dS + E_piped)/I` understates retention by `G/I`, leaving ground-lost water
wrongly flowing downstream (mass not conserved). `1 − O_terrain/I` is correct
regardless of where the non-re-emitted capture went. Do **not** clamp `X` to
`[0,1]` — it legitimately leaves that range:

- `X ≤ 1` always holds (`O_terrain ≥ 0`).
- `X ≥ 0` **fails**: a store filled under heavy rain keeps overflowing after
  inflow drops (`compute_outflow` depends on storage `S`, not on `I`), so
  `O_terrain > I` and `X < 0` — a **positive** correction, the store releasing
  downstream, which is physical. Concretely: reach steady state at `P = 100`
  (storage high), then rain eases to `P = 10`; the store still overflows `≈ 100`,
  so `X = (10 − 100)/10 = −9`. Fires on the first step of a lighter weather
  period because `nbs_state` carries storage across periods.

The product `−X·Q_c` stays **bounded** even at large-negative `X`, because
`Q_c ≤ I`, so `|X·Q_c| ≤ |I − O_terrain|`. The only breakdown is `I = 0` (pure
drainage): `X` is singular *and* every `Q_c = 0`, so `−X·Q_c` is `0·∞` — handled
by the §4.2 equal-split fallback, not by clamping.

Clamping `X` to `[0,1]` would zero the drainage release and destroy mass.
**A brief comment at the correction site must record this** so no one
"helpfully" re-adds a clamp. (The piped discharge `E_piped` needs none of this —
it is `≥ 0` and injected as-is.)

## 5. Woven into the dynamic router, not a one-shot pass

`c_in(t)` is dynamic (the store saturates and drains) and interacts with the
`V_c` caps nonlinearly, so it is evaluated inside the network's per-step routing,
not precomputed. Crucially, a large enough retention can drop a downstream
**full** trap's net inflow below its losses and flip it to `:unspill` — a
topology event. So the correction is part of the trap inflow the ODE sees each
step, routed through traps and their spill paths by the same topology
`_route_flow` already handles. The extra input the router needs is the oblivious
`V_c` per path cell — read straight off the `watercourses` `runoff` grid.

## 6. NBS in series

The oblivious field routes NBS-A's *full* pass-through into NBS-B's footprint, so
B's oblivious capture over-counts. A's correction must land on B's footprint
cells and lower B's actual `I`, which changes B's `O_terrain` and hence B's own
correction. This is causal, resolved in topological order — the same coupling
the old `nbs_into` slots encoded, now expressed as signed corrections landing on
downstream footprints.

## 7. Division of labour

**`watercourses` / `_compute_initial_rateinfo` (NBS-oblivious):**
- compute the plain flow field; expose the `runoff` grid (`V_c`). No `nbs`
  argument, no sink, no `nbs_inflow` return. This is a net deletion.

**`build_network.jl` (topology):**
- **enforce** zero infiltration on footprints (§2);
- identify the four cell-sets per NBS;
- seed outlets, outflow-boundary external `ds` cells, and internal-accumulation
  cells; grow the monolithic network; split into components.

**Solver / rate layer (flux):**
- aggregate capture `I` from the oblivious field (`runoff` at
  `footprint_inflow_cells` + footprint precip);
- run the layer ODE to get `O_terrain(t)` and `E_piped(t)` separately;
- inject the signed terrain correction and the piped discharge (§4.2) and
  propagate them through the router by
  the §4.1 rule, respecting downstream trap-state flips (§5) and NBS-series
  ordering (§6).

## 8. Component splitting

`split_network_into_connected_components` (in `network_utils.jl`) is unchanged:
an NBS bridges otherwise-disjoint regions and must land its coupled elements in
one component. Edges from terrain flow + merges, culvert inlet/outlet owners,
and per NBS: emission paths (`nbs_outlets`), outflow-deduction paths (seeded at
`footprint_outflow_cells`), and the accumulation trap (mapped via
`regions[c] → supertraps_of`, asserted unique).

## 9. Cycle safety

Piped outlets, like culverts, are not constrained to terrain flow direction; an
outlet upstream of the NBS's own inflow could loop and break the
"water flows downstream" invariant. Handled by **user contract**: outlets must
not sit above the footprint elevation. No automatic detection (matches the
deferred reverse-culvert handling).

## 10. What changes vs the current code

**Deleted:**
- the `nbs` argument, footprint-as-sink overlay, and `nbs_inflow` return in
  `watercourses`;
- the positive re-emit delivery-slot machinery in the solver (`nbs_actual`
  slots as a separate positive delivery path, `nbs_into` as positive slots,
  `nbs_path_events` / `nbs_trap_outlets` as the sole delivery mechanism);
- `_nbs_exit_weights` entirely — the split is `Q_c` read from the field, and the
  `I ≈ 0` fallback is a plain equal split over the endpoint cells (no weights, no
  traversal).

**Kept (shrunk):**
- the endpoint cell-sets (`footprint_outflow_cells` / `internal_accumulation_cells`)
  — now used to locate where to read `Q_c`;
- piped-outlet seeding and injection;
- the layer ODE (`NBSSystem` / `compute_outflow`) and cross-period `nbs_state`
  persistence.

**Added:**
- the oblivious `runoff` grid threaded to the solver (`rateinfo.runoff` →
  `DynNetworkContext.runoff` → `_build_nbs_plan`);
- `_build_nbs_plan` (endpoints via `_nbs_endpoints`/`_walk_to_trap`, outlets via
  `_nbs_outlets`), `_propagate_correction`, `_apply_nbs_corrections!`;
- zero-infiltration on footprints (enforced in `fill_sequence`).

## 11. Deferred / open

- **Submergence** — the containing trap flooding the footprint (surface block merges into
  the flood, terrain re-emit stops).  Removed with the old overlay; re-add cleanly as a
  separate geometric interaction if needed.
- **NBS→NBS series** — an upstream element's correction landing on a downstream footprint
  lowers its capture `I`.  Not modelled: each element's `I` is read from the oblivious
  field, which over-counts an upstream element's full pass-through.
- Per-layer distinct outlets (single shared outlet for now; `nbs_outlets`
  currently a 2-tuple `(placement_ix, position)`, extend to add `layer_ix`).
- Evapotranspiration (explicit `0.0` placeholder in the layer loop today).
- `watercourses`'s `used_infiltration` return is the *oblivious* total (computed
  once from the no-NBS flow field); the NBS correction changes actual downstream
  flow and hence true infiltration, so the two differ. Harmless today — the only
  caller discards it (`flow.jl` `_`). Reconcile only if some future caller reads
  it as an NBS-aware water budget.
- Far-downstream *static* traps keep their oblivious (slightly-high) changetime
  estimate until the frontier reaches them and they join the network — the same
  approximation the network already makes.
