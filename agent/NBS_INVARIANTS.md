# NBS module: terms, model, code map and invariants

This is the authoritative description of how Nature-Based Solutions (NBS) are coupled to SWIM. Code comments refer here rather than repeating it. Older design documents are in `agent/archive/` and are obsolete.

New to the NBS code? Start with `NBS_INTRO_OVERVIEW.md` (main logic and workflow), then `NBS_INTRO_DETAILED.md` (function-by-function walkthrough).

## Terms

Each term below is defined using only base SWIM vocabulary and terms defined above it.

**The placement**

- **NBS placement** (`DynNBSPlacement`): an NBS installed on the terrain. It consists of a footprint, a layer stack, `n_terrain` and outlets, each defined just below.
- **footprint**: the grid cells the placement covers. Always "NBS footprint" here; a trap's footprint is called *trap footprint*.
- **layer stack** (`NBSSystem`): storage layers, top first. Water enters the top layer.
  - **storage** `S`: the water a layer holds, in mm (`V` when expressed as a volume).
  - Each layer has two thresholds: `Smax` controls its overflow and `Smin` its percolation.
- **overflow**: what a layer sheds above its `Smax`: `Kout·(S − Smax)^nout` (coefficient `Kout`, exponent `nout`), zero at or below `Smax`.
- **percolation**: what a layer passes to the layer below above its `Smin`: `Kinf·(S − Smin)^ninf` (coefficient `Kinf`, exponent `ninf`), zero at or below `Smin`. The bottom layer's percolation leaves the model ("to ground").
- **`n_terrain`**: the number of top layers whose overflow goes back onto the terrain.
- **piped layer**: a layer below the top `n_terrain` that has overflow (`Kout > 0`). The list is `_piped_layers(system, n_terrain)`, top first.
- **outlet**: the cell where a piped layer's overflow (its pipe flow) re-enters the terrain. `outlets[k]` belongs to the k-th piped layer.
- **`O_t`**: the total overflow of the top `n_terrain` layers, i.e. what the placement returns onto the terrain. **`E_k`**: the overflow (pipe flow) of the k-th piped layer.
- **`id`**: the placement's position in the `nbs` vector passed to `fill_sequence`. It labels the placement's stored layer volumes.

**The terrain and the baseline**

- **baseline** (*oblivious* in code): the flow field and trap inflows base SWIM computes. It treats footprints as ordinary terrain.
- **runoff value** `R`: the baseline runoff grid holds, per cell, either positive surface flow or, if negative, unused infiltration capacity.
- **inflow cells**: cells outside a footprint whose downstream cell (in the flow graph) is inside it.
- **exits**: cells outside a footprint that a footprint cell drains to.
- **ponds**: trap bottoms inside a footprint.
- **leaf regions** of a trap: the lowest-level spill regions it contains (just its own spill region, if it has no sub-traps).
- **`O_0`**: the footprint's baseline throughput: rain on the footprint plus the baseline flow entering it through its inflow cells.
- **share** (`share_e`): the fraction of `O_0` that leaves the footprint through endpoint `e`, where an endpoint is an exit or a pond.
  - **For an exit:** the baseline flow of the footprint cells that drain straight out to that exit, divided by `O_0`. The exit cell's own runoff is not used, because it also includes water that never crossed the footprint.
  - **For a pond:** the baseline flow arriving at the pond, divided by `O_0`.
  - The shares sum to 1. If there is (almost) no throughput, they are set equal.

**The dynamic network**

- **culvert** (`DynCulvert`): a pipe moving water from an inlet cell to an outlet cell, at a rate set by the water levels at both ends. Like an NBS, it needs the dynamic network because its rate changes continuously with the water levels.
- **seed**: a cell from which flow is traced. Seeds are each NBS's exits, outlets and ponds, each culvert's two ends, and one cell (the first footprint cell) of each trap named in `fill_sequence`'s `dyn_traps` argument; these cells are passed internally as `dyn_coords`.
- **network path** (`DynFlowPath`): a run of cells traced downhill from a seed, or from a full trap's spillpoint, to the next trap (or out of the domain).
  - **departure point:** the cell the tracing started from (the seed or the spillpoint). The path's cells begin with it, except that cells inside the target trap's footprint are left out, and a path whose first cell was already claimed by an earlier path gets no cells and just connects its departure point to that path. So a pond seed's path, or one starting inside a trap footprint, is empty.
  - **head:** the start of the path, where flow added to the path (e.g. NBS output) enters.
  - A path can end where it joins another path; it is then a **tributary** of that path.
- **network trap**: a trap reached by such tracing, where tracing continues through full traps along their spill paths. Its volume is an ODE state. The set changes during a run:
  - it grows when a network trap fills and its spill path is traced further;
  - traps drop out when they are no longer reached.
- **covered trap**: a network trap, or a full sub-trap inside one. Base SWIM leaves covered traps alone: they are not in its spill graph, and their changetime estimates are replaced by the network's prediction (an exact time for the one trap each component predicts to change next, `Inf` for the others).
- **component** (a `DynNetwork`): a set of network paths, network traps, culverts and NBS placements that are linked (by flow, culverts or NBS) and so are solved as one ODE. One run has several independent components. When a trap stops spilling (`apply_unfill!`), a component that becomes disconnected is not split: it stays one ODE, which gives the same result. When a component is regrown (a merge, or a parent trap absorbing or releasing its sub-traps), it is split again into connected pieces.
- **carrier path**: the network path departing from an exit or an outlet. NBS output enters at its head.
- **routing pass**: one evaluation of "what arrives where" in a component (`_route_flow`), done at every evaluation of the ODE right-hand side. It visits paths and traps from upstream to downstream.
- **network flow**: what the routing pass carries along network paths. It is always an *increment on top of the baseline*. Examples are the spill of a full network trap (base SWIM's spill graph excludes covered traps, so this spill exists only here), a culvert's flow, and NBS corrections. It is signed.
- **per-cell rule**: how network flow `d` changes along a path. At each cell with baseline runoff value `R`, `d ← max(R+d, 0) − max(R, 0)`. Unused infiltration capacity absorbs added water, and a negative `d` can remove at most the flow that is there. For added water this is the rule base SWIM's `_track_flow!` uses when it routes a trap's spill.
- **draw**: the network flow leaving a placement's inflow cells into the footprint (after each inflow cell's own per-cell rule) during one routing pass. It is signed: negative when an upstream NBS is holding water back.
- **injection**: what a placement adds to a routing pass: `(O_t − O_0)·share_e` at the head of each exit's carrier path or into each pond's trap, and `E_k` at the head of each outlet's carrier path.

**Running the network**

- **driver**: the object that owns all components of one weather period, their committed volumes and `nbs_state`.
- **`nbs_state`**: the stored layer volumes of every placement, by `id`. It is carried across touches and weather periods.
- **context** (`DynNetworkContext`): one component together with its *committed* state (trap volumes and layer volumes at a known time) and its rates (`params`).
- **touch**: what happens when an event in `fill_sequence` affects any component (one of its traps changed status, or the baseline inflow to one of its network traps' leaf regions changed). For **all** components of the driver, not only the affected one:
  - **commit:** integrate each context's ODE up to the event time;
  - **harvest:** copy the committed trap volumes into the driver (`vol_by_trapix`) and the layer volumes into `nbs_state`;
  - apply the change to the affected component (it may grow, lose traps, or merge with another);
  - **reconcile** base SWIM's spill graph with the new coverage (see S6);
  - **rebuild:** build a new context for every component, with rates from the current baseline;
  - **predict:** integrate a copy of each ahead to find its next event.
- **frozen rates**: a context's `params` are built at a touch and not re-read from the live baseline grids until the next touch.
- **finalize**: commit every context to the end of a weather period.

## Model

**Base SWIM stays NBS-oblivious.** `watercourses`, `spillanalysis` and the constant-rate `fill_sequence` machinery never see an NBS.

**`external_inflow`**: each network trap's baseline inflow. It is base SWIM's `rateinfo.trap_inflow`, summed over the trap's leaf regions (`_external_inflow`). It includes, after infiltration on the way, the throughput `O_0` of any footprint upstream, because the baseline routes water straight through footprints. The routing pass starts every network trap from `external_inflow` and adds network flow on top.

**Per placement, at every evaluation of the ODE right-hand side:**

```
Q_in  = O_0 + draw  (≥ 0: per inflow cell c it is max(R_c + d_c, 0))   input to the top layer
dV_l  = in_l − overflow_l − percolation_l            in_1 = Q_in,  in_l = percolation_{l−1}
O_t   = overflow of the top n_terrain layers         (back onto the terrain)
E_k   = overflow of the k-th piped layer             (pipe flow, to outlets[k])
```

**Injections into the routing pass:**
- **Exits:** `(O_t − O_0)·share_e` at the head of each exit's carrier path.
- **Ponds:** `(O_t − O_0)·share_e` into each pond's trap.
- **Outlets:** `E_k` at the head of each outlet's carrier path.

The baseline `O_0·share_e` reaches the trap below exit `e` reduced by the infiltration along the way, and is already in that trap's `external_inflow`. The injection travels the same path under the per-cell rule. Together, the trap receives `O_t·share_e` reduced by the same path, in place of `O_0·share_e` reduced by it: what the NBS actually releases.

## Code map

| Concern | File | Key functions |
|---|---|---|
| Layer model (mm, no terrain) | `dynamics/nbs_elements.jl` | `NBSLayer`, `NBSSystem`, model constructors, `compute_outflow`, `_layer_fluxes`, `_nbs_layer_rates!`, `_nbs_terrain_outflow` |
| Placement (input) | `dynamics/elements.jl` | `DynNBSPlacement`, `_piped_layers`, `NBSLink` |
| Placement ↔ terrain | `dynamics/nbs_terrain.jl` | `_prepare_nbs!`, `_nbs_link`, `_nbs_baseline`, `_pond_trap`, `_nbs_zero_infiltration!` |
| Network membership | `dynamics/build_network.jl`, `dynamics/network_utils.jl` | seeds in `setup_network`; `_intersecting_on_path`; `_nbs_coupled_nodes` |
| Solver | `dynamics/networksolver.jl` | `_build_nbs_plan` → `NBSPlan`/`NBSElement`; `_nbs_injections`; `_route_flow`; `dynNetworkRateFunction!` |
| Persistence | `dynamics/network_context.jl` | `_nbs_layer_block` (restore), `_store_nbs_state!` / `_nbs_layers` (save) |
| Flow field at a timepoint | `dynamics/network_watercourses.jl` | `network_watercourses` |

## Invariants

Each invariant names where it is enforced.

**Inputs**
- **I1.** Footprints are non-empty and do not overlap. *`_validate_network_inputs` (error).*
- **I2.** No culvert end, `dyn_traps` cell or sink (a cell where SWIM removes water) lies inside a footprint. A sink would end flow at a cell the placement never sees. *`_validate_network_inputs` (error).*
- **I3.** Infiltration is zero on every footprint cell, because the layers own that water: the footprint is impermeable to base SWIM, and any infiltration there happens inside the NBS element. Everything computed from that grid (the runoff field, and from it the network traps' ponding losses, S12) therefore also treats the footprint as impermeable. *`_nbs_zero_infiltration!`, applied wherever an infiltration grid is built.*
- **I4.** `length(outlets)` equals the number of piped layers. *`DynNBSPlacement` constructor (error).*
- **I5.** `id` is stable for a whole `fill_sequence` call. *Set by `_prepare_nbs!`.*

**Terrain link and baseline**
- **L1.** `O_0 = rain on footprint + Σ_{c ∈ inflow cells} max(R_c, 0)`. In the baseline the footprint has no infiltration (I3), so this is exactly what leaves it through its exits and ponds.
- **L2.** The shares are ≥ 0 and sum to 1. *Asserted in `_nbs_baseline`: the unnormalised weights must sum to `O_0`.*
- **L3.** Every exit and outlet is the departure point of a network path (so it has a carrier path), and every pond lies in exactly one network trap. *Asserted in `_build_nbs_plan` and `_pond_trap`.*

**Network membership**
- **C1.** A placement is in the same component as the carrier paths of its exits and outlets, every path crossing its inflow cells, and its ponds' traps. This is what makes its draw complete within one routing pass. *`_nbs_coupled_nodes`.*

**Solver**
- **S1.** The ODE state is laid out as the component's `nt` trap volumes, then each placement's layer volumes in the component's placement order (`net.nbs`), top layer first. *`_build_nbs_plan` records each placement's offset (`state_base`); `_layer_range` reads it; `_nbs_layer_block` builds the state in the same order.*
- **S2.** Injections depend only on the layer volumes, and draws feed only the layers' rates. So no placement needs to be routed before another, and there is no circular dependency within a routing pass. *Call order in `_routed_inflow`.*
- **S3.** Per placement, `Σ dV = Q_in − O_t − Σ E_k − bottom percolation`. *Unit test "NBS: layer cascade conserves mass".*
- **S4.** A layer that is empty, or at or below a threshold, emits nothing through that flux. No layer rate is clamped ("floored") after the fact; the only `max(·, 0)` is on `Q_in`, which is ≥ 0 by construction and only absorbs rounding. *`compute_outflow`, `_layer_fluxes`, `dynNetworkRateFunction!`.*
- **S5.** A commit uses the context's frozen rates. *`DynNetworkContext.params`; fixed in `d0e296a` (review item C14); test "a network commit is not altered by the event that triggers it".*
- **S6.** Base SWIM's spill graph contains exactly the full traps that are not covered, so a covered trap's spill is routed only by its network, never twice. *`_reconcile_spillgraph!` at every touch that changes coverage, and asserted in `fill_sequence` each time a weather period's network is built. In the first period, a covered trap could still be in the spill graph only if it was already full when the network was built; that means zero volume, and `spillregions` merges zero-volume traps into their downstream region. In later periods the spill graph comes from the previous period's touches while coverage comes from a fresh build, so the assertion also checks that the two agree.*
- **S7.** In the solver a trap spills if and only if it is full in its network (it has a spill path), whatever its volume. A filling trap keeps its own rate through capacity until the `:fill` event, so the event is located exactly and pinning the trap to capacity adds no water. *`_routed_inflow`; fixed 2026-11-13 with review item C4.*
- **S8.** A commit (and a timepoint sample) never stops at a steady state, and does not integrate the same trajectory twice: it starts from the last saved step of the context's stored prediction (`ctx.prediction`, kept by `_predict_network!`) at or before its time and integrates on from there (with `stop_at_steady = false`), usually a step or two. It never interpolates the prediction: the solver's dense output over a step that crosses a rate jump or that an event cut short was seen 1 % of a trap's capacity off (vakaas with infiltration), adding or removing water at every commit there. Within the prediction's last step the short solve may meet the predicted event slightly early (solver tolerance); any other event, e.g. beyond a prediction that stopped at a steady state, is an error. A commit must not lie past the context's predicted event. *`_state_at`, `_commit_network!`, `_sample_network`; asserted in both (review items C3, C4).*
- **S9.** A trap's net rate below rounding level of its volume (`10⁻¹²` of its capacity per time unit) counts as zero, in the rate function, the `:unspill` condition, the t = 0 fast path and the steady-state tests. Otherwise a full trap whose inflow balances its losses up to the last digits flips full / not full forever on the noise's sign; seen on vakaas with green roofs holding their water in the dry period, which never ended. *`_net_rate`, `_rate_tol` in `networksolver.jl`; test "rounding-level net flow does not flip a trap".*
- **S10.** A trap is full exactly when it has a spill path, and then holds exactly its capacity. Any other trap holds `0 ≤ V < C`, strictly below `C` because a callback condition that starts at exactly 0 never fires, so `:fill` could not be detected; the driver seeds such traps at `prevfloat(C)`. A zero-capacity trap (a parent whose subtraps all spill at its own level, common on flattened terrain) has no volume between full and empty: it enters at `prevfloat(0.0)` and the t = 0 fast path fills it (net inflow > 0) or, for a parent, empties it, exposing its subtraps. Left below 0, it could never fire `:empty` and its subtraps' water never infiltrated (vakaas with infiltration). *`_validate_network` (error), `_t0_fast_path`; test "a zero-capacity parent fills or empties at once".*
- **S11.** A leaf reaching empty is set to exactly 0 by a non-terminating `:floor` condition. Its rate drops to 0 there (the physical floor), and without the condition the solver stepped across that jump at abstol-sized steps for the rest of the solve: a million steps per solve on vakaas with infiltration, leaving traps at negative volumes. A solve that stops for any other reason than an event, a steady state or `tmax` (e.g. `maxiters`) is an error, never read as a steady state. *`_event_conditions`, `_build_event_callback`, `solveDynNetwork!`; test "a leaf draining to empty stops at exactly 0".*
- **S12.** A network trap's ponding loss per wetted cell is the infiltration capacity that cell's rain and through-flow leave unused, `-min(runoff, 0)`, as in base SWIM's `_ponding_infiltration`; the runoff field has already taken the rest off the trap's inflow. Charging the full rate counted the wetted cells' infiltration twice during rain: network traps filled up to 1.2 min late on vakaas at 20 mm/h. The losses follow the runoff frozen at each touch (S5), so a reused context recomputes them. *`_ponding_loss`, `_with_losses`, `_with_flows` in `networksolver.jl`; test "fill_sequence dyn_traps parity: rain with infiltration".*

**Persistence**
- **P1.** Layer volumes cross touches and weather periods only through `nbs_state[id]`. They are written at every harvest and finalize, and read at every rebuild.

**Flow field (`network_watercourses`)**
- **F1.** The field splits each placement's terrain outflow `O_t` over its exits with the solver's shares, i.e. from the baseline without covered traps' spills. At each exit it removes exactly what its own field's footprint delivers there (`_exit_throughput` on that field), so the correction never takes away more than the field carries. *Test "network_watercourses: NBS exit split matches the solver".*
- **F2.** A change the field propagates through a full trap passes on `max(D + d, 0) - max(D, 0)`, where `D` is the trap's net output (what reaches its bottoms less the infiltration room on its submerged footprint), as the solver's full trap spills `inflow - footprint loss` only when positive. A trap short of water takes added flow first; `D` accumulates over the changes through it. Passing the change through unchanged overstated the flow below such traps (vakaas, roof water just after the rain: ~44 m3 too much outflow over 20 min). *`_through_trap!`, `_trap_net_output` in `network_watercourses.jl`; test "network_watercourses: a full trap short of water absorbs added flow first".*

## Known simplifications (marked `# @@@` in code)

- **Cell area.** The footprint area is the cell count, i.e. 1 m² per cell.
- **Evapotranspiration (ET).** Not modelled. The layers' `EVCoeff`/`EVS11` are carried but unused.
- **Ponds and the field.** In `network_watercourses`, a pond's share stays in its trap, even when that trap is full and spilling.

## Open follow-ups (as of 2026-11-12, to be picked up)

1. **Field vs solver: how NBS output is split over exits.** *Fixed 2026-11-13; see F1.* The field's own baseline routes full network traps' spills through the footprint, while the solver's does not (they reach the NBS as draw). The field used to take its shares from its own baseline, so a spill entering a footprint unevenly skewed the field's split away from the solver's. It now takes the solver's shares.
2. **Pond shares in the field.** See "Known simplifications": the share should continue along a full pond trap's spill path.
3. **Remaining NBS items from the 2026-11 review.** They are in `agent/reports/perf_memory_review_2026-11.md`, §6.0 (NBS track), and §4 for the details.
   - C3, C4, C5, C6, C7 and C14 are fixed (2026-11-13; see the report for what each fix did). C4's fix also found and fixed early `:fill` events, which created water at each network fill.
   - Still open from C3: in a *prediction*, a trap in a culvert/NBS network whose rate crosses zero smoothly is taken as settled, so a later event of that network can go unpredicted until the next touch. Commits are exact.
   - **D2:** a touch recommits and rebuilds every component; only the touched ones need it.
   - **D3:** RHS allocations, and a dense finite-difference Jacobian.
   - D1 is done; D2's per-context params caching is done.
4. **Known-broken testsets.** `Pkg.test()` now runs the whole suite (2026-11-13): `MaskingTests` (deprecated) is wrapped so its failures don't stop the run, and the bay-grid drift (1091 vs 1119 traps) and the `raise_buildings!` API mismatch are `@test_broken`. The bay drift still needs explaining.
   - Meanwhile, run the NBS/dynamics tests with `julia --project=. agent/reports/nbs_debug/run_tests.jl`. It isolates the known-broken testsets.
5. **Repro material** is in `agent/reports/nbs_debug/` (see its `README.md`):
   - `setup.jl`, `run1.jl`, `run2.jl`: the SUrbArea synthetic case, run outside the webapp, with a mass-balance check;
   - `zero_vol.jl`: the zero-volume trap count behind S6.
   - The case's inputs are in `<SUrbArea>/debug/` (`synthetic_dtm*`, `synthetic_userdefined*.geojson`), which is not in git.
6. **A change through a full trap in the field: resolved 2026-11-16 (F2).** `_propagate_field_delta!` used to pass a change through a full trap unchanged. That was wrong both for an added flow reaching a trap short of water (it should fill the deficit first) and for a removal larger than the trap's outflow. Both now follow the per-cell rule with the trap's net output as its value.
