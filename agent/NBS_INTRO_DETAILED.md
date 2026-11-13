# The NBS module: a detailed walkthrough

This follows one `fill_sequence(...; nbs = ...)` call from the entry point to the results, naming every function involved. It assumes you know SWIM's static analysis and `fill_sequence`.

- For the short version, read `NBS_INTRO_OVERVIEW.md` first.
- **Terms in *italics* on first use are defined either in the glossary of `NBS_INVARIANTS.md` or right where they appear.** The glossary is ordered so that each term uses only earlier ones.
- Invariant labels in parentheses, e.g. (L1), refer to the invariants listed there.

Function names are stable; line numbers are not given because they drift.

---

## 0. Notation used in this walkthrough

| Symbol | Meaning |
|---|---|
| `R` | A cell's baseline runoff value: positive = surface flow, negative = unused infiltration capacity. In the solver it is read from `rateinfo.runoff` (base SWIM's rates, which do not route covered traps' spills); in `network_watercourses` from a plain `watercourses` field (§7.2). |
| `S`, `V` | A layer's storage, in mm (`S`) or as a volume (`V`). The ODE state vector is also called `V` in code, but it is written "the state" here. |
| `A` | Footprint area = number of footprint cells (1 m² per cell). Not the `NBSLayer.A` field, which only SUrbArea's calibration code uses. |
| `O_0` | Baseline throughput of a footprint: rain on it plus baseline flow entering it. |
| `O_t` | Overflow of the top `n_terrain` layers, i.e. what the NBS returns onto the terrain. |
| `E_k` | Overflow (pipe flow) of the k-th piped layer, delivered at `outlets[k]`. |
| `share_e` | Fraction of `O_0` that leaves through exit or pond `e` (they sum to 1). |
| `Q_in` | Input to the top layer, `O_0 + draw`. |
| `nt`, `np` | Number of traps and of paths in a component. |

**Units.**
- Grid quantities follow base SWIM: rain and runoff are a height per time per cell, so a sum over cells is a volume per time with 1 m² cells.
- Layer thresholds are in mm, and `_layer_fluxes` is the only place mm and volume meet.

---

## 1. Inputs

### 1.1 The layer stack: `nbs_elements.jl`

- **`NBSLayer`:**
  - `Smax` (mm): above it, the layer overflows at rate `Kout·(S − Smax)^nout`.
  - `Smin` (mm): above it, the layer percolates to the layer below at rate `Kinf·(S − Smin)^ninf`.
  - The evapotranspiration (ET) fields `EVCoeff`, `EVS11` are carried but unused. `A` is used only by calibration code in SUrbArea, the application that uses SWIM.
- **`compute_outflow(K, n, threshold, S)`:** computes either rate. It is `K·(S − threshold)^n` above the threshold and 0 at or below it (S4).
- **`NBSSystem`:** the layers, top first. Constructors:
  - `puddle`: one layer that overflows;
  - `elhadiGreenRoof`: a soil layer that only percolates, over a drainage layer that overflows;
  - `mantillaRRmodel`: a surface layer that overflows and percolates, over a soil layer that percolates, over a drainage layer that overflows.

### 1.2 The placement: `elements.jl`

```julia
DynNBSPlacement(system, footprint::Vector{Int}, outlets::Vector{CartesianIndex{2}})      # n_terrain = 1
DynNBSPlacement(system, footprint, n_terrain, outlets)
```

- **`n_terrain`:** the number of top layers that overflow onto the terrain.
- The *piped layers* are `_piped_layers(system, n_terrain)`: the lower layers with `Kout > 0`. `outlets[k]` receives the k-th one's overflow.
- The constructor checks that the counts match (I4). For the shipped models with `n_terrain = 1`:
  - **puddle:** 0 outlets;
  - **green roof:** 1 outlet, for the drainage layer;
  - **Mantilla:** 1 outlet, for the drainage layer.

Two fields are filled in later by `_prepare_nbs!` (§3):
- **`link::NBSLink`:** the placement's *inflow cells*, *exits* and *ponds*.
- **`id`:** its position in the `nbs` vector. It labels the stored layer volumes and the `nbs_layers` results.

---

## 2. Entry: `fill_sequence` (`fill_sequence/fill_sequence.jl`)

```julia
fill_sequence(tstruct, weather_events; infiltration, dyn_traps, culverts, nbs)
```

`culverts` and `dyn_traps` (traps you want solved dynamically) are the other inputs that create a dynamic network. They are handled by the same machinery and appear below only where they meet NBS.

1. **Infiltration.** The infiltration grid is copied, and `_nbs_zero_infiltration!` sets it to 0 on every footprint (I3). This is the only change NBS makes to base inputs.
2. **Shared tables and state.**
   - `_compute_z_vol_tables(tstruct)`: per trap, a table mapping stored volume to water level. Computed once. The network uses it for wetted areas and culvert water levels.
   - `nbs_state = Dict{Int,Vector{Float64}}()`: the stored layer volumes, by `id` (P1).
3. **For each weather period:**
   1. `rateinfo = compute_flow(sgraph, rain, infiltration, tstruct)`: plain base SWIM, giving the *baseline*.
   2. `driver = build_network_driver(...; nbs_placements = nbs, precipitation = rain, nbs_state)` (§3).
   3. `_set_initial_changetime_estimates(...)` gives base SWIM's estimates for every trap. Each *component* predicts only its single next event, so `_apply_network_changetimeest!` then sets, for *covered traps*: the exact time for the trap predicted to change, and `Inf` for the others.
   4. A full `SpillEvent` snapshot is pushed. Then `_fill_sequence_for_weather_event!` runs the event loop (§5).

---

## 3. Building the network: `build_network_driver` → `setup_network`

**`build_network_driver`** (`network_driver.jl`) creates the *driver*, the object that owns the network for one weather period:
1. it calls `setup_network`;
2. it records each *network trap*'s committed volume (`vol_by_trapix`);
3. it builds one *context* per component (`_rebuild_contexts!`, §4).

**`setup_network(tstruct, full_traps; dyn_coords, culverts, nbs)`** (`build_network.jl`). Here `full_traps` lists the currently full traps, and `dyn_coords` holds one cell per `dyn_traps` entry. The steps:

1. **`_validate_network_inputs`.** Errors if a footprint is empty or overlaps another (I1), or contains a culvert end, a `dyn_coords` cell or a sink, i.e. a cell where SWIM removes water (I2).
2. **`_prepare_nbs!(nbs, tstruct)`** (`nbs_terrain.jl`). Sets `id` and computes `link = _nbs_link(tstruct, footprint)` from the flow graph:
   - **inflow cells:** the predecessors of footprint cells that lie outside the footprint;
   - **exits:** the successors of footprint cells that lie outside it. If there are none it warns, since a closed basin needs pipes;
   - **ponds:** the `trap_bottoms` inside the footprint.

   This costs O(footprint), and it is the only place a placement is mutated.
3. **Seeds** (*seed*): the `dyn_coords`, both ends of each culvert, and for each NBS its outlets, ponds and exits.
4. **Trace.** `_grow_network_from_seed!` runs from each seed. It follows the flow to the next trap and adds that trap and the connecting *network path* (departing at the seed). If the trap is full, it continues from its spillpoint.
   - A path that reaches a cell already on another path stops there and is recorded as a *tributary* of it.
   - For every new path, `_intersecting_on_path` records each NBS inflow cell the path crosses, as `(placement, position)` in `DynFlowPath.nbs_inlets`. These become `:nbsin` stops (§6).
   - A pond seed traces straight into the trap under the footprint, so that trap becomes a network trap.
5. **Split into components.** `split_network_into_connected_components` builds an undirected *coupling graph* (`_coupling_graph`) and takes its connected pieces. Its vertices are paths and traps. Edges:
   - a trap to its spill path, and a path to its target trap;
   - a tributary to the path it joins;
   - each culvert to the paths and traps holding its ends;
   - for each NBS, everything `_nbs_coupled_nodes` lists (C1): the *carrier paths* of its exits and outlets, every path crossing its inflow cells, and the trap under each pond (`_pond_trap`).

   C1 guarantees that every source of a placement's *draw* lies in its component.

---

## 4. A component's rates: `_make_context` and the NBS plan

`_rebuild_contexts!` creates a context per component with `_make_context(net, tstruct, rateinfo, state0, cur_time, infiltration, z_vol_tables, precipitation, nbs_state)`.

1. **The state (S1).**
   - The trap volumes: `state0(trap)` (`_driver_state0`). A full trap gets its capacity. A trap already in the network gets its committed volume. A trap newly brought in has no committed volume yet, so its last recorded amount is projected to the current time with `fill_trap_until`. A trap that is not full is capped just below capacity.
   - Then `_nbs_layer_block(net, nbs_state)`: each placement's layer volumes in `net.nbs` order, top first. A layer never stored starts empty.
2. **`params = _build_rate_params(...)`.** These are the context's *frozen rates* (S5):
   - per-trap geometry: volume↔level tables, infiltration per cell of the trap footprint;
   - `external_inflow`: per network trap, the baseline `rateinfo.trap_inflow` summed over its *leaf regions*;
   - each path's baseline runoff values `R`, used by the *per-cell rule*;
   - the routing order and each path's ordered stops (`_path_event_templates`):
     - `:trib`: a tributary joins;
     - `:cvin` / `:cvout`: a culvert inlet / outlet;
     - `:nbsin`: an NBS inflow cell;
   - the NBS plan.
3. **`_build_nbs_plan`** (`networksolver.jl`) builds one `NBSElement` per placement:
   - **`O_0` and the shares,** from `_nbs_baseline(nb, tstruct, runoff, precipitation)` (`nbs_terrain.jl`):
     - `O_0 = rain on footprint + Σ max(R, 0)` over the inflow cells (L1);
     - **exit weights** (`_exit_throughput`): the baseline runoff of the footprint cells that drain straight to each exit;
     - **pond weights:** the runoff arriving at each pond;
     - the weights must sum to `O_0` (asserted, L2); dividing by their sum gives the shares.
   - **`terrain_paths`:** `(carrier path, share)` per exit. The carrier is the path whose departure point is the exit (asserted, L3).
   - **`terrain_traps`:** `(trap, share)` per pond, via `_pond_trap`.
   - **`piped_paths`:** `(carrier path, layer)` per outlet.
   - **`state_base`:** where this placement's layers start in the state.
4. **`_predict_network!`** integrates a copy of the state to find the component's next event (§6).

**Why the rates are frozen.** The event loop updates `rateinfo` for an event *before* it commits the interval leading up to that event (§5). Reading the live grid would integrate the past with the event's new flows. That was review item C14 in `agent/reports/perf_memory_review_2026-11.md`, fixed in `d0e296a`.

---

## 5. The event loop: `_fill_sequence_for_weather_event!`

Each iteration:

1. **Next change.** `_identify_next_status_change!` picks the earliest trap status change from the estimates. Network traps take part through their predicted times. Covered traps that are not network traps (full sub-traps inside one) never do.
2. **Base SWIM update.** The changed traps' filled status is set. For traps *not* covered by a network, the spill graph is updated (`update_spillgraph!`), and `_update_flow!` updates the baseline rates. Both are plain base SWIM.
3. **Touch.** `_touch_affected_networks!`: a component is *touched* if one of its traps changed status, or the baseline inflow to one of its network traps' leaf regions changed. If any is touched, `_touch_networks_driver!` does the following:
   1. **Commit** every context to the current time (`_commit_network!`, with its frozen rates). A commit always integrates to that time (`stop_at_steady = false`): it never stops at a steady state, so slow flow is not dropped. Only the event that triggered the touch may fire, at that time; anything earlier is an assertion error (S8).
   2. **Harvest:** trap volumes go into `vol_by_trapix` and layer volumes into `nbs_state` (`_harvest!` → `_store_nbs_state!`).
   3. **Apply the change.** The solver reports each event as `:fill`, `:unspill` or `:empty` (§6), handled by `apply_fill!`, `apply_unfill!` and `apply_empty!` (`network_updating.jl`) respectively. A component can:
      - grow along a newly active spill path;
      - detach traps that are no longer reached;
      - or merge with another component it now reaches.

      Merging, and the restructuring `apply_empty!` and `apply_fill!` need when a parent trap absorbs or releases its sub-traps, rebuild the component from its seeds (`_regrow`); placements keep their `link` and `id`. `apply_unfill!` does not split a component that becomes disconnected: it stays one ODE, which gives the same result. A regrow does split it into connected pieces.
   4. **Reconcile the spill graph** (`_reconcile_spillgraph!`). Base SWIM's spill graph must contain exactly the full traps that are not covered. A trap that has just become covered has its spill withdrawn from the baseline, and one that has just left coverage gets it back.
   5. **Rebuild and predict.** Every context is rebuilt, with new frozen rates from the current `rateinfo`, and predicts again.
4. **Record.** `_collect_amount_updates` gathers the new volumes: those of non-covered traps whose rates (inflow, Smax or Smin) changed (projected at constant rate, under the rates they had until now), those of traps that just changed status (at their capacity), and those of network traps (from the driver). They are stored as incremental updates in a new `SpillEvent`. NBS layer volumes are **not** stored in `SpillEvent`.

**At the end of the period:**
- `_settle_noncovered_at_endtime!` brings the non-covered traps to the end time.
- `_finalize_networks!` commits every context to the end, writes `nbs_state` (P1), and writes the network traps' volumes into `cur_amounts`.

The next period's driver restores the layers from `nbs_state`.

---

## 6. Inside the solver: `solveDynNetwork!` with NBS

`solveDynNetwork!(state, tstruct, p; tmax)` integrates `dynNetworkRateFunction!` from the committed state.

- **Integrator:** `AutoTsit5(Rodas5P(...))`. It uses an explicit method, but switches automatically to an implicit one when the problem becomes *stiff*. NBS layers that drain fast make it stiff.
- **Stopping events** come from traps only:
  - `:fill`: a trap reaches capacity;
  - `:empty`: a parent trap drains to the level of its sub-traps;
  - `:unspill`: a full trap's inflow drops below its losses.
- **Steady state** (only when *predicting*, `stop_at_steady = true`; a commit or a timepoint sample always runs to `tmax`). The integration stops when nothing changes any more:
  - every trap that is not full has a rate of exactly 0, or its rate has changed sign since the start. The sign test catches a trap that settles at a level below capacity, where its wetted infiltration jumps back and forth across one step of the terrain. A small but nonzero rate does *not* count: the trap is still filling or draining, and its event must be predicted;
  - every NBS layer has a rate below a relative tolerance (`_NBS_SETTLED_RTOL` × the largest starting rate). A layer approaches its equilibrium exponentially and never reaches a zero rate; a sign change does not count for it, since a layer can drain and then refill.
  - Known gap: in a culvert or NBS network a trap's rate can cross zero smoothly (a real turn, not a stall) and is then also taken as settled, so a later event can go unpredicted until the next touch.

**One evaluation of the right-hand side** (`dynNetworkRateFunction!(dV, state, p)`):

1. **Route** (`_routed_inflow(state, p)`):
   1. A trap is *spilling* if it is full in this network, i.e. it has a spill path (a *full* trap, in `fill_sequence`'s terms, is one that has filled and not yet unspilled): it then passes its surplus on. Spilling follows the topology, not the volume: a trap that is filling keeps its own rate through capacity until the `:fill` event ends the solve, so that event is located exactly (S7). If there are culverts, each trap's water level is computed.
   2. **Injections**, from the layer volumes only (S2). `_nbs_injections` gives, per placement:
      - `d = O_t − O_0`;
      - `path_inject[carrier] += d·share` per exit;
      - `trap_inject[trap] += d·share` per pond;
      - `path_inject[carrier] += E_k` per outlet.
   3. **The *routing pass*.** `_route_flow` walks paths and traps from upstream to downstream, so whatever feeds a path or trap has been computed by the time it is reached. Every trap starts from `external_inflow + trap_inject`.
      - **At a trap:** apply its culverts; an empty trap's culverts draw at most what reaches it. If it is spilling, push its surplus into its spill path: inflow minus the infiltration over its whole trap footprint.
      - **At a path:** start from `q_head = (flow pushed into it) + path_inject[p]` at the path's head. `_path_delivered!` then walks the path's stops in order:
        - between stops, the *network flow* `d` changes cell by cell as `d ← max(R+d, 0) − max(R, 0)` (`_attenuate_range`);
        - `:trib` adds the tributary's delivered flow;
        - `:cvin` takes the culvert's flow, capped by what is there, and `:cvout` delivers it;
        - **`:nbsin`: the inflow cell's own rule is applied, and what leaves it, i.e. enters the footprint, is added to that placement's draw.** The path carries 0 beyond it.
      - It returns the total inflow per trap (*routed trap inflow*) and `nbs_draw` per placement.
2. **Trap rates.**
   - Spilling: `inflow − trap-footprint infiltration − spill`, which is ≈ 0 while the trap is fed.
   - Otherwise: `inflow − infiltration over the submerged part of the trap footprint` (*wetted infiltration*).
3. **Layer rates.** For each placement, `_nbs_layer_rates!(dV[layers], system, state[layers], A, max(O_0 + draw, 0))`:
   - the top layer gets `Q_in`, and each lower layer gets the percolation of the one above;
   - each layer loses its overflow and its percolation, and the bottom layer's percolation leaves the model;
   - the `max` only absorbs rounding: per inflow cell `c`, its part of `O_0` and its part of the draw add up to `max(R_c + d_c, 0) ≥ 0`, so `Q_in ≥ 0`. The per-placement mass balance is S3.

**Why the downstream trap gets the right amount.**
- Its `external_inflow` already contains the baseline `O_0·share_e`, reduced by the infiltration along the path to the trap.
- The carrier path adds `(O_t − O_0)·share_e` under the same per-cell rule. Together the trap receives `O_t·share_e` reduced by that path, in place of `O_0·share_e` reduced by it.
- The pipe flow arrives separately, through its own carrier path.

---

## 7. Results

### 7.1 Trap and layer states over time: `utils.jl`

**`network_states_at_timepoints(tstruct, seq, times; nbs, culverts, dyn_traps, infiltration)`.** `seq` holds no layer volumes, so this replays the network forward through `seq`:
- at each weather-period start, it rebuilds a driver from that event's full snapshot (`_rebuild_replay_driver`), carrying the layer volumes (`_carry_nbs!`);
- at events that touch a component, it commits, then rebuilds the network from scratch with `setup_network` (not with the incremental `apply_*` updates of the run). Component boundaries can therefore differ from the run's (the run postpones splitting a component after a trap stops spilling); the covered traps and paths should be the same, so the volumes agree up to solver tolerance;
- it samples each requested time with `_sample_network`, which integrates a copy of the state.

The result, per time, is `trap_volumes` (network traps) and `nbs_layers` (per `id`, top layer first).

**`all_states_at_timepoints`** combines the above with base SWIM's `trap_states_at_timepoints`, giving per time:
- `filled`;
- `amounts`, where network traps' volumes replace the base estimates;
- `nbs_layers`;
- `ix`, the index of the last `SpillEvent` at or before that time.

### 7.2 Flow field at a time: `network_watercourses` (`dynamics/network_watercourses.jl`)

1. **Baseline.** Build the baseline `watercourses` field: one copy is modified (`runoff_live`), one is kept untouched (`runoff_orig`). When there are placements, also build the solver's baseline (`runoff_solver`): `watercourses` with only the full traps that no network covers, as in `rateinfo.runoff`.
2. **Network.** Rebuild the network for the given state with `setup_network`.
3. **Culverts.** Take culvert flow `q` at the inlet, capped by the flow there, and add it at the outlet.
4. **NBS.** For each placement:
   - take the shares from `_nbs_baseline` on `runoff_solver`, so the field splits `O_t` over the exits exactly as the solver does (invariant F1);
   - take what the footprint delivers to each exit in this field, `w_e`, from `_exit_throughput` on `runoff_orig`, so placements do not affect each other's baseline. Unlike `runoff_solver`, this field routes full network traps' spills, so `w_e` includes what the solver receives from them as draw;
   - add `O_t·share_e − w_e` starting at each exit cell, and `E_k` starting at each outlet cell (as the solver injects at the head of the carrier path);
   - leave pond shares in their trap. This is a known simplification: if that trap is full, its share should continue along the trap's spill path.
5. **Carrying changes downstream.** `_propagate_field_delta!` applies the per-cell rule. It passes through full traps along their spill paths, and stops when it reaches a footprint cell (after applying the inflow cell's rule, where the solver's `:nbsin` captures) or when nothing is left.
6. **Footprints.** Footprint cells are set to `NaN`.

### 7.3 Maps

`interpolate_timeseries(...; nbs, nbs_color)` samples via `all_states_at_timepoints` and paints the footprints with `nbs_color`.

---

## 8. Things to keep in mind when changing this code

- **Base SWIM must stay NBS-free.** NBS effects belong in the network solver, `nbs_terrain.jl` and `network_watercourses`, not in `watercourses`, `compute_flow` or the spill graph.
- **One owner per physics term.**
  - layer fluxes: `_layer_fluxes`;
  - throughput and shares: `_nbs_baseline`;
  - the per-cell rule: `_attenuate_range` in the solver, and `_propagate_field_delta!` for the field, both with the same formula.

  Don't recompute these elsewhere.
- **A commit must use the context's frozen rates** (S5).
- **Anything that changes what feeds a placement must keep C1 true.** That is, update `_nbs_coupled_nodes` along with it.
- **Known simplifications** are listed at the end of `NBS_INVARIANTS.md`.
