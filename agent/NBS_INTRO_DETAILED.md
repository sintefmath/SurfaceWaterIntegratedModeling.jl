# The NBS module: a detailed walkthrough

This follows one `fill_sequence(...; nbs = ...)` call from the entry point to the results, naming every function involved. It assumes you know SWIM's static analysis and `fill_sequence`.

- For the short version, read `NBS_INTRO_OVERVIEW.md` first.
- For the authoritative statement of the model and its invariants (I1…P1, cited below), see `NBS_INVARIANTS.md`.

Function names are stable; line numbers are not given because they drift.

---

## 0. Vocabulary

| Term | Meaning |
|---|---|
| **baseline** (*oblivious* in code) | The flow field and trap inflows that base SWIM computes, routing water through footprints as if no NBS existed. `rateinfo.runoff` and `rateinfo.trap_inflow` are baseline quantities. |
| **network** / **component** | A `DynNetwork`: traps and flow paths solved together as one ODE. A component is one connected piece. |
| **network trap** (*node*) | A trap in some `DynNetwork.traps`: one reached by tracing downhill from a seed, continuing through full traps along their spill paths. Its volume is an ODE state, and its events come from the network solver, not from base SWIM's estimates. The set grows when a network trap fills (its spill path is traced further) and shrinks when traps are no longer reached. |
| **covered trap** | A network trap, or a full sub-trap inside one. Base SWIM excludes covered traps from its own changetime estimates and spill graph. |
| **seed** | A cell from which a network is traced downhill (NBS exits, piped outlets, ponds; also culvert ends and `dyn_traps`). |
| **carrier path** | The network path that departs from an exit or outlet cell. NBS output enters at its head. |
| **inflow cells** | Outside cells whose flow-graph successor lies in the footprint. |
| **exits** | Outside cells that footprint cells drain to. |
| **ponds** | Trap bottoms inside the footprint. Water arriving there belongs to that trap. |
| **piped layers** | `_piped_layers(system, n_terrain)`: the layers below the top `n_terrain` with `Kout > 0`. Layer `k` of this list discharges at `outlets[k]`. |
| **`O_0`** | Baseline throughput of a footprint: rain on it plus baseline flow entering it. |
| **`external_inflow`** | Per *network trap*: its baseline inflow, i.e. `rateinfo.trap_inflow` summed over the trap's leaf regions (`_external_inflow`). It includes the throughput of any footprint upstream (after infiltration on the way), because the baseline routes water through footprints. Not the flow entering a footprint. |
| **draw** | Network flow (signed) reaching a placement's inflow cells during one routing pass. |
| **injection** | What a placement adds to the routing: `(O_t − O_0)·share` at exits/ponds, and pipe flow at outlets. |
| **context** | A `DynNetworkContext`: one component with its committed state, the time of that state, and its frozen rates (`params`). |
| **touch / commit / predict** | When an event affects a network: the network is *committed* (its ODE advanced to the event time), updated, rebuilt, and its next event *predicted*. |

**Units.**
- Grid quantities are SWIM's: height per time, summed per cell. This amounts to volume per time with 1 m² cells.
- Layer parameters (`Smin`, `Smax`) are in mm.
- The network stores each layer as a volume over the footprint area `A`, which is the cell count. `_layer_fluxes` is the only conversion point.

---

## 1. Inputs

### 1.1 The layer stack: `nbs_elements.jl`

- **`NBSLayer`:** thresholds `Smin`, `Smax` (mm), and power-law coefficients for overflow (`Kout`, `nout`) and percolation (`Kinf`, `ninf`). The ET fields (`EVCoeff`, `EVS11`) are carried but unused, and `A` is used only by SUrbArea's calibration code.
- **`compute_outflow(K, n, threshold, S)`:** returns `K·(S − threshold)^n` above the threshold and 0 at or below it (S4).
- **`NBSSystem`:** the layers, top first. Constructors:
  - `puddle`: one layer, overflow;
  - `elhadiGreenRoof`: soil, which percolates, over drainage, which overflows;
  - `mantillaRRmodel`: surface, which overflows and percolates, over soil, which percolates, over drainage, which overflows.

### 1.2 The placement: `elements.jl`

```julia
DynNBSPlacement(system, footprint::Vector{Int}, outlets::Vector{CartesianIndex{2}})      # n_terrain = 1
DynNBSPlacement(system, footprint, n_terrain, outlets)
```

The constructor checks `length(outlets) == length(_piped_layers(system, n_terrain))` (I4).

For the shipped models with `n_terrain = 1`, the counts work out as follows:
- **puddle:** 0 outlets; it overflows onto the terrain.
- **green roof:** 1 outlet. The soil layer has `Kout = 0`, and the drainage layer is piped.
- **Mantilla:** 1 outlet. The surface layer overflows onto the terrain, the soil layer has `Kout = 0`, and the drainage layer is piped.

Two fields are filled in later by `_prepare_nbs!`:
- **`link::NBSLink`:** the placement's inflow cells, exits and ponds.
- **`id`:** its position in the `nbs` vector. This keys its stored layer state, and it is what `nbs_layers` results are keyed by.

---

## 2. Entry: `fill_sequence` (`fill_sequence/fill_sequence.jl`)

```julia
fill_sequence(tstruct, weather_events; infiltration, dyn_traps, culverts, nbs)
```

1. **Infiltration.** The infiltration grid is copied, and `_nbs_zero_infiltration!` sets it to 0 on every footprint (I3). This is the only NBS-specific change to base inputs.
2. **Shared tables and state.** `_compute_z_vol_tables(tstruct)` is computed once and shared. `nbs_state = Dict{Int,Vector{Float64}}()` is created; this is the persistent layer storage, keyed by `id` (P1).
3. **For each weather period** (`we`):
   1. `rateinfo = compute_flow(sgraph, we.rain_rate, infiltration, tstruct)`: plain base SWIM. This is the **baseline**.
   2. `driver = build_network_driver(...; nbs_placements = nbs, precipitation = we.rain_rate, nbs_state)` (§3).
   3. `_set_initial_changetime_estimates(...)`: base estimates for every trap. For traps covered by a network they are then replaced by the network's prediction (`_apply_network_changetimeest!`): exact for the predicted trap, `Inf` for the rest.
   4. A full `SpillEvent` snapshot is pushed, then `_fill_sequence_for_weather_event!` runs the event loop (§5).

---

## 3. Building the networks: `build_network_driver` → `setup_network`

**`build_network_driver`** (`network_driver.jl`):
1. It calls `setup_network`.
2. It records each network trap's committed volume (`vol_by_trapix`).
3. It builds one context per component (`_rebuild_contexts!`, §4).

**`setup_network(tstruct, full_traps; dyn_coords, culverts, nbs)`** (`build_network.jl`):

1. **`_validate_network_inputs`.** Errors if a footprint is empty or overlaps another (I1), or contains a culvert end, a `dyn_coord` or a sink (I2).
2. **`_prepare_nbs!(nbs, tstruct)`** (`nbs_terrain.jl`). Sets `id` and computes `link = _nbs_link(tstruct, footprint)`:
   - **inflow cells:** the flow-graph *in*-neighbours of footprint cells that lie outside it;
   - **exits:** the flow-graph successors of footprint cells that lie outside it (if there are none, it warns: a closed basin needs pipes);
   - **ponds:** `trap_bottoms` inside the footprint.

   The cost is proportional to the footprint. This is the only place a placement is mutated.
3. **Seeds.** The `dyn_coords`, both ends of each culvert, and for each NBS its `outlets`, `link.ponds` and `link.exits`.
4. **Trace.** `_grow_network_from_seed!` runs from each seed: follow the flow to the next trap, add the trap and the path, and continue from its spill point if it is full.
   - For every new path, `_intersecting_on_path` records each NBS inflow cell the path crosses, as `(nbs index, position)` in `DynFlowPath.nbs_inlets`. These become `:nbsin` stops for the router.
   - A pond seed traces straight into the trap under the footprint, so that trap joins the network.
5. **Split into components.** `split_network_into_connected_components` → `_coupling_graph` links what must be solved together. For NBS, `_nbs_coupled_nodes` adds (C1):
   - the paths departing from the placement's exits and outlets (its carriers);
   - every path crossing one of its inflow cells (its feeders);
   - the trap under each pond (`_pond_trap`).

   This guarantees that every source of a placement's `draw` is in the same component, so one routing pass collects all of it.

---

## 4. A network's rates: `_make_context` and the NBS plan

`_rebuild_contexts!` creates a context per component with `_make_context(net, tstruct, rateinfo, state0, cur_time, infiltration, z_vol_tables, precipitation, nbs_state)`.

1. **State vector** (S1): trap volumes (`state0(g)` from the committed volumes), then `_nbs_layer_block(net, nbs_state)`, which holds each placement's layers in `net.nbs` order, top first. A layer starts empty if it was never stored.
2. **`params = _build_rate_params(...)`**, frozen on the context (S5). It holds the trap geometry, the external inflow (summed from the baseline `trap_inflow` of the network traps' leaf regions), the per-path baseline runoff used by the attenuation rule, the routing order and stop lists (`_path_event_templates`: tributaries, culverts, `:nbsin`), and the NBS plan.
3. **`_build_nbs_plan`** (`networksolver.jl`) builds one `NBSElement` per placement:
   - **`O_0` and the shares:** from `_nbs_baseline(nb, tstruct, runoff, precipitation)` (`nbs_terrain.jl`).
     - `O_0 = rain on footprint + Σ max(runoff, 0)` over the inflow cells (L1).
     - **Exit shares:** the baseline runoff of footprint cells draining straight to each exit (`_exit_throughput`).
     - **Pond shares:** the runoff arriving at each pond.
     - The weights must sum to `O_0` (asserted, L2); they are then normalised.
   - **`terrain_paths`:** `(carrier path, share)` per exit. The carrier is the path whose `departure_point` is the exit (asserted, L3).
   - **`terrain_traps`:** `(trap, share)` per pond, via `_pond_trap`.
   - **`piped_paths`:** `(carrier path, layer)` per outlet, pairing `outlets` with `_piped_layers`.
   - **`state_base`:** the offset of this placement's layers in the state vector.
4. **`_predict_network!`** integrates a copy of the state to find the next event (§6).

Freezing `params` matters. The event loop updates `rateinfo` for an event *before* it commits the interval leading up to it. Reading the live grid would integrate the past with the new flows; this was bug C14, fixed in `d0e296a`.

---

## 5. The event loop: `_fill_sequence_for_weather_event!`

Each iteration:

1. **`_identify_next_status_change!`** picks the earliest trap status change. Network nodes take part through their exact predicted times. Traps *subsumed* under a network parent are excluded.
2. **Spill graph.** Fill updates are applied. The spill graph is updated only for traps *not* covered by a network (`update_spillgraph!`), and `_update_flow!` updates the baseline rates. Both are plain base SWIM.
3. **`_touch_affected_networks!`.** A network is *touched* if one of its traps fired, or if the baseline inflow to one of its leaf regions changed. If any network is touched, `_touch_networks_driver!` does the following:
   1. Commit every context to the current time (`_commit_network!`, using the frozen `params`).
   2. Harvest the volumes into `vol_by_trapix` and the layer states into `nbs_state` (`_harvest!` → `_store_nbs_state!`).
   3. Apply the structural change: `apply_fill!` / `apply_unfill!` / `apply_empty!` (`network_updating.jl`). A network can grow along a new spill path, detach a trap that stopped spilling, or merge with another. Merging regrows the network from its seeds (`_regrow`); placements keep their `link` and `id`.
   4. Reconcile the spill graph with the new network coverage.
   5. Rebuild every context (new `params` from the current `rateinfo`) and predict again.
4. **Record the event.** Amount updates are collected; network traps' volumes come from the driver. A `SpillEvent` is pushed. NBS layer volumes are **not** part of `SpillEvent`.

**At the end of the period:**
- `_settle_noncovered_at_endtime!` projects the non-network traps.
- `_finalize_networks!` commits every context to the period end, writes `nbs_state` (P1) and the network traps' volumes into `cur_amounts`.

The next period's driver restores the layers from `nbs_state`.

---

## 6. Inside the solver: `solveDynNetwork!` with NBS

`solveDynNetwork!(state, tstruct, p; tmax)` integrates `dynNetworkRateFunction!`.
- **Integrator:** `AutoTsit5(Rodas5P(...))`. The NBS layers can be stiff.
- **Events:** traps only (`:fill`, `:empty`, `:unspill`).
- **Steady-state stop:** `_build_steadystate_callback_nbs` also requires the layer rates to settle.

One rate evaluation:

1. **`_routed_inflow(V, p)`.**
   1. `spilling[i] = V[i] ≥ capacity`; trap water levels are computed if there are culverts.
   2. **`_nbs_injections(V, plan, nt, np)`**, from the layer state only (S2):
      - `diff = _nbs_terrain_outflow(...) − O_0` (the overflow of the top `n_terrain` layers, minus the baseline);
      - `path_inject[carrier] += diff·share` per exit, and `trap_inject[trap] += diff·share` per pond;
      - `path_inject[carrier] += overflow` of each piped layer.
   3. **`_route_flow(...)`** walks `p.order`, so every node's upstream side is final when the node is visited.
      - **Trap:** apply its culverts; if it is spilling, push `inflow − footprint loss` into its spill path.
      - **Path:**
        - `head = flow from upstream + path_inject[p]`, then `_path_delivered!` walks the path's stops;
        - between stops each cell applies `d ← max(V+d,0) − max(V,0)` against the baseline runoff;
        - `:trib` adds a tributary, `:cvin`/`:cvout` move culvert water;
        - **`:nbsin` moves the whole passing flow into `nbs_draw[placement]`** and zeroes the carrier.
      - It returns `(trap_inflow, nbs_draw)`.
2. **Trap rates:** if spilling, `inflow − footprint loss − spill` (≈ 0 while fed); otherwise `inflow − wetted infiltration`.
3. **Layer rates:** for each placement, `_nbs_layer_rates!(dV[layers], system, V[layers], A, max(O_0 + draw, 0))`.
   - Top input is `I_1`, and each lower layer receives the percolation of the one above.
   - Each layer loses its overflow and percolation; the bottom layer's percolation leaves the model.
   - The `max` only absorbs rounding, since `I_1 ≥ 0` analytically. The mass balance per placement is S3.

**Why the downstream trap gets the right amount:**
- Its `external_inflow` already holds the baseline `O_0·share`, as it arrived after the path's infiltration.
- The carrier adds `(O_t − O_0)·share`, giving `O_t·share` in total.
- The pipe flow arrives separately via its own carrier.

---

## 7. Results

### 7.1 Trap and layer states over time: `utils.jl`

- **`network_states_at_timepoints(tstruct, seq, times; nbs, culverts, dyn_traps, infiltration)`.** Since `seq` holds no layer states, this replays the networks forward through `seq`:
  - at each weather boundary it rebuilds a driver from the snapshot (`_rebuild_replay_driver`), carrying the layer storage (`_carry_nbs!`);
  - at touched events it commits and rebuilds;
  - it samples each requested time with `_sample_network`, which integrates a copy of the state.

  The result, per time, is `trap_volumes` (network traps) and `nbs_layers` (per `id`, top layer first).
- **`all_states_at_timepoints`.** Combines the above with base SWIM's `trap_states_at_timepoints`, giving `(filled, amounts, nbs_layers, ix)` per time.

### 7.2 Flow field at a time: `network_watercourses` (`dynamics/network_watercourses.jl`)

1. **Baseline.** Build the baseline `watercourses` field (`runoff_live`) and keep an untouched copy (`runoff_orig`).
2. **Networks.** Rebuild the networks for the given state with `setup_network`.
3. **Culverts.** Draw `q` at the inlet (capped by the flow actually there) and deliver it at the outlet.
4. **NBS.** For each placement:
   - get `O_0` and the shares from `_nbs_baseline` on `runoff_orig`;
   - add `(O_t − O_0)·share` downstream of each exit and the pipe flow downstream of each outlet;
   - pond shares are left in the trap. This is a known simplification: if that trap is full, its share should continue along the trap's spill path.
5. **Carrying changes downstream.** `_propagate_field_delta!` uses the same per-cell rule as the solver. It goes through full traps along their spill paths and stops at another footprint or when nothing is left.
6. **Footprints.** Footprint cells are set to `NaN`.

### 7.3 Maps

`interpolate_timeseries(...; nbs, nbs_color)` samples via `all_states_at_timepoints` and paints the footprints with `nbs_color`.

---

## 8. Things to keep in mind when changing this code

- **Base SWIM must stay NBS-free.** NBS effects belong in the network solver, `nbs_terrain.jl` and `network_watercourses`, not in `watercourses`, `compute_flow` or the spill graph.
- **One owner per physics term.** Layer fluxes live in `_layer_fluxes`; throughput and shares in `_nbs_baseline`; the downstream rule in `_attenuate_range` (solver) and `_propagate_field_delta!` (field), which share the same formula. Don't recompute these elsewhere.
- **A commit must use the context's frozen `params`** (S5).
- **Anything that changes what feeds a placement must keep C1 true.** That is, update `_nbs_coupled_nodes` along with it.
- **Known simplifications:** cell area 1 m², no evapotranspiration, and the pond-share limitation in `network_watercourses`. See the end of `NBS_INVARIANTS.md`.
