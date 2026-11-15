# Large-case profile: the SUrbArea vakaas terrain (2026-11-14)

- **Branch / commit:** `SUrbArea` at `9b956c0`, after the correctness fixes C1–C14 of
  `perf_memory_review_2026-11.md`.
- **Case:** `SUrbArea/WebInterface/model_data/geometry/vakaas_dtm_flattened.jld2`, 3162 × 3689 =
  **11.7 M cells**, as cached by the webapp (small traps flattened).
- **Rain:** 120 mm/h for 30 min, then dry until 120 min (and the last, open-ended period).
- **Machine:** all timings pinned to the P-cores (`taskset -c 0,2,4,6,8,10,12,14`), one run at a
  time.  Peak memory is the process's resident peak (`Sys.maxrss`).
- **Scripts** (`agent/reports/`): `perf_vakaas.jl` (stage 1), `perf_vakaas_webapp.jl` (stages 2–3),
  `perf_vakaas_profile.jl` (sampling profile), `perf_vakaas_networks.jl` (network sizes).  The
  10-roof profile is in `perf_vakaas_profile_K10.txt`.

## 1. Timings

| Stage | Phase | Time | Allocated | Peak RSS |
|---|---|---|---|---|
| 1. bare DTM (11,231 traps, hierarchy depth 98) | `spillanalysis` | 19.1 s | 10.9 GB | 5.1 GB |
| | `fill_sequence` (11,060 events) | 21.1 s | 13.7 GB | 5.7 GB |
| | `trap_states_at_timepoints`, 12 timepoints | 1.6 s | 2.3 GB | – |
| | `watercourses`, one call | 3.0 s | 2.5 GB | – |
| 2. webapp setup (24,046 traps, depth 22, 342 sinks) | rasterize FKB layers | 13.0 s | 27.5 GB | 4.4 GB |
| | `spillanalysis` | 24.5 s | 14.2 GB | 7.0 GB |
| | `fill_sequence`, no NBS (23,638 events) | 31.3 s | 31.2 GB | 7.6 GB |
| 3. green roofs on the largest buildings | `fill_sequence`, 1 roof | 40.5 s | 36.0 GB | 7.9 GB |
| | `fill_sequence`, 10 roofs | **> 50 min, not finished** | | |

Stage 2 is the case as the webapp builds it: its own `_compute_layer_cache` / `_compute_tstruct`
(buildings raised 10 m, waterbodies, courtyard and manhole sinks, culverts).  It has no
infiltration: the vegetation layer (FKB-Grønnstruktur) is hidden in the model data.  The roofs are
`elhadiGreenRoof` with the webapp's default parameters; both layers overflow onto the terrain
(`n_terrain = 2`), so no drawn pipe outlet is needed.

**Base SWIM copes with this size**: about half a minute each for the static and the dynamic part,
peak memory below 8 GB.  None of the memory blow-ups the review warned about appear here:
footprints take 0.03 GB despite the deep hierarchy, and there are only 3 full runoff snapshots
(0.09 GB each, one per weather period; F1 grows with the number of periods).  The flow graph is
the largest resident structure, about 1.4 GB (S1).  Allocation is the main CPU overhead: garbage
collection takes 15–30 % of the time in every phase.

**The network code does not scale with the number of NBS placements**: 1 roof adds 30 % to
`fill_sequence`; 10 roofs make it more than 70 times slower.

## 2. Network sizes with 10 roofs

From `setup_network` on the webapp case (`perf_vakaas_networks.jl`).  The later states take the
filled traps from the run without NBS, so they approximate.  Roof sizes: 2,951–7,642 cells.

At the start (nothing full): **9 networks** (two roofs share one), each with 19–59 traps,
221–569 paths and 446–2,421 path cells.

After the rain (t = 30; 23,836 of 24,046 traps full), they have merged into **4**:

| Network | Roofs | Traps (+ covered subtraps) | Paths | Path cells | Trap footprint cells |
|---|---|---|---|---|---|
| 1 | 8 | 245 (565) | 2,570 | 10,810 | 27,153 |
| 2 | 1 | 103 (297) | 489 | 2,835 | 22,397 |
| 3 | 1 | 19 (63) | 242 | 740 | 1,487 |
| 4 | 0 | 2 (2) | 2 | 4 | 12 |

**Paths far outnumber traps.**  A roof drains over its whole edge, and every cell water leaves the
footprint through is an exit with its own path, so one roof brings 200–570 paths.

## 3. Where the time goes (10 roofs)

Sampling profile of `fill_sequence` over 10 minutes (about 53,500 samples on the main thread).
Inclusive shares:

| Where | Share |
|---|---|
| **Handling a touch of the networks** (`_touch_networks_driver!`) | **≈ 92 %** |
| – re-predicting *every* network (`_rebuild_contexts!` → `_predict_network!`) | 29 % |
| – rebuilding *every* context: rate parameters, NBS plan and baselines (`_make_context` → `_build_rate_params` → `_build_nbs_plan` → `_nbs_baseline`) | 22 % |
| – reconciling base SWIM's spill graph from scratch (`_reconcile_spillgraph!` → `compute_complete_spillgraph`) | 19 % |
| – committing *every* network (`_commit_network!`) | 14 % |
| Rest of `fill_sequence` | ≈ 8 % |

Inside the solves (they make up the predict and commit shares):

| Where | Share of the total |
|---|---|
| Event condition: routes the whole network at every check (`condition` → `_routed_inflow`) | 18 % |
| Walking path cells one by one (`_attenuate_range`) | 12 % |
| Allocation (`GenericMemory`) | 20 % |

**Reading:** almost all of the cost is in how touches are handled, not in the ODE dynamics as
such.  An event touches the networks when one of their traps fires or their inflow changes; with
the terrain filling, that happens at a large share of the 23,000+ events.  Every touch then
commits, rebuilds and re-predicts **all** networks and rebuilds the spill graph, so the cost grows
with (touches × networks × network size).

## 4. What to change, in order

1. **D2, selective touch** (the 92 %).
   - Commit and re-predict only the networks a touch affects; keep the others exactly as they
     are (state, time, prediction, rates).  This is safe: a merge regrows the combined network
     from its seeds and leaves the other network objects untouched; only the owner of the firing
     trap is changed in place.
   - Rebuild a context's rate parameters only when its network changed; when only its inflow
     changed, refresh the external inflow and the NBS baselines of placements whose inflow cells
     changed.
   - Update the spill graph only for the traps whose coverage changed, instead of
     `compute_complete_spillgraph` over all traps.
2. **D3, the solver's per-evaluation cost** (≈ 50 % of the solve time).
   - Event condition: read `:fill` / `:empty` directly from `V`; route only for `:unspill`
     conditions, and only when there are any.
   - Preallocate the rate function's buffers (the 20 % allocation).
   - Composed per-segment attenuation (`max(d + a, b)` composes exactly) instead of walking every
     path cell.
3. **Many exits per roof** (a modelling question, not only speed).  The hundreds of edge paths of
   one roof mostly run side by side to the same few traps.  Merging exits that lead to the same
   place would shrink the networks a lot, but changes how an NBS's outflow reaches the terrain.
4. **Base SWIM** (seconds, not minutes): allocation in `fill_sequence` (≈ 1.3 MB per event; F2, F3,
   F7) and in `spillanalysis` (S2, S5, S6); `watercourses` at 3 s per call (F5: cache its
   topological order), which matters where the viewer calls it per timepoint.

Measure again after each step, with the same scripts.

## 5. Update (later on 2026-11-14): the 10-roof slowdown was a never-ending loop

The 10-roof run did not crawl through the rain: it **never ended**.  In the final, open-ended dry
period, a full network trap (23359) flipped `:unspill` / `:fill` every 0.11 min for days of
simulated time.  Its losses were 0 and its routed inflow 5.5e-16: the roofs hold their water at
their thresholds, so the net rate was rounding noise around zero, and its sign decided the events.
Fixed in `06f8eba`: a net rate below 1e-12 of a trap's capacity per time unit counts as zero
(NBS_INVARIANTS S9).  Found with `fill_sequence`'s extended progress line and
`watch_progress.py`.

With that fix and the D2 steps 1–3, step 6 and commit-from-prediction (pinned, same case):

| `fill_sequence` | Before the work | Now |
|---|---|---|
| no NBS | 31 s | 33 s |
| 1 roof | 40.5 s | 36.7 s (network share 9.2 s → 3.5 s) |
| 10 roofs | never finished | **51.6 s** |

The profile shares of section 3 measured mostly that loop; a profile of the finished 10-roof run is
still to be taken.  Remaining items: `agent/D2_TOUCH_PROPOSAL.md` (step 4, bounded predictions,
the stiff solver's finite-difference Jacobian).

## 6. Update (2026-11-15): base SWIM

A profile of the finished 10-roof run showed base SWIM dominating: the candidate scan in
`_identify_next_status_change!` (≈ 32 %, mostly `all_subtraps_filled`) and `compute_flow` →
`watercourses` (≈ 28 %, topological sort and untyped graph lookups); network touches ≈ 25 %.
Fixed, each with bit-identical event histories (`event_record.jl` / `event_record_compare.jl`):

| Commit | Change | vakaas no NBS | vakaas 10 roofs |
|---|---|---|---|
| (before) | | 34.3 s | 51.5 s |
| `2ee7df3` | candidate scan tests the time bound first; `agglomerations::SimpleDiGraph{Int}` | 21.1 s | 38.6 s |
| `5e0b9f7` | `flowgraph::SimpleDiGraph{Int}` | 12.4 s | 29.3 s |
| `bc33c16` | cell graph sorted once per `fill_sequence` (`cell_order`) | 12.7 s | 25.9 s |

The last one gains little with 3 weather periods; it grows with the number of periods.

`spillanalysis` (2026-11-15): `c051041` builds the spill regions with union-find and the flow graph
directly from each cell's flow target, instead of a cell-sized edge set and
`connected_components`.  vakaas: 12.6 → 7.4 s bare, 15.2 → 8.8 s with the webapp setup.  The
trap structure is identical field by field on 17 cases, including the six flattened SUrbArea grids
bare and with the webapp setup (`tstruct_record.jl` / `tstruct_record_compare.jl`).
