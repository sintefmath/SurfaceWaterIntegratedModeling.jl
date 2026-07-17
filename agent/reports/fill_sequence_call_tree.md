# `fill_sequence` Call Structure — Dynamic Networks & NBS

Call tree of `fill_sequence`, emphasising the parts tied to **dynamic networks**
and **NBS** (nature-based solutions).

**Legend**

| Marker | Meaning |
|--------|---------|
| ◆ | dynamic-network machinery |
| ★ | NBS-specific |
| *(plain)* | shared / static engine |

The deepest shared node is `solveDynNetwork!` (the multi-trap ODE). It is entered
from three places — **predict**, **commit**, **finalize** — and is drawn once at
the bottom rather than repeated inline.

> Functions are cited by name only, never by line number — line refs rot on
> every edit. Use jump-to-definition.

---

## Shape of the thing

`fill_sequence` walks weather periods; within each period it repeatedly asks
"what changes next, and when?", advances to that moment, and emits a
`SpillEvent`. Traps not touched by a dynamic network follow a cheap constant-rate
analytic path. Traps inside a network are **covered**: their volumes come from
the ODE instead, and they are deliberately excluded from the spillgraph and the
changetime machinery.

The network side is owned by a single `NetworkDriver`, built once per weather
period and then **mutated in place** per event (`apply_fill!` / `apply_unfill!` /
`apply_empty!`). The old full-retrace path (`_build_dyn_networks` /
`_touch_networks!`) is retired.

---

## Tree

```
fill_sequence(tstruct, weather_events; time_slack, infiltration,
              dyn_traps, culverts, nbs, verbose)          ── fill_sequence/fill_sequence.jl
│
├─  _compute_z_vol_tables                     volume→level tables, once for all traps
├─  compute_complete_spillgraph                                    ‹spillgraph.jl›
├─  compute_flow                                                   ‹flow.jl›
│   ├─ _compute_initial_rateinfo ─► watercourses                   ‹watercourses.jl›
│   ├─ _is_parent                            merge vs terrain spill
│   └─ _track_flow! ─► _update_Smin_Smax! ─► _ponding_infiltration
│
├─◆ _dyn_seeds                               dyn_trap footprints + culvert ends   ‹network_context.jl›
├─◆ build_network_driver                                           ‹network_driver.jl›
│   ├─  setup_network                                              ‹build_network.jl›
│   │   ├─  _validate_network_inputs          no culvert/seed inside an NBS footprint
│   │   ├─★ _compute_nbs_inflow_outflow_cells!   ── mutates `nbs`: stamps id + cell lists
│   │   ├─  _grow_network_from_seed!          [per seed] trace downstream
│   │   │   ├─ _trace_to_next_trap
│   │   │   ├─ _intersecting_culverts_and_nbs_outlets   (trap: bare ids)
│   │   │   ├─ _intersecting_on_path                    (path: id + position)
│   │   │   └─ _update_pathmap!               truncate where paths meet
│   │   ├─  init_in_counts!                                        ‹network_updating.jl›
│   │   └─  split_network_into_connected_components                ‹network_utils.jl›
│   │       ├─ _coupling_graph ─► _culvert_owner_nodes, _nbs_coupled_nodes ★
│   │       └─ _build_subnetwork ─► _localize_path / _localize_trap
│   └─◆ _rebuild_contexts!
│       ├─ _driver_state0 ─► fill_trap_until    project a newly-absorbed trap to now
│       ├─ _make_context ─► ★_nbs_layer_block   restore layer state from nbs_state
│       └─ _predict_network!  ───────────────────────────┐  first event prediction
│
├─  _set_initial_changetime_estimates                    │
│   ├─ _compute_changetime_estimate            [per trap]│
│   └─◆ _apply_network_changetimeest!          covered traps → exact / Inf
├─  push SpillEvent (period start, full snapshot)        │
│                                                        │
└─  for each WeatherEvent: _fill_sequence_for_weather_event!
    │                                                    │
    ├─  _identify_next_status_change!          branch and bound over candidates
    │   ├─ _update_changetime_estimates! ─► _compute_changetime_estimate
    │   ├─ all_subtraps_filled                 gate: children must be full first
    │   └─ _compute_exact_changetime ─► fill_trap_until ─► _setup_dvdt
    │                                          (exact only when min != max)
    ├─◆ _expand_empty_network_fill_updates     :empty parent → expose draining children
    ├─  update_spillgraph!                     covered traps excluded — ODE owns them
    ├─  setsavepoint!                                              ‹rateinfo.jl›
    ├─  _update_flow! ─► _propagate_amount! ─► _track_flow!        ‹flow.jl›
    │
    ├─◆ _touch_affected_networks!              touched = a member fired, or inflow changed
    │   └─◆ _touch_networks_driver!                                ‹network_driver.jl›
    │       ├─ _capture_fired_kinds            read predictions BEFORE mutating
    │       ├─ _commit_network!  ──────────────────────────────────┤ commit to now
    │       ├─ _harvest! ─► ★_store_nbs_state!  volumes + layers survive the change
    │       ├─ _apply_fired_boundaries!        :empty→0, :unspill→prevfloat(C)
    │       ├─ apply_fill! / apply_unfill! / apply_empty!           ‹network_updating.jl›
    │       │   ├─ apply_fill!   ─► _index_components, _fill_subsumes,
    │       │   │                   grow_spill! ─► _grow_network_from_seed! ↑,
    │       │   │                   _fuse_components! ─► _regrow
    │       │   ├─ apply_unfill! ─► detach_spill! ─► _reroot_or_kill_path!
    │       │   │                                    ─► _promote_tributary! / _kill_path!
    │       │   │                                    ─► _compact!
    │       │   └─ apply_empty!  ─► _fuse_components! ↑   de-subsume a drained parent
    │       ├─ _covered_of                     reads comps (valid pre-rebuild)
    │       ├─ _reconcile_spillgraph! ─► _update_flow! ↑
    │       │      invariant: spillgraph == FULL ∧ ¬COVERED
    │       ├─ _rebuild_contexts! ↑ ─► _predict_network! ──────────┤ re-predict
    │       └─ _apply_network_changetimeest! ↑
    │
    ├─  _collect_amount_updates
    │   ├─ _update_affected_amounts ─► _compute_exact_fill ─► fill_trap_until ↑
    │   └─◆ _network_amount_updates            node → ODE volume; subsumed → capacity
    ├─  _apply_updates!
    └─  push SpillEvent (incremental)
    │
    ├─  _settle_noncovered_at_endtime! ─► _compute_exact_fill ↑
    └─◆ _finalize_networks!  ──────────────────────────────────────┤ advance to endtime
        ├─ _commit_network!                                        │
        └─★ _store_nbs_state!            layer storage → next period
                                                                   │
                                                                   ▼
                              ◆ solveDynNetwork!                   ‹networksolver.jl›
                                (see network_solver_callgraph.md)
```

---

## Notes

- **Three entries to the ODE.** `_predict_network!` solves on a *copy* to
  forecast the next event; `_commit_network!` advances the committed state in
  place; `_finalize_networks!` settles at the period boundary. All three run
  under the context's *cached* `extern_inflow`, which is what makes them
  order-independent of any intervening spillgraph update.

- **Covered vs node.** A network *node* is a trap the ODE integrates. A
  *covered* trap is a node or any descendant subsumed under one — those sit
  static at capacity while their parent is a node. `_net_covered_set` is the
  exclusion set for the constant-rate machinery; `_covered_of` computes the same
  thing from the components when the contexts are mid-rebuild.

- **`state0` seeding.** A newly-absorbed trap's `cur_amounts` entry is generally
  stale, so `_driver_state0` projects it forward with `fill_trap_until`. A
  non-full seed is clamped to `prevfloat(C)` to preserve the solver's
  transitory invariant (`V < capacity` ⟺ `spill_path == 0`).

- **NBS is one type end-to-end.** `DynNBSPlacement` carries the layer model, the
  footprint and the outlets; `setup_network` fills in its cell lists and stamps
  its `id`. Layer state lives in the driver's `nbs_state`, keyed by that id, so
  it survives both structural change and weather-period boundaries.

- **`time_slack`** is threaded through `_fill_sequence_for_weather_event!` but
  not yet read.
