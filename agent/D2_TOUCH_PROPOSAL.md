# Proposal: cheaper network touches (review item D2, and related)

- **Date:** 2026-11-14.  **Status:** proposal, for review; nothing implemented.
- **Basis:** the vakaas profile, `agent/reports/perf_vakaas_2026-11-14.md`.  With 10 green roofs,
  `fill_sequence` is more than 70 times slower than with one, and about 92 % of its time is spent
  handling *touches* of the dynamic networks.
- **Goal:** make each touch cheap, so that the cost holds up when networks merge into one large
  network (as on vakaas, where 8 of 10 roofs end up in one); and touch only the affected
  networks, for when placements sit in separate catchments.

## 1. What a touch does today

An event *touches* the networks when one of their traps fires or the baseline inflow to one of
their leaf regions changes (`_touch_affected_networks!`).  Then `_touch_networks_driver!`, for
**every** network, not just the touched one:

| Step | Code | Share of the 10-roof run |
|---|---|---|
| commit it to the event time | `_commit_network!` | 14 % |
| apply the structural change (grow, detach, merge) | `apply_fill!` / `apply_unfill!` / `apply_empty!` | small |
| rebuild base SWIM's spill graph from scratch, if coverage changed | `_reconcile_spillgraph!` → `compute_complete_spillgraph` | 19 % |
| rebuild its context: rate parameters, NBS plan and baselines | `_make_context` → `_build_rate_params` | 22 % (14 % in `_nbs_baseline`) |
| predict its next event, integrating to the end of the weather period | `_predict_network!` | 29 % |

The rate parameters (`_build_rate_params`) fall into three groups:

| Group | Contents | Depends on |
|---|---|---|
| structure | routing order, merge targets, path event templates, culvert plan, trap geometry, footprint infiltration | the network's structure and the period's infiltration grid |
| runoff | per-path oblivious residuals (`path_runoff`), NBS plan (`O_0`, exit and pond shares) | the oblivious runoff field along the paths and around the footprints |
| inflow | `external_inflow` per network trap | `rateinfo`, per leaf region |

## 2. Proposed steps, in order

Each step is a separate commit, with its tests, and is measured on vakaas with 1 and 10 roofs
(`perf_vakaas_webapp.jl`, pinned) before the next one starts.

### Step 1. NBS baselines without per-call overhead (part of the 14 %)

`_exit_throughput` builds a `Set` of the footprint and a `Dict` of exit cells, and walks every
footprint cell's downstream neighbour, on every call; with roofs of 3,000–7,600 cells and hundreds
of exits that is most of `_nbs_baseline`.  The mapping *footprint cell → exit* depends only on
the terrain, so compute it once per placement (in `_nbs_link`, stored on the `NBSLink`) and sum
over it per call.

- **Results:** identical (the same sums, in the same order).
- **Risk:** low.

### Step 2. Reuse a context's rate parameters when its network has not changed (≈ 22 %)

On a touch, a network whose *structure* did not change (no trap of it fired; it was not grown,
detached or merged) keeps its structure group.  Only the runoff and inflow groups are refreshed:

- `external_inflow`: recomputed (cheap);
- `path_runoff`: only for paths that cross a cell whose runoff changed in this event
  (`getrunoffupdates`);
- NBS plan: only for placements whose inflow cells or exit routes cross such a cell.

A network whose structure changed is rebuilt as today.

- **Results:** identical (the refreshed values equal what a rebuild computes).  A test compares
  the refreshed parameters with a full rebuild after a series of events.
- **Risk:** low to medium: the change detection must not miss a cell.  Checked by that test, and
  by an assert in debug runs that compares with a full rebuild.

### Step 3. Update the spill graph only where coverage changed (≈ 19 %)

`_reconcile_spillgraph!` computes a complete spill graph over all traps and diffs it with the
current one.  Only traps whose coverage changed (entered or left a network), and their spill
edges, can differ: update those edges directly and pass the same diff to `_update_flow!`.

- **Results:** identical (the same edge diff).  Invariant S6 is already asserted at period start.
- **Risk:** low to medium (spill edges of parents/children of a trap whose coverage changed).

### Step 4. Touch only the affected networks (helps when networks stay separate)

Commit, change and re-predict only:
- networks that were touched (a trap fired, or their inflow changed);
- networks merged into a touched one (committed before they are regrown together).

All others keep their state, time, prediction, window and rates unchanged.  This is safe: a
merge regrows the combined network from its seeds and leaves the other network objects untouched;
only the owner of the firing trap is changed in place (checked in `apply_fill!` /
`_fuse_components!`).

- **Results:** an untouched network is no longer committed at every touch, so its solver takes
  different steps; results change at the solver's tolerance.  The network tests compare with
  tolerances already (`PARITY_TOL`).
- **Risk:** medium: every place that reads a context's state must commit it first when it is
  behind (`_finalize_networks!` and `_sample_network` already do).

### Step 5. Predict only as far as needed (≈ 29 %, and most of the commit share)

Today a prediction integrates to the end of the weather period, which for the last period is
practically forever, but the next event elsewhere usually comes soon, and the prediction is then
thrown away.  Proposed: predict only up to a *horizon*, the earliest lower bound of the next
constant-rate event (the branch-and-bound's candidate minimum).  If the network's own event
comes before that, it is found as today; if not, the prediction stops at the horizon with a
**re-predict marker**: an event of its own at the horizon, which (if nothing touched the network
before) commits it there and predicts on.

- **Results:** change at solver tolerance (different windows); the commit-span invariant (S8)
  holds, since commits use the prediction's window.
- **Risk:** medium: the marker is a new kind of event in the event loop
  (`_identify_next_status_change!` and the driver), and the horizon must be a true lower bound.
- **Likely the largest single saving once networks have merged.**

### Step 6. Cheaper solver evaluations (D3; ≈ 50 % of what remains of the solves)

- **Event condition:** `:fill` and `:empty` read `V` directly; route the network only for
  `:unspill` conditions, and only if there are any (today: every check routes, 18 %).
- **Allocation:** preallocate the rate function's buffers in the parameters (20 %).
- **Path walking:** compose each path segment's per-cell rule into one `max(d + a, b)` map,
  exact up to rounding, instead of walking every cell (12 %).

- **Results:** identical up to rounding.
- **Risk:** low to medium.

## 3. Exits per roof: deferred (2026-11-14)

A roof drains over its whole edge, and every edge cell water leaves through has its own path down
the terrain: 200–570 per roof on vakaas, 2,570 paths in the merged network.  Paths that reach the
same cell are already joined there (the later one becomes a tributary), so they stay separate
only while they really run on different cells.

Merging exits would send the water of several edge cells down one path: only that path's cells
would infiltrate, and the overland flow would change.  **Decision: not now; possibly revisited
later.**  The cost of many paths is addressed instead by step 6, which makes a path cost a fixed
amount per segment rather than per cell, with the same routing and infiltration.

## 4. What stays the same

- Mass conservation and the NBS model; the exactness of commits (S8).
- `TrapStructure` is not touched (no new fields).
- Base SWIM stays NBS-oblivious.

## 5. Order and decision points

Updated 2026-11-14, after steps 1–3 and a second profile of the 10-roof run.  Steps 1–3 removed
what they targeted (about 35 points of the old time: baselines, context rebuilds, spill-graph
reconcile); the 10-roof run still takes more than 93 minutes, and the ODE solves are now 64 % of
the time (predict 45 %, commit 24 %), with routing (`_routed_inflow`) at 46 %, 28 % of it from
the event condition.

1. Steps 1–3: done.
2. **Step 6 (D3)**: cheaper routing and solver evaluations - now the largest share.
3. **Commit from the stored prediction** (new): a commit integrates the same trajectory the
   prediction already computed (same window and steps, S8).  Keep the prediction's solution on
   the context and read the committed state at `T_commit` from it, instead of solving again;
   timepoint sampling likewise.  Roughly halves the solves, and replaces the `horizon` /
   `stop_at` machinery.
4. Step 4: selective touch (prepared on branch `d2-step4`); matters for separate networks.
5. **Kept on the list, deferred:** bounded, on-demand predictions - give each network a cheap
   lower bound on its next event and solve only when that bound could be the next event, as for
   constant-rate traps.  The largest potential saving on predictions, but it reaches into the
   branch-and-bound and needs bounds that are valid for every event kind; step 5 (a prediction
   horizon) is the simpler variant.  Decided after measuring 2–4.

Each step needs your go-ahead.
