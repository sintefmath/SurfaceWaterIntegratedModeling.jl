# `fill_sequence` — model, logic, and contracts

## 1. What problem it solves

`fill_sequence` answers the question: given a terrain, a series of rain events,
and an infiltration map, *when* does each topographic trap fill or empty, and how
does the overland flow pattern change each time one does?

The output is a compact event log — a `Vector{SpillEvent}` — that records only
the moments when something changes.  Between consecutive events, all rates are
constant, so the water volume in every trap can be computed analytically by
linear interpolation from the nearest recorded `FilledAmount`.

---

## 2. Terrain model

The terrain is a raster grid of elevation values stored in
`TrapStructure.topography`.  Water flows cell-to-cell along the steepest-descent
`flowgraph`.  Cells that have no downstream neighbor are either trap bottoms,
sinks, or domain boundaries.

### Traps and the hierarchy

A **trap** is a topographic depression: a connected set of cells that must fill
to their spillpoint elevation before water can escape.  Traps are numbered
1 … `numtraps(tstruct)`.

Two kinds exist:

| Kind | Indices | Has children? |
|---|---|---|
| **Lowest-level** (= regions) | 1 … `numregions` | No |
| **Parent** | `numregions+1` … `numtraps` | Yes (≥ 2 children) |

The `agglomerations` field is a DAG (actually a forest of trees) encoding the
hierarchy.  `subtrapsof(tstruct, t)` returns the immediate children of trap `t`;
`parentof(tstruct, t)` returns its parent (or `nothing` for top-level traps).

### Volumes

Each trap has two volume quantities:

* `trapvolumes[t]` — total water volume from trap bottom to spillpoint elevation,
  including any volume occupied by sub-traps.
* `subvolumes[t]` — the sum of `trapvolumes` of the immediate children of `t`
  (zero for lowest-level traps, which have no children).

The **own volume** of a trap is `trapvolumes[t] - subvolumes[t]`: the volume of
the trap's *exclusive* area — the portion that sits above all child spillpoints
and below the parent spillpoint.  For lowest-level traps this equals `trapvolumes[t]`
directly.  Only the own volume is tracked by `fill_sequence`; each child is
tracked separately.

### Spillpoints

`spillpoints[t]` records:
* `elevation` — the height at which water escapes the trap.
* `downstream_region` — the index of the trap that receives the overflow
  (positive), or ≤ 0 if it exits the domain.

---

## 3. Physical rates

### Precipitation and infiltration

`WeatherEvent` specifies a rain rate (m/time-unit per cell, uniform or spatially
varying) that remains constant until the next event.  Infiltration is a
time-invariant grid supplied by the caller.

The **runoff** at each cell is:

```
runoff[cell] = precipitation[cell] - infiltration[cell]
```

When runoff is positive, the surplus flows downhill.  When negative, the cell
still has unused infiltration capacity and absorbs water from upstream.

### RateInfo

`RateInfo` is the mutable state object that tracks the current rates everywhere:

| Field | Meaning |
|---|---|
| `runoff[cell]` | Net overland flow (positive) or remaining infiltration capacity (negative) |
| `Smax[t]` | Total unfilled infiltration capacity within trap `t`'s footprint |
| `Smin[t]` | For parent traps: sum of `Smax` of immediate children; zero for lowest-level traps (see below) |
| `trap_inflow[t]` | Net water delivery rate to trap `t` (see below) |

`Smax` is recomputed via `_update_Smin_Smax!` whenever the runoff field changes.

`Smin` is non-zero only for parent traps.  The reasoning is geometric: a
lowest-level trap's wetted footprint can shrink all the way to zero as the water
level falls to zero, so its full `Smax` infiltration capacity is genuinely
available as an uncertainty range.  A parent trap's wetted footprint, by contrast,
always includes the footprints of its children (as long as those children are
active), so the children's infiltration capacity is already committed and must not
be counted twice.  Setting `Smin[t] = sum(Smax[children])` subtracts that
already-counted portion, leaving only the parent's exclusive-area infiltration as
the uncertain part.

### Inflow to a trap

For a **lowest-level trap** `t`, `trap_inflow[t]` is the integral of all positive runoff
draining into its watershed region, accumulated by `watercourses`.

For a **parent trap** `t` that is currently full (i.e. overflowing into its
parent), `trap_inflow[t]` includes the signed overflow from its children:

```
signed_outflow(child) = trap_inflow[child] - Smax[child]
```

A positive `signed_outflow` means the child is overflowing; a negative value
means the child absorbs more via infiltration than it receives, and that deficit
propagates to the parent as a negative contribution.

### Net inflow bounds

Because infiltration within a trap's footprint depends on how much of the
footprint is currently submerged (more water → more infiltration capacity used),
the *exact* net inflow is volume-dependent.  Two bounding rates bracket it:

```
min_net_inflow(t) = trap_inflow[t] - (Smax[t] - Smin[t])
max_net_inflow(t) = trap_inflow[t]
```

`max_net_inflow` assumes all infiltration is already saturated (no additional
capacity).  `min_net_inflow` assumes maximum infiltration operates (all remaining
capacity is available).  These bounds drive the `ChangeTimeEstimate` intervals.

---

## 4. SpillGraph

The `SpillGraph` is a dictionary of directed edges on trap indices.  An edge
`t → target` means trap `t` is full and its overflow goes to `target`.

Rules for what `target` is:

1. If all siblings of `t` are also full → `target = parent(t)`.  Flow is
   accumulated inside the composite parent before escaping.
2. Otherwise → `target = spillpoints[t].downstream_region` (the direct spillpoint
   destination), or `numtraps+1` if outside the domain.

Only **filled** traps have edges.  Unfilled traps have no out-edge.

`update_spillgraph!` keeps the graph consistent incrementally when fill states
change, including redirecting siblings when all of them become full (or when one
empties and siblings must revert to their own spillpoints).

---

## 5. The z-vol tables

For each trap, `_compute_z_vol_tables` precomputes a piecewise-linear lookup
from water volume to water-surface elevation.  This is used inside the ODE solver
to compute the time-varying infiltration rate.

For parent traps, the "floor" of the lookup is the spillpoint elevation of the
highest child (since the child traps must fill first, and the parent's own volume
only starts accumulating above that level).

---

## 6. Top-level algorithm: `fill_sequence`

```
initialise filled_traps, cur_amounts, SpillGraph
for each WeatherEvent:
    compute RateInfo for current fill graph + new rain rate
    compute initial ChangeTimeEstimates
    push full-snapshot SpillEvent
    loop until end of weather period:
        _identify_next_status_change!   → (next_time, fill_updates)
        break if no change before period end
        update filled_traps
        update SpillGraph                → graph_updates
        update RateInfo (flow rerouting) → runoff_updates, inflow_updates
        compute amount_updates for affected traps
        push incremental SpillEvent
finalize cur_amounts at period end
```

### Full vs incremental SpillEvents

The **first** event of each weather period is a full snapshot:

| Field | Type |
|---|---|
| `amount` | `Vector{FilledAmount}` — one entry per trap |
| `filled` | `Vector{Bool}` — one entry per trap |
| `inflow` | `Vector{Float64}` — one entry per trap |
| `rain_rate` | Current rain rate (scalar or matrix) |
| `runoff` | `Matrix{Float64}` — full grid |

All **subsequent** events within a period are incremental:

| Field | Type |
|---|---|
| `amount` | `Vector{IncrementalUpdate{FilledAmount}}` — only changed traps |
| `filled` | `Vector{IncrementalUpdate{Bool}}` — only changed traps |
| `inflow` | `Vector{IncrementalUpdate{Float64}}` — only changed traps |
| `rain_rate` | `Nothing` |
| `runoff` | `Vector{IncrementalUpdate{Float64}}` — only changed cells |

The query functions `amount_at`, `filled_at`, `inflow_at`, `runoff_at`, and
`rainrate_at` reconstruct the full state at any event index by walking backwards
to the nearest full snapshot and then replaying incremental updates.

---

## 7. Change-time estimation: `_identify_next_status_change!`

### Step 1 — Update estimates for affected traps

Only traps whose `trap_inflow` changed (tracked via `RateInfo`'s savepoint
mechanism) need re-estimation.  `_compute_changetime_estimate` returns a
`ChangeTimeEstimate(filling, min, max)`:

* `filling = true` → trap is not yet full and is expected to fill.
* `filling = false` → trap is full and is expected to empty.
* `min`/`max` bracket the time of the transition.

The estimate for a **filling trap** is:
```
min_time = remaining_volume / max_net_inflow
max_time = remaining_volume / min_net_inflow
```
(Infinity if `min_net_inflow < 0`, since the trap may never fill.)

The estimate for a **full trap emptying** depends on whether the trap is
submerged by its parent.  If submerged, the trap cannot start emptying until the
parent drains below the child spillpoint, so the estimate defers to the parent's
drain time.

### Step 2 — Branch-and-bound candidate selection

Candidates are all traps whose `min` changetime is less than the current best
known time, with the additional filter that all sub-traps must be full (parent
traps can only change status once their children are full).

Candidates are examined in order of increasing `min`, each refined to an exact
time via `_compute_exact_changetime` (which integrates the ODE if needed).  As
soon as a candidate's exact time beats the current best, it becomes the new best
and candidates with higher `min` are discarded.

### Step 3 — Temporary fill-flip for subtrap re-estimation

When a candidate trap `c` is confirmed to change status, its children's
`ChangeTimeEstimate`s are recomputed immediately, because their emptying
behaviour depends on whether `c` is submerged.  To compute this correctly, the
function *temporarily* flips `filled_traps[c]` to its new value, recomputes the
child estimates, then flips it back.  The actual update to `filled_traps` happens
in the caller.

---

## 8. Integrating the ODE: `fill_trap_until`

When `Smax ≠ Smin` (i.e. infiltration within the footprint is not yet saturated
and the rate depends on how much of the footprint is wetted), the volume ODE

```
dV/dt = inflow - infiltration(z(V))
```

must be solved numerically.  `fill_trap_until` sets up and solves this ODE using
`DifferentialEquations.jl` with a `VectorContinuousCallback` that terminates
integration at whichever of three conditions fires first:

1. Trap full: `trapvolume - V = 0`
2. Trap empty: `V = 0`
3. Rate sign change: `(dV/dt at t) * (dV/dt at t=0) < 0` — the rate crossed
   zero, meaning the trap is asymptotically approaching a steady-state volume
   and will neither fill nor empty within this weather period.  Because all rates
   are constant between events, this steady state genuinely persists until the
   next weather change, so no further status-change event is registered for this
   trap in the current period.

When `Smax == Smin`, the rate is constant and `fill_trap_until` returns an
analytic answer without invoking the ODE solver.

---

## 9. Flow rerouting: `_update_flow!`

After the `SpillGraph` changes, `_update_flow!` adjusts `RateInfo` by:

1. **Subtracting** the old flow along each changed edge (using `_propagate_amount!`
   with the old edges recorded before they were updated).
2. **Adding** the new flow along each new edge.

`_propagate_amount!` walks the SpillGraph edge-by-edge.  When it crosses a
terrain edge (trap → downstream region, not parent), it calls `_track_flow!` to
adjust the per-cell `runoff` values along the flow path from the trap spillpoint
to the receiving trap bottom.  `_track_flow!` stops early if the added/removed
flow is cancelled by residual infiltration capacity.

---

## 10. Contracts and invariants

### Data structure invariants

**Trap status consistency** — no filled trap may have a non-filled sub-trap.  The
helper `_valid_trap_status` checks this; `compute_complete_spillgraph` asserts it.

**SpillGraph membership** — only filled traps have out-edges.  An unfilled trap
has no entry in `SpillGraph.edges`.

**SpillGraph sibling rule** — if all siblings of a filled trap are also full, all
of them must point to their common parent, not to their individual spillpoints.

**`cur_amounts` timestamp** — at the start of each weather period, all
`cur_amounts[t].time` equal the period start time.  The assertion on line 65 of
`fill_sequence.jl` checks this.

**z-vol monotonicity** — `_compute_z_vol_tables` produces strictly increasing
(z, V) pairs by discarding duplicate z values; the Interpolations lookup requires
this.

### Behavioural contracts

**Event ordering** — the returned `Vector{SpillEvent}` is strictly ordered by
`timestamp`.  Each event's timestamp is ≤ the end of its weather period.

**Snapshot rule** — the first `SpillEvent` of every weather period carries full
(non-incremental) copies of `amount`, `filled`, `inflow`, and `runoff`, plus
a non-`Nothing` `rain_rate`.  All subsequent events within the period are
incremental, with `rain_rate = nothing`.

**`_identify_next_status_change!` mutation contract** — this function mutates
`changetimeest` (refined estimates are stored back), temporarily flips
`filled_traps[cand]` while recomputing child estimates, then restores it.  On
return, `filled_traps` is identical to on entry.  The `cur_amounts` vector is
*not* mutated by this function; amount updates happen in the caller.

**Mass conservation** — water that leaves one element of the network must arrive
at exactly one other element.  The signed outflow of a full trap equals
`trap_inflow[t] - Smax[t]`; this is the quantity routed downstream by
`_track_flow!` / propagated to the parent by direct inflow addition.

**No event at period boundary if period ends** — the loop in
`_fill_sequence_for_weather_event!` breaks as soon as `cur_time > endtime` or
`fill_updates` is empty, so no event is registered beyond the weather period.

---

## 11. Notes and known uncertainties

**`tvolume` override for lowest-level traps in `fill_trap_until`** — the function
first computes `tvolume = trapvolumes[t] - subvolumes[t]`, then for lowest-level
traps overwrites it with `tvolume = trapvolumes[t]`.  Because `subvolumes[t] = 0`
for lowest-level traps, the override is a no-op; the two expressions are
equivalent.

**Stagnation callback** — when condition 3 fires (rate sign crosses zero), no
status-change event is registered for the trap.  Its volume at the moment the
callback fires is recorded as the final state for that weather period.  This is
correct because all rates are fixed between weather events: if the rate crossed
zero, the trap is at a stable steady state and will remain there until the next
`WeatherEvent`.

**`time_slack` parameter** — documented as a tolerance for merging events that
are close in time, but the source includes the note `@@@ NB: Support for this
currently unimplemented`.  It is unclear whether this is a planned future feature
or a stale design artefact.

**No terminal SpillEvent at period end** — `_fill_sequence_for_weather_event!`
finalises `cur_amounts` at `endtime` after the loop but deliberately does not
push a new event.  The finalised amounts are carried into the next weather
period's initial full-snapshot event.  Callers that need amounts at exactly
`endtime` must project forward from the last event using the then-current inflow
rates.
