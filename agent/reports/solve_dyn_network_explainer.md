# `solveDynNetwork!` — model, logic, and contracts

## 1. What problem it solves

`solveDynNetwork!` answers the question: given a network of topographic traps,
connecting flow paths, and culverts — all carrying a constant external inflow —
*when* does the first topology-changing event occur, and which trap triggers it?

A *topology change* is any moment that requires rebuilding the network: a trap
fills and starts spilling, a spilling trap runs dry and stops, or a full trap's
net inflow turns negative and it starts draining.  Between such events, the
routing graph is fixed and the ODE governing each trap's volume is smooth (or at
least piecewise-smooth), so it can be integrated numerically.

`solveDynNetwork!` integrates the trap-volume ODE forward, stopping at the first
event or at steady state.  The caller then rebuilds the network for the new
topology and calls again.

This is the lowest-level solver: it knows nothing about weather events, fill
sequences, or the broader `fill_sequence` machinery.  Its caller owns the event
loop.

---

## 2. Data model

### DynNetwork

A `DynNetwork` bundles three element types:

| Field | Type | Role |
|---|---|---|
| `flow_paths` | `Vector{DynFlowPath}` | Terrain paths water travels between traps |
| `traps` | `Vector{DynTrap}` | Topographic depressions that accumulate water |
| `culverts` | `Vector{DynCulvert}` | Engineered pipes crossing between elements |

All indices within a `DynNetwork` are *network-local*: path 1 is `flow_paths[1]`,
trap 1 is `traps[1]`.  The global trap index into `TrapStructure` is stored as
`DynTrap.trap_ix`.

### DynFlowPath

A `DynFlowPath` is an ordered sequence of grid cells from a source (a trap spill
point or an external inflow location) to a target trap (or the domain boundary):

| Field | Meaning |
|---|---|
| `cells` | Ordered grid cells along the path |
| `target_trap` | Network-local index of the receiving trap (0 = exits domain or merges into another path) |
| `culvert_inlets` | `(culvert_index, cell_position)` pairs where culverts draw from this path |
| `culvert_outlets` | `(culvert_index, cell_position)` pairs where culverts deliver to this path |
| `merges` | `(tributary_path_index, junction_cell_index)` for tributary paths that merge in |

Flow is routed in segments between junctions (tributaries, culvert inlets/outlets),
so each tributary is charged infiltration only over the cells it actually travels
through after joining.

### DynTrap

A `DynTrap` wraps one trap from the `TrapStructure`:

| Field | Meaning |
|---|---|
| `trap_ix` | Global index into `TrapStructure.trapvolumes` etc. |
| `spill_path` | Network-local index of the flow path this trap spills into when full; `0` if the trap is not yet full (no spill path); `-1` if it is full and spills straight out of the domain (out-of-domain sentinel, no in-network successor) |
| `culvert_inlets` | Culvert indices that draw from this trap |
| `culvert_outlets` | Culvert indices that deliver into this trap |

### DynCulvert

A `DynCulvert` represents a physical pipe:

| Field | Meaning |
|---|---|
| `inlet`, `outlet` | Grid cell `CartesianIndex` of each pipe end |
| `r` | Barrel radius |
| `Cd`, `Ke`, `Kf`, `Cw` | Hydraulic coefficients (see `culvert_rate`) |

The flow through a culvert is computed by `culvert_rate` (HDS-5 model) from the
current water-surface elevations at the inlet and outlet.  Only downhill flow is
supported; the solver passes `allow_reverse=false`.

### TrapGeometry

`TrapGeometry` is precomputed once per solve from the `TrapStructure` and held
in `DynNetworkRateParams`:

| Field | Meaning |
|---|---|
| `trap_ix` | Global trap index |
| `capacity` | Own storage, net of subtraps |
| `footprint` | Linear indices of the trap's footprint cells |
| `bottom` | Terrain elevation at each footprint cell (raised to sub-trap spillpoint for parent traps) |
| `infil` | Per-cell infiltration rate |
| `zmin` | Minimum bottom elevation (water level of an empty trap) |
| `v2z` | Volume → water-level interpolation object |

The `capacity` and `bottom` conventions mirror `fill_trap_until`: for a parent
trap, the floor is raised to the uppermost child's spillpoint because the parent
can only accumulate water above its children.

---

## 3. Three-state contract

Every trap in the network must be in exactly one of three states on entry to
`solveDynNetwork!`:

| State | Condition | `DynTrap.spill_path` | Trap type |
|---|---|---|---|
| **FULL** | `V == C` exactly | `!= 0` (`> 0`, or `-1` out of domain) | Any |
| **TRANSITORY** | `0 < V < C` strictly | `== 0` | Any |
| **EMPTY** | `V == 0` exactly | — | Leaf only (`trap_ix <= numregions`) |

The rationale:

- A FULL trap is actively spilling, so it must have somewhere for its overflow to
  go: either an in-network downstream path (`spill_path > 0`) or the domain
  boundary (`spill_path == -1`, the out-of-domain sentinel).  Only `spill_path == 0`
  is forbidden for a FULL trap — that is the signature of a trap that filled
  without the network being rebuilt to give it a successor, and `_validate_network`
  rejects it.  The routing (`_route_flow`) and ordering (`_network_order`) layers
  treat any `spill_path <= 0` identically as "no in-network successor" (a FULL
  trap's spill is simply dropped when it exits the domain), so the sentinel matters
  only to state validation.  A parent trap can only be FULL if all its in-network
  children are also FULL (water only accumulates in the parent's own volume once
  the children's spillpoints are submerged).
- A TRANSITORY trap is still filling; it has no outgoing spill path yet.
- An EMPTY state is physically meaningful only for leaf traps.  A parent trap
  at `V = 0` would imply its children are not submerged, meaning the parent should
  not be in the network at all.

`_validate_network` enforces this contract with informative error messages before
the ODE runs.

---

## 4. Geometry helpers: `water_level` and `wetted_infiltration`

Two per-trap geometry queries drive the rate function:

**`water_level(g, V)`** — water-surface elevation at stored volume `V`:
- `V <= 0` → `g.zmin` (empty trap, water at its deepest point)
- `V >= capacity` → `Inf` (sentinel: full trap, entire footprint submerged)
- otherwise → interpolated from `g.v2z`

The `Inf` sentinel is intentional: a full trap's whole footprint is treated as
submerged, simplifying the wetted-infiltration computation below.  For culvert
head calculations, where real elevations are needed and `Inf - Inf = NaN`, the
internal `_surface_level` function clamps to the table top instead.

**`wetted_infiltration(g, V)`** — total infiltration loss rate at volume `V`:
sums `g.infil[k]` over only the footprint cells `k` where `g.bottom[k] <=
water_level(g, V)`.  A partially-filled trap infiltrates only through its wetted
footprint.  A full trap (level `Inf`) infiltrates through its entire footprint.
This matches the `infilfun(trap_bottom .<= z)` pattern in `fill_trap_until`.

---

## 5. Flow routing: `_route_flow`

`_route_flow` computes the total inflow arriving at each trap given the current
spilling status of every trap.  It processes the network in topological order
(upstream to downstream), handling three element types at each node:

**Trap nodes**: accumulate external inflow.  If the trap is spilling (FULL), it
emits `max(inflow - footprint_infiltration, 0)` into its spill path.  If culverts
are present, trap-owned culvert inlets draw from the trap's inflow, and
trap-owned culvert outlets add their delivered flow.

**Path nodes**: route the head flow (from a trap spill or external path inflow)
downstream through segments divided at junctions.  For each segment:

```
delivered_out = max(flow_in - infiltration_on_segment, 0)
```

Tributaries add their delivered output at the exact junction cell, so each
tributary is charged only the post-junction path cells' infiltration — not the
whole path.  Culvert inlets on the path draw from the live flow; culvert outlets
on the path deliver culvert output.

**Mass conservation** is enforced structurally: the actual flow drawn at a culvert
inlet (`culvert_actual[ci]`) equals the flow delivered at its outlet — the same
variable is written at the inlet and read at the outlet.

The routing plan (topological order, tributary lists, culvert ownership) is
precomputed in `DynNetworkRateParams` and reused every rate-function call.

---

## 6. Rate function: `dynNetworkRateFunction!`

`dynNetworkRateFunction!` is the in-place ODE right-hand side.  For each trap `i`:

```
if V[i] >= capacity:                         # FULL (spilling)
    loss  = footprint_infiltration[i]
    spill = max(inflow[i] - loss, 0)
    dV[i] = inflow[i] - loss - spill         # ≈ 0 while adequately fed; < 0 when draining

else:                                         # TRANSITORY (accumulating)
    loss  = wetted_infiltration(geom[i], V[i])
    spill = 0
    dV[i] = inflow[i] - loss
```

`inflow[i]` is the routed total from `_route_flow`.

A floor guard clamps `dV[i]` to 0 when `V[i] <= 0` and the rate is negative
(e.g. a culvert trying to drain an empty trap).  This prevents the ODE from
pushing volumes below zero and breaking the `v2z` interpolation.

---

## 7. Event detection

Two callbacks stop the integration:

### 7.1 VectorContinuousCallback (`cb_topo`) — topology events

A single `VectorContinuousCallback` (with `LeftRootFind`) monitors one condition
per registered event:

| Kind | Condition | Fires when |
|---|---|---|
| `:fill` | `capacity - V[i]` | An evolving trap reaches full |
| `:empty` | `V[i]` | An evolving parent trap reaches empty |
| `:unspill` | `inflow[i] - footprint_infil[i]` | A full trap's net inflow turns negative |

Which traps get which conditions:
- Every **evolving** (non-full) trap: `:fill`
- Evolving **parent** traps (`trap_ix > numregions`): also `:empty` (leaf traps
  draining to zero is their physical floor, not a topology change)
- Every **full** trap: `:unspill`

`LeftRootFind` means a condition that starts at exactly 0 (a full trap with zero
initial net inflow) does not trigger immediately — it waits for a genuine
downward crossing.  This prevents degenerate root-finding at `t=0`.

When a condition fires, the callback records `(kind, trap)` into the shared
`DynNetworkEvent` object and calls `terminate!`.

### 7.2 DiscreteCallback (`cb_ss`) — steady-state detector

A `DiscreteCallback` fires at the first *accepted* ODE step where
`|dV/dt| < abstol` for all evolving traps.  It uses a `DiscreteCallback` (not
`ContinuousCallback`) because `wetted_infiltration` is a step function of `V`:
evaluating the condition at dense-output (interpolated) states mid-step would
detect spurious crossings.

Spilling-veto: if any evolving trap has `V >= capacity` at this accepted step,
the steady-state check is skipped — the VCB `:fill` event takes priority.

Steady state is reported as `(time = Inf, kind = :none)`.

The two callbacks are combined in a `CallbackSet(cb_topo, cb_ss)`.

---

## 8. Driver: `solveDynNetwork!`

```
validate three-state contract
build DynNetworkRateParams (geometry, routing plan)
compute du0 at t=0

if any FULL trap already has du0 < 0:
    return (time=0.0, kind=:unspill)      # t=0 fast path

evolving = [traps with V < capacity]

if evolving is empty:
    return (time=Inf, kind=:none)

if all |du0[i]| <= abstol for evolving:
    return (time=Inf, kind=:none)

build cb_topo (VectorContinuousCallback)
build cb_ss   (DiscreteCallback)
solve ODEProblem(dynNetworkRateFunction!, V0, (0, Inf), p);
       callback=CallbackSet(cb_topo, cb_ss)

state .= sol.u[end]           # update in-place

if kind == :fill:  state[trap] = capacity   # clamp to exact threshold
if kind == :empty: state[trap] = 0.0

return (time, trap=global_trap_ix, kind)
```

### t=0 fast path

Before running the ODE, `solveDynNetwork!` evaluates `du0` once.  If any FULL
trap already has `du0 < 0` (net inflow below its footprint losses at the very
start), it returns `:unspill` immediately.  The ODE is not needed: the event is
already certain at `t=0`.

### State clamping after an event

When the ODE detects `:fill`, floating-point residual may leave `state[trap]`
infinitesimally below `capacity`.  `solveDynNetwork!` clamps it to exactly
`capacity` so the validator passes on the next call without requiring the caller
to handle this residual.  Similarly `:empty` clamps to exactly `0.0`.

FULL traps that did not trigger an event may accumulate small ODE-induced drift
away from exact capacity.  The caller (not `solveDynNetwork!`) is responsible for
clamping those, since the authoritative list of full traps lives in the caller.

---

## 9. Caller protocol after each event

| Event | What `solveDynNetwork!` leaves in `state` | Action before next call |
|---|---|---|
| `:fill` | `state[trap] == capacity` | Add trap to `full_traps`, rebuild network via `setup_network` |
| `:unspill` | `state[trap] == capacity` (unchanged) | Remove trap from `full_traps`, set `state[trap] = prevfloat(capacity)`, rebuild network |
| `:empty` | `state[trap] == 0` | Remove parent + all its children from `full_traps`; set each child `state[child] = prevfloat(C_child)`; rebuild network |
| `:none` | steady state | No further event in this weather period |

The `prevfloat` step moves a just-emptied or just-unspilling trap from the FULL
state to the TRANSITORY state: `V < capacity` strictly, consistent with
`spill_path == 0` in the rebuilt network.

---

## 10. Contracts and invariants

### Three-state contract (entry)

Enforced by `_validate_network` before every call:
- FULL traps have `spill_path != 0` — `> 0` (in-network path) or `-1` (spills out
  of the domain); parent FULL traps have all in-network children FULL.
- TRANSITORY traps have `spill_path == 0`.
- EMPTY traps are leaf traps only.

### Mass conservation

The flow drawn at a culvert inlet is recorded in `culvert_actual[ci]` and is the
*only* quantity delivered at the outlet.  `_route_flow` never uses the *requested*
rate when less water is available — the actual drawn amount is always
`min(culvert_rate(...), current_flow)`.  The spill emitted by a full trap equals
`max(inflow - footprint_infiltration, 0)`, so no water is created or destroyed.

### Topological ordering

The routing graph must be a DAG.  `_network_order` (via `Graphs.topological_sort_by_dfs`)
detects cycles and raises an error.  Culvert edges are included in the DAG with
the constraint that the inlet owner is processed before the outlet owner
(downhill-culvert assumption).  Uphill / reverse culverts would introduce cycles
and are deferred (flagged with `@@@`).

### No mutation except `state`

`solveDynNetwork!` mutates only the `state` vector passed in.  All other inputs
(`tstruct`, `net`, `infiltration`, `inflow`) are read-only.  `DynNetworkRateParams`
is built fresh each call from those inputs.

### Event ordering in `CallbackSet`

`CallbackSet(cb_topo, cb_ss)` places the `VectorContinuousCallback` first.  The
spilling-veto in `cb_ss` additionally ensures that if a trap reaches exactly
`V == capacity` at an accepted step, the topology callback's `:fill` event takes
priority over the steady-state report.

---

## 11. Notes and known limitations

**Culverts are downhill-only.**  The routing assumes that `culvert_rate` returns a
non-negative value (inlet → outlet).  `allow_reverse=false` is passed throughout.
Uphill / reverse-flow culverts are deferred; a cycle in the routing DAG would be
detected by `_topological_order` and raise an error.

**Grid resolution assumed 1 m/cell.**  `DynCulvert`'s convenience constructor
computes barrel length from cell displacement assuming 1 m/cell and flags this
assumption with `@@@`.

**`path_inflow` defaults to zeros.**  External inflow directly onto flow-path cells
(e.g. rainfall falling on path cells rather than trap footprints) can be supplied
via `path_inflow`.  The default of zeros is correct when all rainfall is pre-summed
into the `inflow` vector for each trap.

**`zvt` caching.**  The volume-level tables `_compute_z_vol_tables` are computed
once per `TrapStructure`.  Callers that repeatedly call `solveDynNetwork!` on the
same structure should pass the cached `zvt` to avoid redundant computation.

**No terminal event at steady state.**  When `cb_ss` fires, `solveDynNetwork!`
returns `(time=Inf, kind=:none)`.  The caller must not push a further topology
event for this; it should instead note the steady-state volumes and advance time
to the next weather period boundary.

**Parent traps and the EMPTY event.**  An `:empty` event on a parent trap means
its own volume (above all children's spillpoints) has drained.  At that point all
children must be considered at their individual capacities, not empty — their
accumulated water is still there.  The caller's protocol (`state[child] =
prevfloat(C_child)`) reflects this: children are set TRANSITORY so they can
resume evolving, not EMPTY.
