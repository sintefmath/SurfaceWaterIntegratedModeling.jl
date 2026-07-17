# `networksolver.jl` — call graph

The multi-trap ODE solver: evolve a `DynNetwork`'s trap volumes forward to the
**first topology-changing event** (`:fill` / `:empty` / `:unspill`), to steady
state, or to a `tmax` cutoff. Network *construction* is separate
(`build_network.jl` traces it, `network_updating.jl` mutates it); this file only
*solves* an already-built network.

**Exports:** `solveDynNetwork!`, `dynNetworkRateFunction!`.
**External callers:** `_predict_network!` / `_commit_network!`
(`network_context.jl`) call `solveDynNetwork!` once per predict / commit.

> Functions are cited by name only, never by line number — line refs rot on
> every edit. Use jump-to-definition.

---

## Big picture

```
        _predict_network! / _commit_network!          ‹network_context.jl›
                        │
                        ▼
                 solveDynNetwork!                      ── the driver
        ┌───────────────┼────────────────┬──────────────────┐
        ▼               ▼                ▼                  ▼
 _validate_network  _build_rate_params  _t0_fast_path    callbacks
 (entry contract)   (static, once)      (skip the ODE)   (event + steady state)
                        │                                   │
                        └───────────────┬───────────────────┘
                                        ▼
                           dynNetworkRateFunction!           ── ODE RHS, every step
                                        │
                                        ▼
                                 _routed_inflow              ── shared hot path
                                        │
                                   _route_flow               ★ core router ★
```

`_routed_inflow` is the single place where "what is arriving where" is decided.
The rate function, the event conditions and the steady-state check all funnel
through it, so they can never disagree about whether a trap is spilling.

---

## Setup — `_build_rate_params`

Everything static for one solve, bundled into `DynNetworkRateParams` so per-step
work is just routing plus a loop over traps.

```
_build_rate_params
├─ _build_trap_geometry ──► _compute_z_vol_tables      ‹fill_sequence.jl›
│                           numregions, subtrapsof     ‹TrapStructure.jl›
│                           «TrapGeometry»
├─ _network_order ────────► _add_culvert_edges!         topological sort over
│                                                       paths (1:np) + traps (np+1:np+nt)
├─ _merge_targets                                       tributary → host path
├─ _build_culvert_plan ───► «CulvertPlan»               endpoint owners + inverts
├─ _path_event_templates                                per-path ordered stops
├─ _footprint_infiltration                              full-trap loss per trap
├─ _path_cell_infiltration / _path_cell_runoff          per-path residual reference
└─ _build_nbs_plan ───────► «NBSElement» «NBSPlan»       ★ NBS only
```

A net with no culverts gets `cvplan = nothing`; a net with no NBS gets
`nbsplan = nothing`. `runoff` is read **only** to build the NBS plan.

---

## The router — `_route_flow`

```
_routed_inflow
├─ water_level                     (culverts only — heads at each end)
├─ _nbs_routing ──► compute_outflow ‹nbs_elements.jl›   ★ signed output diffs
│                   «NBSRouting»
└─ _route_flow                     ★ core ★
   ├─ _path_delivered!             walk one path in event order
   │  ├─ _attenuate_range          the one flow-tracking rule
   │  └─ _culvert_flow ──► culvert_rate ‹culvert_rate.jl›
   │                       └─► _directional_capacity
   └─ _route_trap_node!            culvert deliver/drain, then spill if full
      └─ _culvert_flow ↑
```

Both nodes are visited in `_network_order`, so everything upstream of a node is
final by the time it is reached.

**Mass conservation.** `_path_delivered!` records what a culvert inlet *actually
drew* in `culvert_actual`, and that exact amount — never the requested capacity
— is what the outlet later delivers. A full trap emits
`max(inflow - footprint_infil, 0)` and nothing else.

**Attenuation.** `_attenuate_range` charges flow against the *oblivious
residual* `path_runoff`, not raw infiltration, so a cell the background already
saturated is not charged twice. It is the network-side mirror of `_track_flow!`
(`flow.jl`), which builds the very grid it reads.

---

## Rate function and events

```
dynNetworkRateFunction!
├─ _routed_inflow ↑
├─ wetted_infiltration ──► water_level     accumulating trap: wetted area only
├─ _nbs_state_count
└─ compute_outflow ‹nbs_elements.jl›       ★ per NBS layer

_build_event_callback ──► _event_conditions ──► «EventCondition»
│                         _routed_inflow ↑
└─ «DynNetworkEvent»                       VectorContinuousCallback, LeftRootFind

_build_steadystate_callback      ──► _routed_inflow ↑, wetted_infiltration
_build_steadystate_callback_nbs  ──► dynNetworkRateFunction! ↑    ★ NBS only

_t0_fast_path ──► _routed_inflow ↑, wetted_infiltration
```

Three conditions are monitored and all three must stay: `:fill`
(`capacity - V`), `:empty` (`V`, parents only — a leaf at 0 is just its floor),
`:unspill` (`inflow - loss`).

The steady-state detector is a `DiscreteCallback`, not a continuous one, because
`wetted_infiltration` is a step function of `V` and interpolated mid-step states
would give spurious crossings. `_build_steadystate_callback_nbs` replaces it when
NBS are present: the cheap `_routed_inflow` check cannot express the layer
dynamics, so it evaluates the full rate function instead.

`_t0_fast_path` returns an immediate event when the entry state cannot be held (a
full trap already draining, a parent at its floor with negative inflow), letting
the caller skip integration entirely.

---

## Types

| Type | Role |
|---|---|
| `TrapGeometry` | per-trap capacity, footprint, bottom, infil, `v2z` |
| `CulvertPlan` | per-culvert endpoint owners, inverts, diameter |
| `NBSElement` / `NBSPlan` | ★ per-placement layer model, state offset, carrier paths |
| `NBSRouting` | ★ per-step scratch: signed diffs in, `nbs_draw` out |
| `DynNetworkRateParams` | everything static for one solve (the ODE `p`) |
| `EventCondition` / `DynNetworkEvent` | one monitored event; the event that fired |

## External leaves

| Callee | Defined in |
|---|---|
| `_compute_z_vol_tables` | `fill_sequence/fill_sequence.jl` |
| `compute_outflow` | `dynamics/nbs_elements.jl` |
| `culvert_rate` | `dynamics/culvert_rate.jl` |
| `numregions`, `subtrapsof` | `TrapStructure.jl` |
| `solve`, `VectorContinuousCallback`, `DiscreteCallback` | DifferentialEquations.jl |
| `Graphs.topological_sort_by_dfs` | Graphs.jl |
