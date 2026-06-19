# Dynamic Network Solver — Implementation & Testing Report

**Scope:** implementation of `solveDynNetwork()` and its supporting machinery for
"solving" a `DynNetwork` (evolving the water content of terrain traps until a
topology-changing event), as specified in `prompt_refined.org`.

**Location:** all new code lives in `src/dynamics/networksolver.jl` (≈ 617 lines),
included from the package module after `dynamics/elements.jl` and `utils.jl`.
Durable tests are in `test/dynamics_test.jl`.

**Status:** complete. All task-tree items in `prompt_refined.org` are `DONE`. The
public entry points `solveDynNetwork` and `dynNetworkRateFunction!` are exported.

---

## 1. What the solver does

A `DynNetwork` describes the *transient* part of terrain flow: a set of lakes
(traps) connected by flow paths ("rivers"), where some elements' state (water
content) is not constant because of inflow/outflow. "Solving" the network means
integrating the lakes' stored volumes forward in time until the **first event**
that changes the network topology, at which point an outer loop (out of scope
here) would rebuild the network and continue.

In the current scope — **lakes and rivers only, constant inflow** — there are two
live outcomes:

| Outcome | Meaning | Returned |
|---|---|---|
| **fill** | a leaf lake reaches capacity and begins to spill (type-1 event) | `(time, trap, :fill, state)` |
| **steady state** | the network settles with no further event | `(Inf, 0, :none, state)` |

`solveDynNetwork` returns a named tuple `(; time, trap, kind, state)`:
`time` is the event time (`Inf` for steady state), `trap` is the **global** trap
index in the `TrapStructure` (0 for steady state), `kind ∈ {:fill, :empty, :none}`
(plus dormant `:unspill`), and `state` is the trap-volume vector at the event.

---

## 2. Architecture

The driver is intentionally thin; the work lives in four reusable pieces, built
and tested bottom-up so each was verified before anything depended on it:

```
solveDynNetwork()                         driver: params → ODEProblem → solve → map result
├── dynNetworkRateFunction!()             ODE right-hand side (dV/dt)
│   ├── geometry helpers                  capacity, volume↔level, wetted-area infiltration
│   └── flow routing                      total inflow per trap (lossy paths, merges)
└── event-detection callback              first-event VectorContinuousCallback
```

Nesting expresses dependency (a parent depends on its children); siblings are
independent. The geometry helpers and flow routing are genuinely independent:
routing only ever needs *whole-footprint* infiltration (full traps are fully
submerged), so it never touches the volume↔level geometry.

### Reuse of existing machinery

A guiding constraint was to build on what the package already has rather than
duplicate it:

- **`_compute_z_vol_tables`** (fill_sequence) — per-trap volume↔level tables.
- **`_setup_dvdt`'s `infilfun(trap_bottom .<= z)` pattern** — wetted-area
  infiltration (sum over submerged cells).
- **`waterheight`** — used as an independent oracle in tests.
- **`_update_runoff!`'s throughput-capped logic** — the model for lossy path
  routing, lifted from the grid to the network abstraction.
- **`fill_trap_until`'s `VectorContinuousCallback` `[full, empty, stagnation]`
  triple** — generalised from one trap to the network.

---

## 3. Components in detail

### 3.1 Geometry helpers

`TrapGeometry` precomputes, per trap in the network, everything the rate function
repeatedly needs: `capacity`, `footprint` (linear cell indices), `bottom` (terrain
elevation per footprint cell), `infil` (per-cell infiltration rate), `zmin`, and
`v2z` (the volume→level map). Built by `_build_trap_geometry`.

Two query functions:

- **`water_level(g, V)`** — surface elevation at stored volume `V`: `zmin` when
  empty, `Inf` when full (so the whole footprint reads as submerged), interpolated
  in between. Mirrors the `z` computation in `_setup_dvdt`.
- **`wetted_infiltration(g, V)`** — infiltration loss = the per-cell rate summed
  over only the **currently-submerged** cells (`bottom[k] <= water_level`).

**State convention:** the state is each trap's stored volume *net of its
subtraps*. `capacity = trapvolumes − subvolumes` for a parent trap (and just
`trapvolumes` for a lowest-level region). For a parent trap, the footprint cells
lying over a subtrap have their `bottom` raised to the subtrap's spillpoint
elevation, exactly as `_compute_z_vol_tables` does — so a submerged subtrap reads
as wetted and contributes its infiltration once the parent holds any water.

### 3.2 Flow routing (`_route_flow`)

Produces the total inflow arriving at each trap = the caller's constant external
inflow **plus** everything routed in from upstream. Each *spilling* (full) trap
emits `max(inflow − footprint_infiltration, 0)` into its spill path; paths lose
`path_infiltration` (floored at zero) in transit; the remainder is delivered to
the path's target trap, into the path it merges into, or out of the domain.
Non-spilling (accumulating) traps pass nothing on.

The router is **pure topology + numbers** — it takes `footprint_infil` and
`path_infil` as arguments rather than reaching into the `TrapStructure`, which
makes it directly testable on hand-built networks. Helpers
`_footprint_infiltration` and `_path_infiltration` compute those sums from a real
`TrapStructure` + infiltration grid.

**Combined topological sweep (`_network_order`).** Routing visits a *combined*
graph of trap nodes and path nodes (edges: trap→spill-path, path→target-trap,
tributary-path→merged-into-path) and processes them in topological order. The
trap list alone is ordered upstream→downstream, but that does **not** fix the
trap/path interleaving once tributaries exist (two traps feeding one path are
unordered relative to each other), so the combined graph is sorted explicitly.

### 3.3 Rate function (`dynNetworkRateFunction!`)

The in-place ODE right-hand side, analogous to NBS.jl's
`NBSNetworkRateFunction!`. For each trap, using the unified expression
`dV/dt = inflow − loss − drainage − spill` (drainage = 0 in the current scope):

- **Full trap** (`V ≥ capacity`): `loss = footprint_infiltration`,
  `spill = max(inflow − loss, 0)`. While adequately fed this gives `dV ≈ 0`
  (a steady pass-through); once inflow drops below losses it goes negative
  (onset of draining).
- **Accumulating trap** (`V < capacity`): `loss = wetted_infiltration(V)`,
  `spill = 0`, so `dV = inflow − wetted_infiltration(V)`.

Its parameter bundle `DynNetworkRateParams` (built by `_build_rate_params`) holds
everything static for a solve: geometry, constant inflow, the infiltration sums,
and the routing plan.

### 3.4 Event-detection callback

`_build_event_callback` constructs a `VectorContinuousCallback` whose condition
vector generalises `fill_trap_until`'s triple to every **evolving** trap:

| Condition | Expression | Fires when |
|---|---|---|
| `:fill` | `capacity − V` | trap fills, starts spilling |
| `:empty` | `V` | trap empties, exposes subtraps |
| `:stagnation` | `dv0 · dV` | net rate changes sign → steady state |

`affect!` records which condition fired (handling DiffEq's scalar/vector index
forms) and terminates. **Only evolving traps are monitored** — a full trap
already sits at `capacity − V = 0` and would trigger `:fill` spuriously.

Dormant scaffolds (`:unspill`, culvert-outlet thresholds) plug into the same
vector but generate no entries until a state-dependent drain (culvert/NBS) exists:
`:unspill` is emitted only for traps with culvert outlets, of which there are none
in the current scope.

### 3.5 Driver (`solveDynNetwork`)

1. Build the static `DynNetworkRateParams`.
2. Identify evolving traps (`V < capacity`).
3. Short-circuit the steady-state cases: no evolving traps, or all evolving rates
   already ≈ 0 → return `(Inf, 0, :none, state)`.
4. Integrate `ODEProblem(dynNetworkRateFunction!, V0, (0.0, Inf), p)` with the
   first-event callback; the solver auto-selects the algorithm (as
   `fill_trap_until` does).
5. Map the result: `:stagnation`/no-trigger → `(Inf, 0, :none)`; otherwise
   `(time, global trap index, kind, state)`.

---

## 4. Key design choices and insights

These are the non-obvious points needed to understand *why* the code is shaped
the way it is.

### 4.1 Constant inflow ⇒ steady state means "no further event"

Inflow is a constant rate per trap. With only lakes and rivers, the only state
that evolves is accumulating (leaf) traps; full traps are steady pass-throughs.
A leaf's net fill rate is monotonically non-increasing as it fills (more submerged
cells → more infiltration), so a solve has exactly two outcomes: the leaf fills
(finite-time event) or its rate decays to zero (steady state). Hence **no `tspan`
argument** — steady state is reported as `time = Inf`.

### 4.2 Wetted-area infiltration is a *step* function — so stagnation is reached in finite time

This is the crucial insight that makes the event detection work. `water_level(V)`
is continuous, but `wetted_infiltration` counts cells with `bottom ≤ level`, which
**jumps** each time the level crosses a cell elevation. So the net rate `dV` is
**piecewise constant** in `V`. A filling leaf therefore does not approach its
equilibrium asymptotically; it crosses strata at constant rate and either fills,
or hits a stratum where the rate flips negative — a stable equilibrium pinned
*at* a cell elevation, reached in **finite time**. This is why the sign-change
stagnation detector (`dv0 · dV`, borrowed from `fill_trap_until`) fires at a
finite time rather than chasing an asymptote, and why an unbounded `(0, Inf)`
integration horizon with a terminating callback is safe.

*(Verified empirically: in the controlled single-trap case, the solve stops
exactly at the threshold volume with the equilibrium bracketed —
`wetted_below < inflow < wetted_above`.)*

### 4.3 Leaf-vs-full is emergent, not hardcoded

A trap is treated as full when `V ≥ capacity`; otherwise it accumulates. Nothing
labels traps as "leaf" topologically. This is deliberately forward-compatible:
once culverts/NBS add state-dependent drainage, a non-leaf trap can sit below
capacity (draining downstream without spilling), and the same unified rate
expression already covers it.

### 4.4 Inflow is additive; paths transport but never originate

A trap's total inflow is the caller's constant inflow **plus** routed-in flow from
upstream. Flow paths only transport (and lose water to infiltration); they carry
no source of their own.

### 4.5 Merge approximation

When paths overlap, the later one is truncated and registered as a tributary; the
junction cell is not retained. So a tributary's delivered flow is added at the
**head** of the path it merges into, charging it that path's *full* infiltration
capacity rather than only the post-junction portion. This is **exact without
merges** and a slight over-infiltration with them — documented in the
`_route_flow` docstring as a point to revisit if it ever matters.

### 4.6 Performance posture

The deliberate decision (recorded in `prompt_refined.org`) is to **keep the code
simple and rebuild per solve**, revisiting only if a real bottleneck is measured —
because networks are built and solved frequently. Two consequences:

- **Done now (genuine waste, not premature):** the routing plan (`_network_order`
  + `_merge_targets`) is precomputed once into `DynNetworkRateParams` rather than
  rebuilding the graph and topological sort on every ODE step.
- **Deferred (noted in code/org):** `_build_trap_geometry` calls
  `_compute_z_vol_tables` over *all* traps though only the network's are used;
  `v2z` is stored untyped; the rate function allocates small per-call buffers.
  The intended future optimisation is to cache the structure-invariant per-trap
  geometry once for the whole simulation and have each solve select entries.

---

## 5. Testing

### 5.1 Durable tests (`test/dynamics_test.jl`)

Two testsets were added; both pass in full (`7/7` and `363/363`).

- **`networksolver: flow routing (pure topology)`** — hand-built networks (no
  `TrapStructure` needed), with results checked against values computed by hand:
  no-loss chain, path loss, footprint loss, leaf-doesn't-spill, loss floored at
  zero, **tributary merge** (`[10, 27, 20]`), and spill-exits-domain.

- **`networksolver: geometry / rate / solve on mini.txt`** — integration tests
  against a real `TrapStructure` (the `mini.txt` artifact, 58-trap network):
  - *Geometry:* for every trap, empty level = bottom; full level = `Inf` and
    whole-footprint infiltration; half-full level within the trap and infiltration
    monotone in fill; and `water_level` cross-checked against the independent
    Roots-based `waterheight`.
  - *Rate function:* zero-infiltration case (full traps net to ≈0, leaf
    accumulates) and a starved case (every full trap drains at exactly
    `−footprint_infiltration`).
  - *Driver:* a **fill** event (correct global trap index, finite time, leaf at
    capacity), a **stagnation → `Inf`** case, and the **no-evolving** case (all
    traps full → immediate steady state).

### 5.2 Additional verification performed during development

Beyond the durable suite, each component was validated in isolation as it was
built — including the controlled single-trap **fill vs stagnation** scenarios (the
stagnation equilibrium bracketing in §4.2) and the network-level callback run.
These confirmed the step-function/finite-time stagnation behaviour and the
`(0, Inf)` integration terminating correctly via the callback.

### 5.3 Suite status

The full package test suite passes **except for two pre-existing errors in
`MaskingTests`**, which are unrelated to this work: a stray `t` at
`src/utils.jl:320` (a `# constant` comment split across lines 319–320) makes
`flatten_grid!` reference an undefined variable. This file was not touched by the
network-solver work and the failure predates it; it is left for separate
attention.

---

## 6. Current limitations / future work

- **Dormant events.** `:empty` (expose-subtraps), `:unspill` (stop-spilling), and
  culvert-outlet threshold crossings are scaffolded into the event vector but
  cannot fire in the current lakes-and-rivers, constant-inflow scope; they become
  live once state-dependent drainage (culverts/NBS) is added.
- **Culverts / NBS.** `DynCulvert`s and NBS systems are not yet modelled as
  evolving network components; the rate expression and event taxonomy were shaped
  to accommodate them with minimal change.
- **Merge fidelity** (§4.5) and the **performance caching** (§4.6) are the two
  documented places to revisit if accuracy or speed ever demands it.
- **Units** are gridcells × height (1 m assumed for both); the infiltration grid
  is static between solves but kept as a per-solve input to leave room for a
  future dynamic-saturation extension.
