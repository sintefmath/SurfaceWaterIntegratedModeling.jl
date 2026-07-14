# `networksolver.jl` — call graph

The multi-trap ODE solver: evolve a `DynNetwork`'s trap volumes forward to the
**first topology-changing event** (`:fill` / `:empty` / `:unspill`), to steady
state, or to a `tmax` cutoff. Network *construction* is separate (`elements.jl` /
`build_network.jl`, already reviewed); this file only *solves* a built network.

Line numbers are as of the current `src/dynamics/networksolver.jl`.

**Exports:** `solveDynNetwork!`, `dynNetworkRateFunction!`.
**External callers:** `_predict_network!` / `_commit_network!` (`network_context.jl`)
call `solveDynNetwork!` once per predict/commit.

---

## Big picture

```
                       _predict_network! / _commit_network!   (network_context.jl)
                                        │
                                        ▼
                              solveDynNetwork!                 :1085   ── the driver
                            ╱        │         ╲
              _build_rate_params  dynNetworkRateFunction!  callbacks
              (static, once)      (ODE RHS, every step)    (event + steady-state)
                                        │
                                        ▼
                                 _routed_inflow  ── shared hot path ──┐
                                        │                             │
                                   _route_flow  ★ core router ★   _apply_nbs_corrections!
```

Three things run per ODE step: the **rate function**, the **event condition**,
and the **steady-state condition** — and all three funnel through
`_routed_inflow` → `_route_flow`. That router is the hot path.

---

## 1. Build phase — `_build_rate_params` (static, once per solve)

Produces the immutable `DynNetworkRateParams` (`:654`) the rate function reads.

```
_build_rate_params(tstruct, net, infiltration, inflow; path_inflow, runoff, zvt)  :678
├─ _network_order(net)                                             :155
│  ├─ _add_culvert_edges!(g, net, np)                              :217
│  └─ Graphs.topological_sort_by_dfs
├─ _path_cell_infiltration(net, infiltration)                     :433
│  └─ _infil_prefix(cell_infil)                        [per path]  :440
├─ _merge_targets(net)                                            :177
├─ _build_culvert_plan(net, tstruct)                 [if culverts] :199
├─ _build_nbs_plan(net, tstruct, runoff)             [if NBS]      :597
│  ├─ _nbs_endpoints(nb, tstruct, runoff, trap_of_cell)           :561
│  │  └─ _walk_to_trap(start, tstruct, runoff, trap_of_cell)      :542
│  └─ _nbs_outlets(nb, tstruct, runoff, trap_of_cell, LI)         :582
│     └─ _walk_to_trap(...)                                       :542
├─ _path_event_templates(net)                        [if culverts] :238
├─ _build_trap_geometry(tstruct, net, infiltration; zvt)          :54
│  ├─ _compute_z_vol_tables(tstruct)     [if zvt===nothing]  (fill_sequence.jl)
│  ├─ subtrapsof / numregions                                (TrapStructure / sshierarchy)
│  └─ Interpolations.linear_interpolation           ── builds v2z
└─ _footprint_infiltration(tgeom)                                 :422
```

`runoff` is read **only** to build the NBS plan; omit it for NBS-free nets. A net
with no culverts gets `cvplan = nothing` and no `path_events`; a net with no NBS
gets `nbsplan = nothing`.

---

## 2. Rate function + routing (hot path, every ODE step)

```
dynNetworkRateFunction!(dV, V, p, t)                              :751
├─ _routed_inflow(V, p)                                           :735
│  ├─ _surface_level(geom[i], V[i])                  [if culverts] :118   (→ v2z)
│  ├─ _route_flow(net, ..., cvplan, trap_level, path_events)  ★   :367   ── CORE ROUTER
│  │  │   one topological pass over path+trap nodes (order from _network_order)
│  │  ├─ _path_delivered!(prefix, head, tribs, events, ...)      :311   [path node]
│  │  │  └─ _culvert_flow(cvplan, net, ci, trap_level)  [culvert] :254
│  │  │     └─ culvert_rate(cv, ...)                         (culvert_rate.jl)
│  │  └─ _route_trap_node!(i, net, trap_inflow, ...)             :343   [trap node]
│  │     └─ _culvert_flow(...)                          [culvert] :254
│  └─ _apply_nbs_corrections!(inflow, V, p, nt)      [if NBS]     :626
│     ├─ compute_outflow(...)                    (nbs_elements.jl)  ── O_terrain, E_piped
│     └─ _propagate_correction(c, path_V)                        :532   ── max(V+c,0)-max(V,0)
│
├─ per-trap dV:  full → inflow - footprint_infil - spill
│                accumulating → inflow - wetted_infiltration(geom, V)   :130 (→ water_level :109)
└─ per-NBS-layer dV  [if NBS]:  compute_outflow (overflow qo + infil qi)  (nbs_elements.jl)
```

Two `_route_flow` methods: the 6-arg convenience wrapper (`:291`) computes the
static routing data on the fly (tests / hand-built nets); the solver always hits
the pre-supplied core form (`:367`) via `DynNetworkRateParams`.

Mass conservation lives in `_path_delivered!`'s segmented loop and
`_route_trap_node!`'s culvert deliver/drain — read those two together.

---

## 3. Event detection (halts the integration)

Two callbacks combined in a `CallbackSet` (a trap-free NBS-only net uses only the
steady-state one):

```
_build_event_callback(p, evolving, V0, nreg)                     :926    ── cb_topo (VCB, LeftRootFind)
├─ _event_conditions(p, evolving, nreg)                          :841    ── :fill/:empty/:unspill list
├─ condition(out, V, t, integrator)   [nested, per step]
│  └─ _routed_inflow(V, p)                                       :735    (:unspill needs routed inflow)
└─ affect!(integrator, ix)            [nested]  ── sets event.{kind,trap}, terminate!

_build_steadystate_callback(p, evolving, abstol, du0)            :863    ── cb_ss, NBS-free
│  └─ condition:  _routed_inflow + wetted_infiltration, veto if any evolving trap spilling
_build_steadystate_callback_nbs(p, ss_indices, abstol, du0)     :895    ── cb_ss, with NBS
   └─ condition:  full dynNetworkRateFunction! into a reused dbuf, settle over traps + layers
```

`cb_ss` is a `DiscreteCallback` (fires only at accepted steps) *by design* —
`wetted_infiltration` is a step function of `V`, so a continuous callback would
detect spurious crossings mid-step.

---

## 4. Driver sequence — `solveDynNetwork!` (:1085)

```
solveDynNetwork!(state, tstruct, net, infiltration, inflow; tmax, path_inflow, runoff, ...)
│
├─ _build_rate_params(...)                              → p        (section 1)
├─ _validate_network(tstruct, net, V0, p.geom)          :964       ── three-state contract
├─ dynNetworkRateFunction!(du0, V0, p, 0.0)                        ── initial rates
│
├─ _t0_fast_path(V0, du0, p, nreg)                      :1010      ── immediate-event checks:
│  ├─ FULL trap with du0<0            → (:unspill, t=0)
│  └─ parent-EMPTY with neg net rate  → (:empty,   t=0)   [uses _routed_inflow]
├─ further t=0 returns (inline):
│  ├─ tmax<=0                          → (:none,    t=tmax)
│  └─ nothing evolving / all rates≈0  → (:none,    t=Inf)  [steady state]
│
├─ evolving = traps with V != capacity;  ss_indices = evolving (+ NBS layer states)
├─ cb_ss = _build_steadystate_callback[_nbs](...)                 (section 3)
├─ cb_topo, event = _build_event_callback(...)          [if nt>0] (section 3)
│
├─ solve(ODEProblem(dynNetworkRateFunction!, V0, (0,tmax_ode), p), Tsit5();
│        callback = CallbackSet(cb_topo, cb_ss))
│     └─ DiffEq repeatedly calls the rate function + the two callback conditions
│
├─ state .= sol.u[end]                                            ── write back in place
└─ classify result:
   ├─ event.kind == :none →  tmax cutoff (t≈tmax_ode) vs steady state (t=Inf)
   └─ :fill/:empty/:unspill →  clamp triggering trap to C or 0, return global trap_ix
```

Explicit `Tsit5` (not the auto-switcher): the RHS has C0 kinks (culvert regime
switches, routing `min`) but isn't stiff — the auto-switcher misreads the kinks
as stiffness and pays for dense Jacobians (~4× slower). See the comment at the
`solve` call.

---

## Leaf / external dependencies

| Symbol | Where | Used for |
|--------|-------|----------|
| `culvert_rate` | `culvert_rate.jl` | HDS-5 barrel capacity, via `_culvert_flow` |
| `compute_outflow` | `nbs_elements.jl` | power-law layer overflow/infiltration (NBS) |
| `_compute_z_vol_tables` | `fill_sequence.jl` | volume↔level tables for `v2z` |
| `subtrapsof` / `numregions` / `numtraps` | `TrapStructure.jl` / `sshierarchy.jl` | hierarchy + counts |
| `Graphs.*` | Graphs.jl | routing-order topological sort |
| `Interpolations.linear_interpolation` | Interpolations.jl | `v2z` interpolant |
| `solve` / `ODEProblem` / callbacks / `Tsit5` | DifferentialEquations.jl | integration |

---

## Reading guide

- **`_route_flow` core (:367) is the heart.** One topological pass; every step
  routes inflow through it. Mass conservation = `_path_delivered!` (:311) +
  `_route_trap_node!` (:343).
- **`_routed_inflow` (:735) is the shared hot path** — rate function, event
  condition, and (NBS-free) steady-state condition all call it, so the router
  runs many times per solve.
- **NBS is a signed correction, bolted on after routing.** The baseline routing
  is NBS-oblivious; `_apply_nbs_corrections!` (:626) adds `-X·Q_c` per terrain
  endpoint and `+E` per piped outlet, propagated with `_propagate_correction`.
  Layer storage is appended to the ODE state after the `nt` trap volumes.
- **Nested closures** (`condition`/`affect!` inside `_build_event_callback` and
  the two steady-state builders) are what DiffEq drives — invisible to a
  top-level grep.
- **Geometry leaves** — `water_level` (:109), `_surface_level` (:118),
  `wetted_infiltration` (:130) — are the volume↔level↔infiltration conversions
  every step relies on.
- **Result classification** hinges on `event.kind` and the stop time: `:none`
  splits into steady-state (`Inf`) vs `tmax` cutoff; a real event returns the
  global trap index and clamps the trigger to its threshold.
```
