# Network solver — call graph

The dynamic network solver lives in `src/dynamics/` and splits into two phases:
**construction** (`elements.jl` — build the `DynNetwork` topology) and **solving**
(`networksolver.jl` — integrate water volumes forward to the next event).
`culvert_rate.jl` is a leaf used by both.

---

## Phase 1 — Network construction (`elements.jl`)

```
setup_network(tstruct, dyn_coords, full_traps; culverts)        elements.jl:181
├─ _subnetwork(tstruct, coord, full_traps)                      :235   [per start coord]
│  ├─ flow_path_from(tstruct, ...)                              watercourses.jl:142
│  ├─ _unfilled_trap_at(tstruct, cell, full_traps)             :263
│  └─ _build_network(paths, traps, starts_with_trap, tstruct)  :275
│     └─ _merge_networks(networks)                             :301 (1-arg → 3-arg)
├─ _expand_with_culverts(tstruct, subnets, culverts, full_traps):210   [only if culverts]
│  ├─ _occupied_cells(tstruct, subnets)                        :191
│  └─ _subnetwork(...)                                          (adds inlet/outlet subnets)
└─ _merge_networks(networks, cv_objs, tstruct)                 :301
   ├─ _combine_subnets(subnets)                                :384
   ├─ _dedup_traps(all_paths, all_traps)                       :328
   ├─ _resolve_cell_overlaps!(all_paths)                       :417
   ├─ _culvert_owners(tstruct, all_paths, all_traps, cv_objs)  :359
   ├─ _components(all_paths, all_traps, inlet_owner, outlet_owner) :454
   │  └─ unite!(x, y)                                          :459 (union-find, nested)
   └─ _build_component(...)                                    :559   [per connected component]
      └─ _topological_order(...)                               :498
         └─ Graphs.topological_sort_by_dfs / is_cyclic
```

Output: a `Vector{DynNetwork}`, each a connected component of `DynFlowPath` +
`DynTrap` (+ `DynCulvert`) nodes.

---

## Phase 2 — Solving (`networksolver.jl`)

`solveDynNetwork` is the driver; it advances to the *first topology-changing
event*, so a full simulation calls it in a loop (see
`examples/verification/solve_dynamic_network.jl`).

```
solveDynNetwork(tstruct, net, infiltration, inflow, state; ...)  :813
│
├─ _build_rate_params(tstruct, net, infiltration, inflow; ...)   :540
│  ├─ _network_order(net)                                        :189
│  │  ├─ _add_culvert_edges!(g, net, np)                         :257
│  │  └─ Graphs.topological_sort_by_dfs
│  ├─ _path_cell_infiltration(net, infiltration)                 :464
│  ├─ _infil_prefix(cell_infil)                       [per path] :471
│  ├─ _merge_targets(net)                                        :213
│  ├─ _build_culvert_plan(net, tstruct)         [if culverts]    :239
│  ├─ _build_trap_geometry(tstruct, net, infiltration; zvt)      :72
│  │  ├─ _compute_z_vol_tables(tstruct)         [if zvt===nothing] (trapvolumes.jl)
│  │  ├─ subtrapsof(tstruct, tix)                                (sshierarchy.jl)
│  │  └─ Interpolations.linear_interpolation
│  └─ _footprint_infiltration(tstruct, net, infiltration)        :450
│
├─ dynNetworkRateFunction!(du0, V0, p)        ── initial rate    :608
│  └─ _routed_inflow(V, p)                                       :592
│     ├─ _surface_level(geom[i], V[i])         [if culverts]     :136
│     │  └─ g.v2z(...)                          (interpolant)
│     └─ _route_flow(net, ...)  ★ CORE ROUTER ★                 :343  (9+ arg form)
│        ├─ _culvert_flow(cvplan, net, ci, trap_level)  [if culverts] :279
│        │  └─ culvert_rate(cv, tstruct; ...)                    culvert_rate.jl:86
│        │     └─ _directional_capacity(...)                     culvert_rate.jl:27
│        └─ (wetted_infiltration / water_level used by caller)
│
│   [the convenience _route_flow at :322 wraps the core one — only the
│    core 9-arg form is on the hot path, via _routed_inflow]
│
├─ _build_event_callback(p, evolving, V0)                        :712
│  ├─ dynNetworkRateFunction!(dv0, V0, p)   ── baseline for sign test
│  ├─ _event_conditions(p, evolving)                             :681
│  ├─ condition(out, V, t, integrator)      ── nested, per ODE step :724
│  │  ├─ _routed_inflow(V, p)                  (→ _route_flow ...)
│  │  └─ wetted_infiltration(geom[i], V[i])                      :153
│  │     └─ water_level(g, V)                                    :127
│  └─ affect!(integrator, ix)               ── nested; sets event, terminate! :749
│
└─ solve(ODEProblem(dynNetworkRateFunction!, V0, (0.0,Inf), p); callback=cb)
   └─ [DifferentialEquations.jl repeatedly calls]
      ├─ dynNetworkRateFunction!  →  _routed_inflow  →  _route_flow
      └─ condition / affect!  (the callback above)
```

---

## How to read this for understanding

- **`_route_flow` (`:343`) is the heart.** Everything funnels through it. It does
  one topological pass over the merged path/trap/culvert node order, propagating
  inflow downstream while charging infiltration per-segment. Read its
  segmented-routing loop (`:367`–`:438`) carefully — that's where mass
  conservation lives.
- **`_routed_inflow` (`:592`) is the shared hot path**, called both by the ODE
  rate function and by the event `condition`. So the router runs many times per
  solve.
- **The two `_route_flow` methods** at `:322` and `:343`: the first is a
  convenience wrapper that computes the static routing data on the fly (used in
  tests); the solver always uses the 9-arg core form via precomputed
  `DynNetworkRateParams`.
- **Three nested closures** inside `_build_event_callback` (`condition`,
  `affect!`) are what DiffEq drives — they're not top-level functions, easy to
  miss when grepping.
- **Leaf geometry helpers** — `water_level`, `_surface_level`,
  `wetted_infiltration` (`:127`–`:162`) — are the volume↔level↔infiltration
  conversions every step relies on.
