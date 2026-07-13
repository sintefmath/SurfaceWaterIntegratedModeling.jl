# `networksolver.jl` call graph + `fill_sequence` hookup

ASCII call graph of the dynamic network solver and how it plugs into `fill_sequence`.
`[NBS]` = NBS-overlay only (Phase-B `DynNBS`→`DynNBSPlacement` targets).
`[cv]` = culvert-only. `[OLD]` = the current (unreviewed) integration, replaced by the gate.

## 1. Hookup: fill_sequence → network_context → solver

```
fill_sequence (main event loop, fill_sequence/fill_sequence.jl)          [OLD integration]
│
├─ _build_dyn_networks              (network_context.jl L196)   period start: build contexts
│    ├─ _dyn_seeds                                              seeds = dyn_traps + culvert ends
│    ├─ _nbs_elements               [NBS]  makes DynNBS from placements   ← B3 deletes
│    ├─ setup_network (elements.jl, OLD 3-arg)                  ← D retires (new build_network is the replacement)
│    ├─ _make_context               builds DynNetworkContext (state, global_ix, extern_inflow, …)
│    │    ├─ _nbs_layer_block       [NBS]
│    │    ├─ _external_inflow       leaf-sum of trap_inflow
│    │    └─ _inflow_sources
│    └─ _predict_network!  ─────────────────────────────┐  next-event prediction on a COPY
│                                                        │
├─ _expand_empty_fill_updates       (network_context.jl L530)   :empty parent ⇒ expose children
│                                                        │       (the S2 hazard guard)
├─ _touch_networks!                 (network_context.jl L451)   per-event structural update
│    ├─ _affected_contexts / _capture_fired_kinds        │
│    ├─ setup_network_cached (elements.jl, SubnetCache)   │      ← C/D replace with apply_fill!/
│    │    └─ _subnetwork / _merge_networks (+ [NBS])      │        apply_unfill!/apply_empty!
│    ├─ _reconcile_spillgraph!       full recompute+diff  │      ← C keeps a minimal version
│    ├─ _reuse_plan / _contexts_to_commit                 │
│    ├─ _commit_contexts!  ──────────────────────────────┤
│    │    ├─ _commit_network!  ──────────────────────────┤
│    │    ├─ _clamp_full_traps!                           │      state-handoff (reusable, §S1)
│    │    └─ _store_nbs_state!        [NBS]               │
│    ├─ _apply_fired_boundaries!      :empty/:unspill     │      state-handoff (reusable)
│    └─ _assemble_contexts                                │
│         ├─ _make_context (rebuilt) / _predict_network! ─┤
│         └─ fill_trap_until (project newly-absorbed)     │      stale-state hazard (§S1.2)
│                                                         │
├─ _network_amount_updates          (network_context.jl L554)   node vol / subsumed=cap → cur_amounts
│                                                         │
└─ _finalize_networks!              (network_context.jl L505)   period end: settle to endtime
     ├─ _commit_network!  ──────────────────────────────┤
     └─ _store_nbs_state!            [NBS]               │
                                                         │
   _predict_network! / _commit_network! ─────────────────┘
        └─ solveDynNetwork!   ◀── the one entry point into networksolver.jl
```

The whole left column (`[OLD]`) is the integration layer the gate rewrites (Phases C/D). The
**only** contract with the solver is `solveDynNetwork!(state, tstruct, net, infiltration,
inflow; tmax, nbs_placements, nbs_inflow, …) -> (time, trap, kind)`.

## 2. `solveDynNetwork!` internals (networksolver.jl)

```
solveDynNetwork!                                                   (L1492)
├─ _build_rate_params                                             (L860)  static per-solve inputs
│   ├─ _network_order            path/trap topo order (+ _add_culvert_edges! [cv])
│   ├─ _path_cell_infiltration → _infil_prefix
│   ├─ _build_culvert_plan       [cv]
│   ├─ _build_nbs_plan           [NBS]  ← B2 REWORK (reads nb.placement_ix today; must use
│   │    └─ _nbs_exit_weights    [NBS]        nb.system / nb.id from DynNBSPlacement)
│   ├─ _path_event_templates     [cv]/[NBS]
│   ├─ _build_trap_geometry  → _compute_z_vol_tables (fill_sequence)
│   └─ _footprint_infiltration
├─ _validate_network                                              (L1369)  three-state contract
├─ _reconcile_submergence!       [NBS]                            (L1166)
├─ dynNetworkRateFunction!       (t=0 fast-path probes: :unspill / :empty)
├─ _build_event_callback         → VectorContinuousCallback (LeftRootFind)   (L1297)
│   ├─ _event_conditions         (:fill / :empty / :unspill / :submerge[NBS])(L1135)
│   ├─ condition(out,V,t): _routed_inflow  (per step)
│   └─ affect!: _apply_submergence! [NBS] | else set event + terminate!
├─ _build_steadystate_callback   (+ _build_steadystate_callback_nbs [NBS])
├─ solve(ODEProblem, …)          DifferentialEquations
└─ returns (; time, trap=global_ix, kind ∈ :fill|:empty|:unspill|:none)
```

## 3. Rate function / routing tree (hot path, per ODE step)

```
dynNetworkRateFunction!(dV, V, p)                                 (L987)
├─ _routed_inflow(V, p)                                           (L928)
│   ├─ (per full trap) spilling[i] = V[i] >= capacity
│   ├─ _surface_level              [cv] trap heads
│   ├─ _nbs_fill_actual!           [NBS] layer overflow → delivery slots
│   ├─ _route_flow (core, L442)    topological routing over paths+traps
│   │    ├─ _path_delivered!       segmented infil + tribs (+ events [cv]/[NBS])
│   │    │    └─ _culvert_flow      [cv]  → culvert_rate
│   │    └─ _route_trap_node!      apply culverts, spill surplus into spill_path
│   │         └─ _culvert_flow      [cv]
│   └─ _nbs_captured / _nbs_saturated_draw   [NBS] submerged exchange
├─ (per trap) dV = inflow − loss − spill    (wetted_infiltration for accumulating)
└─ (per NBS layer) [NBS] compute_outflow cascade (dry / submerged)
```

## 4. What Phase B touches

- **B2 (solver, self-contained):** `_build_nbs_plan` (+ `NBSPlan.placement_ix`,
  `_build_rate_params`'s `placements`/`nbs_inflow` args, `_nbs_captured`,
  `_reconcile_submergence!`, `_apply_submergence!`) — switch from `DynNBS.placement_ix`
  indirection to `DynNBSPlacement.system` / `.id` read straight off `net.nbs`. Everything
  else in §2/§3 is type-agnostic.
- **B3 (delete `DynNBS`):** the struct (`elements.jl`) + `_nbs_elements`
  (`network_context.jl`) + the OLD `elements.jl` NBS path (`_expand_with_nbs`,
  `_nbs_union_edges`, `_build_component` `comp_nbs`, `_merge_networks` NBS args) +
  `nbs_overlay_test` / `nbs_dynamic_test`. Same dead code Phase D removes.
```
