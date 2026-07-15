# Plan: NBS as a routing node emitting a signed (O_1 − O_0) correction

## Context

The current NBS handling is an oblivious-correction model whose correction is
computed from the *static* oblivious runoff, so it (a) cannot abstract
network-sourced inflow crossing the footprint, (b) does not cascade past the
first trap (applied after `_route_flow`), and (c) carries a pile of bespoke
machinery (`NBSLayerParams`, `NBSDelivery`, `_walk_to_trap`,
`_propagate_correction`, precomputed `capture`).

**Constraint (settled):** `watercourses` must stay **NBS-oblivious** — NBS
elements cannot own watersheds; they intersect existing ones. So the oblivious
overland flow still passes through the footprint and is already routed downstream.

**Design:** make the NBS a **node in the network `_route_flow` traverses**, and
have it emit only the *delta* from what oblivious routing already delivered.

- `O_0` — **oblivious (background) output**: sum of `runoff` at the footprint's
  drainage endpoints — terrain exit cells (`footprint_outflow_cells`) and internal
  depressions (`internal_accumulation_cells`) — read from the grid. The grid is
  oblivious to the network (spillgraph = `FULL ∧ ¬COVERED`), so this is the
  **background** pass-through only, exactly what oblivious routing delivered
  downstream. Static, precomputed. (Piped outlets have no `O_0`.)
- `I_1` — **input**: background (`Σ runoff[input_cell]`, grid) **+** router
  deliveries to the NBS (network-trap spills + upstream outlet discharge) **+**
  precip. Grid and router are complementary (grid excludes the network), so no
  double-count. Drives the layer ODE. Measured at the input cells (not output-side)
  — avoids double-charging the footprint's own infiltration.
- `O_1` — **real output**: the layers' `O_terrain + E` at their current state.
  Dynamic.
- **Correction `O_1 − O_0`** issued at the output cells. Downstream ends at
  `O_0 + (O_1 − O_0) = O_1`. No double-count; `watercourses` untouched.

The correction can be **negative**, so the flow tracking gains a signed component.

## Design

### Input `I_1` — background from the grid, network/outlet from the router
Crucial fact: `compute_flow` (flow.jl) tracks runoff only for traps **in the
spillgraph**, and the spillgraph is `FULL ∧ ¬COVERED` (`_reconcile_spillgraph!`
withdraws network traps). So the `runoff` grid is genuinely **oblivious to the
network** — it contains background + non-network spills, *not* network-trap spills.
The grid and `_route_flow` are therefore **complementary** (no overlap), and can be
summed:

```
I_1 = Σ runoff[input_cell]        # background (grid, non-network)
    + router deliveries to the NBS # network spills + upstream outlet discharge
    + precip_on_footprint
```

`I_1` drives the layer ODE. The router deliveries reach the NBS because the topo
sort routes upstream spills/outlets along their `DynFlowPath`s and they **terminate
at the footprint** (see build_network change) — so no consuming/zeroing and no
tracking flow inside the footprint. No double-count (the grid never carried the
network part). Measuring at the input cells (not output-side `Q_c`) also avoids
double-charging the footprint's own infiltration.

### Output correction — emitted at the output cells, cascading
The output stays fully dynamic: `O_1 = O_terrain(V) + E(V)` from the evolving layer
state. The NBS emits the delta from what oblivious routing already sent (`O_0`),
via `:nbsout` events on the output paths (reuse the culvert `:cvout` machinery) into
the path's **signed** component:
- **terrain**: `(O_terrain(V) − O_0)·ratio_e`, `ratio_e = O_0[e]/O_0_total`
  (downstream total `O_terrain·ratio_e`). Guard `O_0_total ≈ 0` → even split
  (mirrors the current `O_terrain/nend` fallback). Internal-depression endpoints
  deposit their share straight into the accumulation trap (no path).
- **piped**: `+E[l](V)` at each outlet `l`, direct.
Totals to `O_1 − O_0`. The terrain split math equals today's `−X·Q_c`
(`Q_c = O_0[e]`); the change is the input-side `I_1` and the cascade.

### Signed flow tracking (the accepted addition)
- `_path_delivered!` carries, alongside the positive `current`, a **signed**
  correction component attenuated **V-relative**: `c_out = max(V+c,0) − max(V,0)`,
  `V` = oblivious runoff at the cell. This is the correct rule for a correction to
  the oblivious baseline (which already paid each cell's infiltration) — NOT the
  positive `max(current − infil, 0)` clamp. It's exactly what `_track_flow!` already
  does (its negative branch); `_route_flow` gains the same.
- The correction is delivered to the path's target trap, so a full trap spills the
  corrected amount and it **cascades through full traps** — closing the earlier
  "doesn't cascade past the first trap" gap.

### Structs
- Type `NBSLayer` fields `Float64` → delete `NBSLayerParams`
  (`NBSSystem.layers::Vector{NBSLayer}` already concrete).
- Delete `NBSDelivery`, `_walk_to_trap`, and the standalone `_apply_nbs_corrections!`
  (its work moves into the router's NBS-node case). `_propagate_correction`'s rule
  becomes the router's signed-attenuation step.
- `NBSElement`/`NBSPlan` slim to the per-placement plan: `system` ref,
  `state_base`, `n_terrain`, precomputed **background** `I_bg = Σ runoff[input_cell]`
  and `O_0` (with per-endpoint `ratio_e`), and the output routes. The dynamic router
  deliveries are added to `I_bg` at solve time (see build_network).
- Reuse `DynNBSPlacement.footprint_inflow_cells` / `footprint_outflow_cells` /
  `internal_accumulation_cells` / `outlets`; `DynFlowPath.nbs_outlets` and the
  `_add_accumulation_traps!` region→trap map (network_utils.jl) for the output
  routes and internal-depression deposits.

### Build-network change (NBS as a routing node)
- Paths must **terminate at NBS footprints** (target = the NBS), the same way
  build_network already strips trap-footprint cells (build_network.jl:133), so an
  upstream spill/outlet path delivers its flow to the NBS instead of running through
  it. `DynFlowPath.target_trap` gains an NBS-target form (or a parallel field), and
  `_network_order` includes NBS nodes so the topo sort processes each NBS after its
  delivering paths. The NBS node's received deliveries + `I_bg` + precip = `I_1`.

## Open implementation questions (small)
- Exact threading of the signed correction component through `_path_delivered!`
  (second accumulator vs a tagged pass) and the trap cascade.
- Placing `:nbsout` at output-cell positions on paths (reuse the `nbs_outlets` /
  `_intersecting_on_path` bookkeeping already added for culverts).

## Critical files
- `src/dynamics/networksolver.jl` — `_network_order` (NBS nodes), `_route_flow`
  (NBS-node case: receive deliveries → `I_1`, emit `O_1 − O_0`; signed component in
  `_path_delivered!`), the NBS plan, the rate-function NBS block; delete the old
  correction machinery.
- `src/dynamics/build_network.jl` — paths terminate at NBS footprints (NBS-target);
  `src/dynamics/elements.jl` — `DynFlowPath` NBS-target form.
- `src/dynamics/nbs_elements.jl` — type `NBSLayer`.
- Reuse: `DynNBSPlacement` cell lists (elements.jl); network_utils mappings;
  `DynFlowPath.nbs_outlets`. `watercourses.jl` unchanged (stays oblivious — the
  spillgraph excludes network traps, so the grid is complementary to the router).

## Verification
- `test/nbs_routing_test.jl` — retention/pass-through must still hold
  qualitatively; values will change (dynamic `I_1` vs oblivious). Verify physically
  and update pinned numbers with justification.
- New **cascade** test: NBS feeding a full trap that spills into a second trap —
  the delta must reach the second trap.
- New **network-inflow** test: NBS fed mainly by an upstream network spill (small
  background) — it must abstract that flow (old model would not).
- New **upstream-outlet** test: a culvert/NBS outlet upstream whose discharge reaches
  the footprint via a `DynFlowPath` — the NBS must capture it (grid-only read would
  miss it; the router delivery must reach the NBS node).
- Signed-tracking unit test: the V-relative rule on a hand path (mirrors the
  existing `_propagate_correction` test).
- Mass conservation: net surface water removed over a step = layer storage change
  (`−dS`); assert it.
- Run `dynamics_test.jl`, `dynamic_membership_test.jl`, `network_driver_test.jl`,
  the `Sequencing` set; package loads clean.
