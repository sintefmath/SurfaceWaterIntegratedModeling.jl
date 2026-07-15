# Plan: NBS as a signed-diff router element (culvert-style `:nbsin`/`:nbsout`)

## Context

The current NBS handling is an oblivious-correction model whose correction is
computed from the *static* oblivious runoff, so it (a) cannot abstract
network-sourced inflow crossing the footprint, (b) does not cascade past the
first trap (applied after `_route_flow`), and (c) carries a pile of bespoke
machinery (`NBSLayerParams` [removed], `NBSDelivery`, `_walk_to_trap`,
`_apply_nbs_corrections!`, `_propagate_correction`, precomputed `capture`).

**Constraint (settled):** `watercourses` must stay **NBS-oblivious** — NBS
elements cannot own watersheds; they intersect existing ones. So the oblivious
overland flow still passes through the footprint and is already routed downstream.

**Design (settled):** the NBS is intercepted at **`:nbsin`/`:nbsout` cell events on
the paths crossing it** (approach B — like culverts; paths not terminated,
`target_trap` unchanged). It captures its live total input `I_1` at the input cells,
runs the layer ODE, and emits only the **diff** `O_1 − O_0` from what oblivious
routing already delivered — carried downstream by the **regular flow tracker
extended to signed values**, attenuated against the read-only `runoff` grid. Net
downstream = `O_0 + diff = O_1`; `watercourses` untouched; no double-count.

- `I_1` — live total input: `Σ runoff[input_cell]` (background, grid) **+** the
  `:nbsin` draws (network spills + upstream outlet discharge, from the router) **+**
  precip. Grid and router are complementary (grid excludes the network), no
  double-count. Drives the layer ODE.
- `O_0` — oblivious (background) output: `Σ runoff` at the footprint's drainage
  endpoints (terrain exits + internal depressions). Static, precomputed. (Piped
  outlets have no `O_0`.)
- `O_1 = O_terrain(V) + E(V)` — live layer output. The emitted diff `O_1 − O_0` can
  be **negative**, which is why the flow tracker needs signed support.

See the Design section for the exact `I_1`, diff-split, and signed-tracking rules.

## Design

**Chosen realization = approach (B): culvert-style cell events on paths, NOT path
termination.** The NBS is intercepted at `:nbsin` / `:nbsout` cell events on the
paths that cross it (exactly like `:cvin` / `:cvout`). Paths are **not** terminated
and `target_trap` is **not** changed — no ripple into `build_network`'s trace, the
membership layer, or `network_utils`. Tradeoff (accepted): the network structure
won't *show* the NBS as a node; interception lives in the events. Mitigate with
explicit event kinds + comments.

### The oblivious runoff grid (`rateinfo.runoff`) is the reference
`compute_flow` runs on the `FULL ∧ ¬COVERED` spillgraph, so `runoff` is
**network/NBS-oblivious**: per-cell **net overland flow (positive)** or **remaining
infiltration capacity (negative)**, static within a solve (see `RateInfo` docstring).
It is threaded into `_route_flow` (from `ctx.runoff`) and read **read-only** — never
mutated, no per-step copy.

### Input `I_1` (drives the layer ODE) — captured at `:nbsin`
```
I_1 = Σ runoff[input_cell]         # background (grid, non-network)  — precomputed I_bg
    + Σ :nbsin draws               # network spills + upstream outlet discharge (dynamic)
    + precip_on_footprint
```
- `:nbsin` at each input cell **consumes** the passing router flow: `I_1[nbs] += current;
  current = 0`. Paths pass *through* the footprint carrying 0 afterward (footprint has
  zero infiltration by contract, so the interior cells do nothing). `I_1` is accumulated
  during `_route_flow` and returned for the layer-ODE block.
- Grid and router are **complementary** (grid excludes the network), so summing the
  precomputed `I_bg` and the dynamic draws does **not** double-count.
- Measuring at the input cells (not output-side `Q_c`) avoids double-charging the
  footprint's own infiltration.

### Output diff (fully dynamic) — emitted at `:nbsout`
`O_1 = O_terrain(V) + E(V)` from the live layer state. The NBS emits only the **diff**
from what oblivious routing already delivered (`O_0`):
- **terrain**: `diff_e = (O_terrain(V) − O_0) · ratio_e`, `ratio_e = O_0[e]/O_0_total`.
  Guard `O_0_total ≈ 0` → even split across the terrain endpoints. Internal-depression
  endpoints deposit their share straight into the accumulation trap (no path).
- **piped**: `+E[l](V)` at each outlet `l`, direct.
Totals to `O_1 − O_0`. Terrain split math equals today's `−X·Q_c` (`Q_c = O_0[e]`);
what changed is the dynamic input-side `I_1` and that the diff now rides the router.

### Signed diff tracking = "the regular flow tracker, modified for negatives"
The diff emitted at `:nbsout` is a **signed** value carried by the router along the
path to the target trap. At each cell it is attenuated against that cell's
`runoff` value `V` (read-only), by the single rule
```
diff_out = max(V + diff, 0) − max(V, 0)
```
equivalently, per the sign cases (all four verified against the illustration in
`agent/prompts/flowtracking_adjust.org`):

| diff | cell `V`        | effect                                             |
|------|-----------------|----------------------------------------------------|
| `+`  | `+` (flow)      | passes unchanged (`c`)                             |
| `+`  | `−` (capacity)  | remaining capacity `\|V\|` consumes it until it vanishes |
| `−`  | `+` (flow)      | the flow absorbs it, capped at `−V` (can't remove more than is there) |
| `−`  | `−` (no flow)   | moot → 0                                           |

The attenuated diff is delivered to the path's **target trap**; since the target is
downstream (later in topo order), the trap cascade carries it — closing the
"doesn't cascade past the first trap" gap. This is *not* the router's positive-spill
`max(current − infil, 0)` rule (which charges the infiltration *rate*); it is the
runoff-relative rule above, which requires the per-cell `runoff` reference.

### No double-count (option ii)
`external_inflow[trap]` keeps the oblivious baseline — `O_0`'s contribution as it
reached the trap. The diff-tracker delivers only the **attenuated diff** (the change),
so the trap ends at `O_0 + diff = O_1`. Track the diff, not the total; `runoff` stays a
read-only reference.

### Structs
- Type `NBSLayer` fields `Float64` → delete `NBSLayerParams`
  (`NBSSystem.layers::Vector{NBSLayer}` already concrete). **[done]**
- `NBSElement` references `NBSSystem` + footprint `A`. **[done]**
- Delete `NBSDelivery`, `_walk_to_trap`, `_apply_nbs_corrections!`, and
  `_propagate_correction` — the diff work now lives in the router's `:nbsout`
  signed-tracking step.
- `NBSElement`/`NBSPlan` hold the per-placement plan: `system`, `A`, `state_base`,
  `n_terrain`, precomputed **`I_bg`** = `Σ runoff[input_cell]`, and the output routes
  with precomputed `O_0` per endpoint and `ratio_e`.

### `DynFlowPath` + `build_network` (additive only — no termination)
- Add `DynFlowPath.nbs_inlets::Vector{Tuple{Int,Int}}` (nbs index, cell position),
  populated by `build_network` via `_intersecting_on_path` on `footprint_inflow_cells`
  — the same additive bookkeeping already used for `nbs_outlets`.
- `_path_event_templates` emits `:nbsin` (from `nbs_inlets`) and `:nbsout` (from
  `nbs_outlets`) alongside `:cvin`/`:cvout`.
- **No** `target_trap`/NBS-target change, **no** path termination, **no** `_network_order`
  node — the membership/trace machinery is untouched.

### Threading `runoff` into `_route_flow`
- `_route_flow` / `_path_delivered!` gain the per-cell `runoff` for each path (from
  `ctx.runoff`, read-only) and a signed diff accumulator handled by the rule above.
- `_routed_inflow` computes per-NBS `O_terrain(V)` / `E(V)` and the per-route diffs up
  front, passes them + `runoff` into `_route_flow`, and gets back `I_1` (the `:nbsin`
  draws + `I_bg`) for the layer-ODE block.

## Open implementation questions (small)
- Exact carrier for the signed diff through `_path_delivered!` — a second signed
  accumulator advanced per cell by `max(V+diff,0)−max(V,0)` alongside the positive
  `current`, delivered to `target_trap`. (`current` uses the existing infil-prefix
  rule; the diff uses the runoff-relative rule — they are separate quantities.)
- How `runoff` is passed per path: precompute per-path cell-runoff vectors from
  `ctx.runoff` (mirrors `path_infil_prefix`), read-only.

## Critical files
- `src/dynamics/networksolver.jl` — `_route_flow` / `_path_delivered!`
  (`:nbsin` draw → `I_1`; `:nbsout` emit the diff into the signed accumulator; per-cell
  `runoff` reference), `_path_event_templates` (`:nbsin`/`:nbsout`), `_routed_inflow`
  (compute `O_terrain`/`E`/diffs up front, thread `runoff`, return `I_1`), the NBS
  plan, and the rate-function NBS block (layer ODE driven by `I_1`); **delete**
  `NBSDelivery`, `_walk_to_trap`, `_apply_nbs_corrections!`, `_propagate_correction`.
- `src/dynamics/elements.jl` — add `DynFlowPath.nbs_inlets` (additive; **no**
  `target_trap` change). `src/dynamics/build_network.jl` — populate `nbs_inlets` via
  `_intersecting_on_path` on `footprint_inflow_cells` (like `nbs_outlets`); **no**
  path termination.
- `src/dynamics/nbs_elements.jl` — `NBSLayer` typed (**done**).
- Reuse: `DynNBSPlacement` cell lists (elements.jl); `DynFlowPath.nbs_outlets` +
  `_add_accumulation_traps!` region→trap map (network_utils.jl). `watercourses.jl`,
  `_network_order`, `network_updating.jl`, and the membership layer are **untouched**.

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
