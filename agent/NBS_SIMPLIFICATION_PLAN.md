# NBS code: assessment and simplification plan

Status: **steps 1–5 done, 2026-11-12** (`ae22903`, `354757b`, `0cc54e3`, `f922f95`, `030b8d2`). Step 6 (model change) was not pursued, per D-a. The current design is described in `agent/NBS_INVARIANTS.md`; this file records the reasoning.

**Decisions (user, 2026-11-12):**
- **D-a:** restructure first, keeping the signed-diff model. Change the model only if very strong arguments emerge, to keep the line between base SWIM and the bolted-on NBS module sharp.
- **D-c:** keep the evapotranspiration fields (no significant cost).
- **D-e:** archive the outdated NBS documents, clearly flagged as obsolete.
- **D-b:** a single rule for carrying an NBS's flow change downstream, the physically sensible one: the solver's per-cell `max(V+d,0) − max(V,0)`, the same rule base SWIM's `_track_flow!` uses. Added water is absorbed by spare infiltration capacity, and a removal cannot push a cell below zero. `network_watercourses` will share it.
- **D-d:** keep `n_terrain` general; make the outlet-to-layer pairing explicit (§3.2).

---

## 1. What the NBS code does today (one-page model)

An NBS *placement* is three things:
- a footprint (grid cells);
- a stack of storage layers (`NBSSystem`), whose ODE states are appended to the network's trap volumes;
- optional piped outlets.

**Design choice: the rest of SWIM stays NBS-oblivious.**
- `watercourses` and the `fill_sequence` rate machinery route water *through* footprints as if the NBS were not there. This is the *oblivious baseline*.
- The network solver then corrects that baseline. Per placement, at each RHS evaluation:

```
I_1   = O_0 + draw                 live input to the top layer
          O_0  = oblivious throughput: rain on footprint + Σ runoff at footprint-inflow cells
          draw = network flow intercepted at the inflow cells while routing (signed)
layers: dS_l/dt = in_l − overflow_l − percolation_l      (in_1 = I_1, in_l = percolation_{l−1})
out   = overflow of the top n_terrain layers     → back onto the terrain
        overflow of each lower layer             → its piped outlet
emit  = (out − O_0)·ratio_e  at each terrain exit e   (signed "diff")
        + piped overflow     at each piped outlet
```

**How the diff reaches the network.**
- Each exit's diff is injected at the head of the flow path leaving that exit. It is attenuated along that path against the baseline runoff (`max(V+d,0) − max(V,0)` per cell) and delivered to the path's target trap.
- An internal depression (a trap bottom inside the footprint) gets its share straight into the trap.
- Because the trap's `external_inflow` already contains the baseline `O_0`, the trap ends up with `O_0 + diff = out`.

**Network construction.**
- `setup_network` computes and stores each placement's inflow cells, exit cells and accumulation cells, and seeds the network at the exits and outlets.
- `_nbs_coupled_nodes` forces everything an NBS touches into one component.
- Layer storage persists across events and weather periods in `nbs_state[id]`.

**Model status.** The model is sound and, after `d0e296a`, mass-consistent on the tested cases. The difficulty is the **code**, not the model.

---

## 2. Why it is hard to read

Where the code lives: about 650 lines of NBS-specific code across 9 files (`nbs_elements`, `elements`, `build_network`, `network_utils`, `network_updating`, `networksolver`, `network_context`, `network_watercourses`, `utils`), plus `fill_sequence`.

| # | Problem | Where |
|---|---|---|
| P1 | **The same NBS terms are computed twice**, in two slightly different ways. `O_0`, the exit weights and ratios, the diff, the piped outlets and the layer fluxes are all rebuilt independently by the solver plan and by `network_watercourses`. The two can drift apart. They already differ: the solver attenuates the diff along paths, while `network_watercourses` adds it unattenuated. | `_build_nbs_plan` + `_nbs_routing` vs. `network_watercourses` |
| P2 | **Unit conversion `V·1000/A` … `·A·1e-3` is written out 7 times.** Meanwhile `NBSLayer.A` exists but is never used; the area comes from the footprint cell count. | solver ×5, `network_watercourses` ×2 |
| P3 | **`NBSRouting` is a two-way scratch struct over three index spaces.** It holds paths in, traps in, and placements out. The docstring has to explain why the two directions must not meet. | `networksolver.jl` |
| P4 | **Derived data is stored on the user's input object.** `setup_network` mutates `DynNBSPlacement`: it writes `id`, inflow, outflow and accumulation cells. So a "setup" call mutates its arguments, and the placement mixes input with cache. | `elements.jl`, `build_network.jl` |
| P5 | **Implicit positional conventions.** The k-th piped outlet belongs to the k-th outflowing layer below `n_terrain`. `state_base` offsets are computed in two places (`_build_nbs_plan`, `_nbs_layer_block`). `id = position in the nbs vector`. | several |
| P6 | **Coupling rules are spread out.** Which paths and traps must share a component with an NBS is decided in `_nbs_coupled_nodes` and `_add_accumulation_traps!`. Seed selection is in `setup_network`. Inflow-cell registration on paths is in `_intersecting_on_path`. Re-rooting is in `_promote_tributary!`. | 4 files |
| P7 | **Dead or half-used code.** `DynFlowPath.nbs_outlets` is never read by the router (coupling only). The `nbs_out` from `_intersecting_culverts_and_nbs_outlets` is discarded. `NBSLayer.A` is unused. Evapotranspiration fields are carried but ET is fixed at 0. `compute_outflow(…; mollifier)` is calibration-only. | various |
| P8 | **Mass leaks at the edges** (report C7). `max(S−Smax, eps())^n` makes an empty layer emit flow. The RHS floor zeroes `dS` at `S≤0` while the outflow still leaves. `max(O_0+draw, 0)` discards a deficit. | solver, `nbs_elements` |
| P9 | **"Zero infiltration on footprints" is a contract enforced by copy-paste** in 4 places, each writing `infil[footprints] .= 0`. | `fill_sequence`, `utils`, `network_watercourses`, `build_network` (validation) |
| P10 | **Docs describe four successive designs.** Only `NBS_ROUTING_NODE_PLAN.md` is current, and it predates the July fixes. The comments carry a lot of design history. | `agent/NBS_*.md`, docstrings |

---

## 3. Proposed structure (keeps today's model and today's results)

The goal is to give each concern **one owner and one representation**. Data flows one way: placement → link → plan → RHS.

### 3.1 `nbs_elements.jl`: the layer model only (pure, no terrain)

- `NBSLayer`, `NBSSystem` and the three constructors, unchanged in meaning.
- Keep `NBSLayer.A` and `compute_outflow`'s signature: SUrbArea's (deprecated) `NBS.jl` and calibration scripts use them. Document `A` as unused by the network solver, which takes the area from the footprint. Keep the ET fields, documented as "not yet used".
- **One flux function.** It works in volumes, so callers never see mm:
  ```julia
  layer_fluxes(L::NBSLayer, V, A) -> (overflow, percolation)   # m³/time; S_mm = V·1000/A inside
  ```
- **One RHS for a stack.** The top layer receives `inflow`. It returns the total overflow of the top `n_terrain` layers plus the per-layer piped overflow:
  ```julia
  nbs_stack_rates!(dV, system, V, A, inflow, n_terrain) -> (terrain_out, piped_out)
  ```
  This becomes the only place that knows the cascade. The solver's RHS and `network_watercourses` both call it.
- Fix P8 here: no `eps()` floor, and an empty layer emits nothing.

### 3.2 `DynNBSPlacement`: input only

- **Fields:** `system`, `footprint`, `n_terrain`, `outlets`, `id`.
- **`id`** is supplied by the caller or assigned once by `fill_sequence`, not by `setup_network`.
- **Outlet pairing is explicit.** The constructor validates the outlet count and stores `piped_layers::Vector{Int}`, the layer index for each outlet. This removes the positional convention (P5).

### 3.3 New `nbs_terrain.jl`: the coupling between a placement and the terrain (pure functions)

```julia
struct NBSLink                      # derived once per placement per tstruct; immutable
    inflow_cells ::Vector{Int}      # outside cells draining into the footprint
    exits        ::Vector{Int}      # outside cells the footprint drains to  (terrain re-emission)
    ponds        ::Vector{Int}      # footprint ∩ trap bottoms               (internal depressions)
end
nbs_link(tstruct, placement) -> NBSLink                      # replaces the mutation in setup_network (P4)
nbs_baseline(link, placement, runoff, precip) -> (O_0, shares::Vector{Float64})   # shares over exits ∪ ponds, Σ = 1
```

- `nbs_baseline` is the **single** home of `O_0` and the exit weights. Both the solver plan and `network_watercourses` call it (P1).
- The links are computed once per `fill_sequence` call. The footprint geometry depends only on `tstruct`, so the grid-wide inverse-flow map (report D4) goes away: inflow cells come from an 8-neighbour scan.

### 3.4 Solver: a plan built from links, and one-way routing

- **Plan.** `NBSElement` holds `A`, `state_base`, `O_0`, `(endpoint, share)` and `(carrier path, layer)`. It is built once per context (already true after `d0e296a`). `state_base` is computed only here; `_nbs_layer_block` reads it from the plan (P5).
- **Split `NBSRouting` (P3):**
  1. `nbs_injections(V, plan) -> (path_head::Vector, trap_add::Vector)`, computed **before** routing from the layer state.
  2. The router returns `draw::Vector` per placement as a plain output, next to the trap inflows.

  The no-algebraic-loop property then follows from the call order, with no explanatory docstring needed.
- **RHS.** Traps, then `nbs_stack_rates!` with `inflow = O_0 + draw`. One place, one conversion.

### 3.5 `network_watercourses`: reuse, don't recompute

- Use `nbs_baseline` and `nbs_stack_rates!` (or a flux-only variant) for the NBS terms.
- **Decision needed (D-b below):** attenuate the propagated diff as the solver does, or document why the snapshot field is additive.

### 3.6 Coupling in one place

- `nbs_coupled_nodes(net, links)` is the only function that knows what an NBS ties together: the paths carrying its exits and outlets, the paths crossing its inflow cells, and its pond traps.
- Seed selection in `setup_network` comes from the same `NBSLink` (exits, outlets, ponds), so seeds and coupling cannot disagree.
- Drop `DynFlowPath.nbs_outlets` and the `nbs_out` byproduct (P7). Coupling reads outlets from the placement.

### 3.7 Footprint infiltration contract (P9)

- **One helper**, `_nbs_zero_infiltration!(infil, nbs)`, called wherever an infiltration grid is built.
- The contract is listed in `NBS_INVARIANTS.md`.

### Expected effect

- Fewer NBS-specific lines (estimate 20–30%, mainly from P1/P2/P3/P7).
- One file per concern: layer model / terrain link / solver plan.
- No duplicated physics.
- Every positional convention made explicit.
- Comments limited to *what* and *why*, with design history moved to `NBS_INVARIANTS.md`.

---

## 4. Decisions needed before starting

**D-a: model, keep or change?**
- **Option 1: keep (recommended first step).** Keep today's signed-diff model and restructure as in §3. Results are bit-for-bit, or within ODE tolerance, of today's, so the existing tests act as the safety net.
- **Option 2: change.** Make the oblivious baseline treat footprints as **sinks**: flow entering a footprint, and rain on it, terminates there. Then NBS output is only ever *added* (positive injections). This removes `O_0`, the diff, signed attenuation and negative draws, which is the conceptually hardest part.
  - **Cost:** `watercourses` / `fill_sequence` flow tracking must accept extra sink cells, so the hand-written core loses its NBS-obliviousness. Exit weights are still needed to split terrain output, but could be purely geometric.
  - **Recommendation:** decide after Option 1, when the code is small enough to judge it.

**D-b: `network_watercourses` NBS diff.** Should the snapshot field attenuate the diff along paths the way the solver does, or stay purely additive?
- Today they differ, so the viewer's flow field near an NBS is not the solver's flow.
- Recommendation: attenuate, sharing the solver's rule.

**D-c: ET fields.** Keep the evapotranspiration fields (documented as unused) or remove them until ET is implemented?

**D-d: `n_terrain` generality.** All three models use the default `n_terrain = 1`. Keep the parameter, or fix "top layer re-emits on terrain, lower outflowing layers are piped"?

**D-e: stale design docs.** Delete or archive `NBS_INTEGRATION_PLAN.md`, `NBS_OPTION1_OVERLAY_PLAN.md`, `NBS_INTERIOR_REGION_PLAN.md`, `NBS_ROUTING_REDESIGN.md` and `NBS_REMAINING.md` once `NBS_INVARIANTS.md` exists?

---

## 5. Draft `NBS_INVARIANTS.md` content

This is the list the restructured code would assert at the point stated, and document in one place.

**Inputs**
- I1. The footprint is a non-empty set of in-domain cells, disjoint from every other placement's footprint.
- I2. There are no culvert endpoints, `dyn_coords` or sinks inside a footprint (a sink would end flow at a cell the plan never sees). *Asserted in `setup_network`.*
- I3. Infiltration is zero on every footprint cell. The layers own that water. *Applied by one helper; asserted in the plan builder.*
- I4. `length(outlets) == number of layers below the top n_terrain with Kout > 0`, and `outlets[k]` serves `piped_layers[k]`. *Constructor.*
- I5. `id` values are unique and stable for the whole `fill_sequence` call. They key `nbs_state`.

**Terrain link and baseline**
- L1. `O_0 = rain on footprint + Σ_{c ∈ inflow_cells} max(runoff_c, 0)`. This is exact because of I3: what enters equals what leaves.
- L2. The shares over `exits ∪ ponds` are non-negative and sum to 1. If `O_0 ≈ 0` they are split evenly. With non-zero throughput, a share is that endpoint's interior throughput divided by `O_0`, not the exit cell's own runoff. *Asserted.*
- L3. Every exit and every piped outlet is the departure point of exactly one network flow path (its carrier). Every pond lies in exactly one network trap. *Asserted in the plan builder.*

**Coupling**
- C1. Every path that crosses a placement's inflow cells, departs from its exits or outlets, or contains its pond traps is in the **same component** as the placement. This is what makes `draw` complete in a single routing pass.

**Solver**
- S1. The state layout is the trap volumes `1:nt`, then each placement's layers in `net.nbs` order, top to bottom, starting at `nt + state_base + 1`.
- S2. Injections depend only on the layer state, and draws feed only `dS`. There is no algebraic loop and no routing order among NBS.
- S3. Mass: `Σ dS_layers = I_1 − terrain_out − piped_out − ground loss`, where ground loss is the bottom layer's percolation. *Unit test.*
- S4. An empty layer emits nothing (`S ≤ 0 ⇒ overflow = percolation = 0`). No flux floors that break S3.
- S5. A context's rates (inflow, baseline, shares) are frozen at its last touch (`d0e296a`).

**Persistence**
- P1. Layer storage is carried across events and weather periods only through `nbs_state[id]`, written at every harvest and finalize and read at every context build.

---

## 6. Suggested sequence (each step behaviour-preserving, tests green after each)

1. Layer model: `layer_fluxes` / `nbs_stack_rates!`, used by the solver RHS and `_nbs_routing` (P2, P8). C7 is fixed here, and the expected values move with the fix.
2. `NBSLink` + `nbs_baseline`; remove the mutation from `setup_network`; switch the solver plan and `network_watercourses` to them (P1, P4, D4).
3. Split `NBSRouting` into injections and draws (P3); make `state_base` single-sourced (P5).
4. Single coupling function; drop `nbs_outlets` and the dead outputs (P6, P7); add the infiltration helper (P9).
5. Write `NBS_INVARIANTS.md`; trim the docstrings to match it; retire the stale design docs (P10).
6. Optionally revisit D-a, Option 2.
