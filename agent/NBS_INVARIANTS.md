# NBS module: model, code map and invariants

This is the authoritative description of how Nature-Based Solutions are coupled to SWIM. Code comments refer here rather than repeating it. Older design documents are in `agent/archive/` and are obsolete.

New to the NBS code? Start with `NBS_INTRO_OVERVIEW.md` (main logic and workflow), then `NBS_INTRO_DETAILED.md` (function-by-function walkthrough).

## Model

**The base stays NBS-oblivious.** `watercourses`, `spillanalysis` and the constant-rate `fill_sequence` machinery never see an NBS. They route water *through* a footprint as if it were bare terrain. This field is the **baseline** (called *oblivious* in the code), and the traps' `external_inflow` comes from it.

**The NBS module corrects the baseline**, inside the dynamic network solver only. For each placement, at every rate evaluation:

```
O_0     baseline throughput: rain on the footprint + baseline flow entering it      (constant per touch)
draw    network flow reaching its inflow cells in this routing pass               (signed)
I_1   = O_0 + draw                      input to the top layer
dS_l  = in_l − overflow_l − percolation_l     in_1 = I_1, in_l = percolation_{l−1};
                                              the bottom layer's percolation leaves the model
O_t   = overflow of the top n_terrain layers                 → back onto the terrain
E_k   = overflow of the k-th piped layer                     → its outlet cell
```

**What the network receives.**
- Each exit `e` gets `(O_t − O_0)·share_e` at the head of the path leaving it.
- Each pond gets the same kind of term directly into its trap.
- Each piped outlet gets `E_k` at the head of its path.

Because `external_inflow` already contains `O_0`, a downstream trap ends up receiving `O_t` instead of `O_0`.

**Carrying a flow change downstream.** Any change of flow `d` travels down a path cell by cell, with `d ← max(V+d,0) − max(V,0)` at each cell (`V` = the cell's baseline value). Spare infiltration absorbs added water, and a removal cannot take a cell below zero flow. This is the same rule as base SWIM's `_track_flow!`, and it is used by both the solver and `network_watercourses`.

## Code map

| Concern | File | Key functions |
|---|---|---|
| Layer model (mm parameters, no terrain) | `dynamics/nbs_elements.jl` | `NBSLayer`, `NBSSystem`, the model constructors, `compute_outflow`, `_layer_fluxes`, `_nbs_layer_rates!`, `_nbs_terrain_outflow` |
| Placement (input) | `dynamics/elements.jl` | `DynNBSPlacement`, `_piped_layers`, `NBSLink` |
| Placement ↔ terrain | `dynamics/nbs_terrain.jl` | `_prepare_nbs!`, `_nbs_link`, `_nbs_baseline`, `_pond_trap`, `_nbs_zero_infiltration!` |
| Network membership | `dynamics/build_network.jl`, `dynamics/network_utils.jl` | seeds in `setup_network`; `_intersecting_on_path` (inflow cells on paths); `_nbs_coupled_nodes` |
| Solver | `dynamics/networksolver.jl` | `_build_nbs_plan` → `NBSPlan`/`NBSElement`; `_nbs_injections`; `_route_flow` (`:nbsin` draws); `dynNetworkRateFunction!` |
| Persistence | `dynamics/network_context.jl` | `_nbs_layer_block` (restore), `_store_nbs_state!` / `_nbs_layers` (save) |
| Flow field at a timepoint | `dynamics/network_watercourses.jl` | `network_watercourses` |

## Invariants

Each invariant names where it is enforced.

**Inputs**
- **I1.** Footprints are non-empty and disjoint from each other. *Enforced by `_validate_network_inputs` (error).*
- **I2.** No culvert endpoint, `dyn_coord` or sink lies inside a footprint. A sink would end flow at a cell the placement never sees. *Enforced by `_validate_network_inputs` (error).*
- **I3.** Infiltration is zero on every footprint cell, because the layers own that water. *Applied by `_nbs_zero_infiltration!` wherever an infiltration grid is built.*
- **I4.** `outlets[k]` is the outlet of layer `_piped_layers(system, n_terrain)[k]`, and the counts match. *Checked by the `DynNBSPlacement` constructor (error).*
- **I5.** `id` is the placement's position in the `nbs` vector, and is stable for a whole `fill_sequence` call. It keys the persistent layer storage. *Set by `_prepare_nbs!`.*

**Terrain link and baseline**
- **L1.** `O_0 = rain on footprint + Σ_{c ∈ inflow_cells} max(runoff_c, 0)`. This is exact because of I3: everything entering also leaves.
- **L2.** The shares over `exits ∪ ponds` are ≥ 0 and sum to 1. For an exit, the share is based on the baseline flow the footprint hands to it, not the exit cell's own runoff. When `O_0 ≈ 0` they are split evenly. *Asserted in `_nbs_baseline`: the unnormalised weights must sum to `O_0`.*
- **L3.** Every exit and every piped outlet is the departure point of a network path (its carrier). Every pond lies in exactly one network trap. *Asserted in `_build_nbs_plan` and `_pond_trap`.*

**Network membership**
- **C1.** A placement shares its component with every path that departs from its exits or outlets, every path that crosses its inflow cells, and the trap under each pond. This is what makes `draw` complete in one routing pass. *`_nbs_coupled_nodes`.*

**Solver**
- **S1.** State layout: trap volumes `1:nt`, then each placement's layers in `net.nbs` order, top layer first. *`_build_nbs_plan` records the offsets (`state_base`); `_layer_range` reads them; `_nbs_layer_block` builds them in the same order.*
- **S2.** NBS injections depend only on the layer state, and draws feed only the layers' rates. So there is no algebraic loop and no routing order among NBS. *Call order in `_routed_inflow`.*
- **S3.** Per placement, `Σ dS = I_1 − O_t − Σ E_k − bottom percolation`. *Unit test "NBS: layer cascade conserves mass".*
- **S4.** A layer at or below its threshold, or empty, emits nothing. There are no flux floors. *`compute_outflow`, `_layer_fluxes`.*
- **S5.** A network context's inflow, baseline and shares are frozen at its last touch. Commits never read live grids. *`DynNetworkContext.params` (fix `d0e296a`; test "a network commit is not altered by the event that triggers it").*

**Persistence**
- **P1.** Layer storage crosses events and weather periods only through `nbs_state[id]`. It is written at every harvest and finalize, and read at every context build.

## Known simplifications (`@@@` in code)

- **Cell area.** Footprint area is the cell count, i.e. 1 m² per cell.
- **Evapotranspiration.** Not modelled. `EVCoeff`/`EVS11` are carried but unused.
- **Ponds and the field.** In `network_watercourses`, a pond's share stays in its trap, even when that trap is full and spilling.
