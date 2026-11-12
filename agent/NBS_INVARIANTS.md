# NBS module: terms, model, code map and invariants

This is the authoritative description of how Nature-Based Solutions (NBS) are coupled to SWIM. Code comments refer here rather than repeating it. Older design documents are in `agent/archive/` and are obsolete.

New to the NBS code? Start with `NBS_INTRO_OVERVIEW.md` (main logic and workflow), then `NBS_INTRO_DETAILED.md` (function-by-function walkthrough).

## Terms

Each term below is defined using only base SWIM vocabulary and terms defined above it.

**The placement**

- **NBS placement** (`DynNBSPlacement`): an NBS installed on the terrain. It consists of a footprint, a layer stack, `n_terrain` and outlets, each defined just below.
- **footprint**: the grid cells the placement covers. Always "NBS footprint" here; a trap's footprint is called *trap footprint*.
- **layer stack** (`NBSSystem`): storage layers, top first. Water enters the top layer.
  - **storage** `S`: the water a layer holds, in mm (`V` when expressed as a volume).
  - Each layer has two thresholds: `Smax` controls its overflow and `Smin` its percolation.
- **overflow**: what a layer sheds above its `Smax`: `Kout·(S − Smax)^nout` (coefficient `Kout`, exponent `nout`), zero at or below `Smax`.
- **percolation**: what a layer passes to the layer below above its `Smin`: `Kinf·(S − Smin)^ninf` (coefficient `Kinf`, exponent `ninf`), zero at or below `Smin`. The bottom layer's percolation leaves the model ("to ground").
- **`n_terrain`**: the number of top layers whose overflow goes back onto the terrain.
- **piped layer**: a layer below the top `n_terrain` that has overflow (`Kout > 0`). The list is `_piped_layers(system, n_terrain)`, top first.
- **outlet**: the cell where a piped layer's overflow (its pipe flow) re-enters the terrain. `outlets[k]` belongs to the k-th piped layer.
- **`id`**: the placement's position in the `nbs` vector passed to `fill_sequence`. It labels the placement's stored layer volumes.

**The terrain and the baseline**

- **baseline** (*oblivious* in code): the flow field and trap inflows base SWIM computes. It treats footprints as ordinary terrain.
- **runoff value** `R`: the baseline runoff grid holds, per cell, either positive surface flow or, if negative, unused infiltration capacity.
- **inflow cells**: cells outside a footprint whose downstream cell (in the flow graph) is inside it.
- **exits**: cells outside a footprint that a footprint cell drains to.
- **ponds**: trap bottoms inside a footprint.
- **leaf regions** of a trap: the lowest-level spill regions it contains (the trap itself, if it has no sub-traps).
- **`O_0`**: the footprint's baseline throughput: rain on the footprint plus the baseline flow entering it through its inflow cells.
- **share** (`share_e`): the fraction of `O_0` that leaves the footprint through endpoint `e`, where an endpoint is an exit or a pond.
  - **For an exit:** the baseline flow of the footprint cells that drain straight out to that exit, divided by `O_0`. The exit cell's own runoff is not used, because it also includes water that never crossed the footprint.
  - **For a pond:** the baseline flow arriving at the pond, divided by `O_0`.
  - The shares sum to 1, and are set equal if `O_0 ≈ 0`.

**The dynamic network**

- **culvert** (`DynCulvert`): a pipe moving water from an inlet cell to an outlet cell, at a rate set by the water levels at both ends. Like an NBS, it needs the dynamic network because its rate changes continuously with the water levels.
- **seed**: a cell from which flow is traced. Seeds are each NBS's exits, outlets and ponds, each culvert's two ends, and one cell of each trap named in `fill_sequence`'s `dyn_traps` argument (passed internally as `dyn_coords`).
- **network path** (`DynFlowPath`): a run of cells traced downhill from a seed, or from a full trap's spillpoint, to the next trap (or out of the domain).
  - **departure point:** where it starts. **head:** its first cell.
  - A path can end where it joins another path; it is then a **tributary** of that path.
- **network trap**: a trap reached by such tracing, where tracing continues through full traps along their spill paths. Its volume is an ODE state. The set changes during a run:
  - it grows when a network trap fills and its spill path is traced further;
  - traps drop out when they are no longer reached.
- **covered trap**: a network trap, or a full sub-trap inside one. Base SWIM leaves covered traps alone: no changetime estimates, and they are not in its spill graph.
- **component** (a `DynNetwork`): a set of network paths, network traps, culverts and NBS placements that are linked (by flow, culverts or NBS) and so are solved as one ODE. One run has several independent components.
- **carrier path**: the network path departing from an exit or an outlet. NBS output enters at its head.
- **routing pass**: one evaluation of "what arrives where" in a component (`_route_flow`), done at every evaluation of the ODE right-hand side. It visits paths and traps from upstream to downstream.
- **network flow**: what the routing pass carries along network paths. It is always an *increment on top of the baseline*. Examples are the spill of a full network trap (base SWIM's spill graph excludes covered traps, so this spill exists only here), a culvert's flow, and NBS corrections. It is signed.
- **per-cell rule**: how network flow `d` changes along a path. At each cell with baseline runoff value `R`, `d ← max(R+d, 0) − max(R, 0)`. Unused infiltration capacity absorbs added water, and a negative `d` can remove at most the flow that is there. This is the same rule base SWIM's `_track_flow!` uses when it routes a trap's spill.
- **draw**: the network flow reaching a placement's inflow cells during one routing pass. It is signed: negative when an upstream NBS is holding water back.
- **injection**: what a placement adds to a routing pass: `(O_t − O_0)·share_e` at each exit's carrier path or pond's trap, and each piped layer's pipe flow at its outlet's carrier path. `O_t` is defined in the Model section below.

**Running the network**

- **context** (`DynNetworkContext`): one component together with its *committed* state (trap volumes and layer volumes at a known time) and its rates (`params`).
- **touch**: what happens when an event in `fill_sequence` affects a component (one of its traps changed status, or the baseline inflow to one of its network traps' leaf regions changed):
  - **commit:** integrate the context's ODE up to the event time;
  - **harvest:** read the committed volumes out;
  - apply the change (the component may grow, lose traps, or merge with another);
  - **rebuild:** build a new context with rates taken from the current baseline;
  - **predict:** integrate a copy ahead to find its next event.
- **frozen rates**: a context's `params` are built once at its touch and never re-read from the live baseline grids until the next touch.
- **finalize**: commit every context to the end of a weather period.
- **`nbs_state`**: the stored layer volumes of every placement, by `id`. It is carried across touches and weather periods.

## Model

**Base SWIM stays NBS-oblivious.** `watercourses`, `spillanalysis` and the constant-rate `fill_sequence` machinery never see an NBS.

**`external_inflow`**: each network trap's baseline inflow. It is base SWIM's `rateinfo.trap_inflow`, summed over the trap's leaf regions (`_external_inflow`). It includes, after infiltration on the way, the throughput `O_0` of any footprint upstream, because the baseline routes water straight through footprints. The routing pass starts every network trap from `external_inflow` and adds network flow on top.

**Per placement, at every evaluation of the ODE right-hand side:**

```
Q_in  = O_0 + draw                                   input to the top layer
dV_l  = in_l − overflow_l − percolation_l            in_1 = Q_in,  in_l = percolation_{l−1}
O_t   = overflow of the top n_terrain layers         (back onto the terrain)
E_k   = overflow of the k-th piped layer             (pipe flow, to outlets[k])
```

**Injections into the routing pass:**
- **Exits:** `(O_t − O_0)·share_e` at the head of each exit's carrier path.
- **Ponds:** `(O_t − O_0)·share_e` into each pond's trap.
- **Outlets:** `E_k` at the head of each outlet's carrier path.

Both the injection and the baseline `O_0·share_e` travel the same path under the per-cell rule. The trap downstream of exit `e` therefore receives `O_t·share_e` instead of `O_0·share_e`: what the NBS actually releases.

## Code map

| Concern | File | Key functions |
|---|---|---|
| Layer model (mm, no terrain) | `dynamics/nbs_elements.jl` | `NBSLayer`, `NBSSystem`, model constructors, `compute_outflow`, `_layer_fluxes`, `_nbs_layer_rates!`, `_nbs_terrain_outflow` |
| Placement (input) | `dynamics/elements.jl` | `DynNBSPlacement`, `_piped_layers`, `NBSLink` |
| Placement ↔ terrain | `dynamics/nbs_terrain.jl` | `_prepare_nbs!`, `_nbs_link`, `_nbs_baseline`, `_pond_trap`, `_nbs_zero_infiltration!` |
| Network membership | `dynamics/build_network.jl`, `dynamics/network_utils.jl` | seeds in `setup_network`; `_intersecting_on_path`; `_nbs_coupled_nodes` |
| Solver | `dynamics/networksolver.jl` | `_build_nbs_plan` → `NBSPlan`/`NBSElement`; `_nbs_injections`; `_route_flow`; `dynNetworkRateFunction!` |
| Persistence | `dynamics/network_context.jl` | `_nbs_layer_block` (restore), `_store_nbs_state!` / `_nbs_layers` (save) |
| Flow field at a timepoint | `dynamics/network_watercourses.jl` | `network_watercourses` |

## Invariants

Each invariant names where it is enforced.

**Inputs**
- **I1.** Footprints are non-empty and do not overlap. *`_validate_network_inputs` (error).*
- **I2.** No culvert end, `dyn_traps` cell or sink (a cell where SWIM removes water) lies inside a footprint. A sink would end flow at a cell the placement never sees. *`_validate_network_inputs` (error).*
- **I3.** Infiltration is zero on every footprint cell, because the layers own that water. *`_nbs_zero_infiltration!`, applied wherever an infiltration grid is built.*
- **I4.** `length(outlets)` equals the number of piped layers. *`DynNBSPlacement` constructor (error).*
- **I5.** `id` is stable for a whole `fill_sequence` call. *Set by `_prepare_nbs!`.*

**Terrain link and baseline**
- **L1.** `O_0 = rain on footprint + Σ_{c ∈ inflow cells} max(R_c, 0)`. This is exact because of I3: everything entering also leaves.
- **L2.** The shares are ≥ 0 and sum to 1. *Asserted in `_nbs_baseline`: the unnormalised weights must sum to `O_0`.*
- **L3.** Every exit and outlet is the departure point of a network path (so it has a carrier path), and every pond lies in exactly one network trap. *Asserted in `_build_nbs_plan` and `_pond_trap`.*

**Network membership**
- **C1.** A placement is in the same component as the carrier paths of its exits and outlets, every path crossing its inflow cells, and its ponds' traps. This is what makes its draw complete within one routing pass. *`_nbs_coupled_nodes`.*

**Solver**
- **S1.** The ODE state is laid out as the component's `nt` trap volumes, then each placement's layer volumes in the component's placement order (`net.nbs`), top layer first. *`_build_nbs_plan` records each placement's offset (`state_base`); `_layer_range` reads it; `_nbs_layer_block` builds the state in the same order.*
- **S2.** Injections depend only on the layer volumes, and draws feed only the layers' rates. So no placement needs to be routed before another, and there is no circular dependency within a routing pass. *Call order in `_routed_inflow`.*
- **S3.** Per placement, `Σ dV = Q_in − O_t − Σ E_k − bottom percolation`. *Unit test "NBS: layer cascade conserves mass".*
- **S4.** A layer that is empty, or at or below a threshold, emits nothing through that flux. No rate is clamped ("floored") after the fact. *`compute_outflow`, `_layer_fluxes`.*
- **S5.** A commit uses the context's frozen rates. *`DynNetworkContext.params`; fixed in `d0e296a` (review item C14); test "a network commit is not altered by the event that triggers it".*

**Persistence**
- **P1.** Layer volumes cross touches and weather periods only through `nbs_state[id]`. They are written at every harvest and finalize, and read at every rebuild.

## Known simplifications (marked `# @@@` in code)

- **Cell area.** The footprint area is the cell count, i.e. 1 m² per cell.
- **Evapotranspiration (ET).** Not modelled. The layers' `EVCoeff`/`EVS11` are carried but unused.
- **Ponds and the field.** In `network_watercourses`, a pond's share stays in its trap, even when that trap is full and spilling.

## Open question (to verify)

- **Possible double count at a weather-period start.** Base SWIM's spill graph should exclude covered traps (`_reconcile_spillgraph!`). That reconciliation runs when coverage changes at an event, but not when a weather period's network is first built. If a trap is already full at that moment and covered by a network, its spill may be routed both by the baseline and by the network. An example would be a zero-volume trap downstream of an NBS, which counts as full from the start.
