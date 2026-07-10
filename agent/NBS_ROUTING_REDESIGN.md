# NBS routing redesign — design notes

Status: design settled, implementation in progress.
Scope: how Nature-Based Solution (NBS) installations couple to the dynamic
surface-flow network. Supersedes the overlay-element approach in
`NBS_OPTION1_OVERLAY_PLAN.md` for the routing/coupling layer (the NBS *layer
model* primitives — `NBSSystem` / `NBSLayer` / `compute_outflow` — are retained).

The new network-build entry point is `src/dynamics/build_network.jl`
(`setup_network`); the flux/deduction/layer-routing consumer is being
re-implemented separately in the solver / rateinfo layer.

---

## 1. Core idea

An NBS is **invisible to the surface-flow network**. The fixed overland flow
(the baseline that sets regular trap fill times) is routed across the NBS
footprint exactly as if no NBS were present. The NBS then acts *only* as a
correction:

- it captures ("infiltrates") some water on its footprint,
- stores it in its layered model (the `NBSSystem` ODE),
- re-emits it at explicit **outlets**, from where ordinary flow paths carry it
  downstream.

The practical effect of an NBS is therefore entirely a **deduction applied as a
correction** to the fixed flow — negative flow where water is captured, positive
flow (re-emission) at the outlets — never a rewrite of the surface network
topology.

This mirrors how dynamic flow paths already correct the fixed flow elsewhere: a
static baseline, plus a dynamic correction where the network is active.

## 2. Water balance

For one NBS over a solve window:

```
capture  C  =  precip on footprint  +  inflow from inflow-cells
storage  ΔS =  change in layer storage
emission E  =  Σ outlet discharge (from the layer model)

C = ΔS + E                         (NBS internal balance)
```

The correction applied to the surface:

- **deduct** the capture `C` from the fixed flow, spread over the NBS's
  **outflow-boundary cells** and **internal-accumulation cells** (see §4);
- **add** the emission `E` at the **outlet cells**, injected into the flow paths
  that are seeded there.

Net water removed from the surface at any instant is `C − E = ΔS`, i.e. exactly
what the NBS is currently holding. At steady state `ΔS = 0`, so deduction equals
emission and mass closes. During filling, more is deducted than emitted (water
held back); during draining, the reverse. Mass is conserved at all times.

## 3. The four cell-sets (computed in `build_network.jl`)

`_compute_nbs_inflow_outflow_cells!` classifies each footprint's coupling cells:

| Set | Definition | Role |
|-----|-----------|------|
| `footprint_inflow_cells`  | cells feeding *into* the footprint from outside (inverse-flow neighbours of footprint cells, minus the footprint itself) | where capture `C` is sourced |
| `footprint_outflow_cells` | footprint cells whose downstream neighbour lies *outside* the footprint | a deduction location; also seeded (at the external `ds` cell) so re-routing continues downstream |
| `internal_accumulation_cells` | footprint ∩ trap-bottoms | a deduction location; feeds any trap whose footprint intersects here |
| `outlets` (user-supplied) | piped-discharge cells for the outflowing layers below the top `n_terrain` | where emission `E` re-enters the terrain |

Contract: **NBS footprints have zero infiltration** (see NOTES-TO-SELF in
`build_network.jl`) so the footprint contributes no losses of its own and the
only abstraction on it is the NBS capture.

## 4. The deduction (correction) — heuristic

The NBS internals are not spatially modelled, so the capture `C` is distributed
back onto the surface by a simple heuristic:

> Split the negative correction **equally** across the union of
> (`footprint_outflow_cells` ∪ `internal_accumulation_cells`).

**Per-cell clamp is mandatory.** An equal share can exceed the baseline flow
actually present at a thin-flow cell. Deducting more than is available there
would drive negative surface flow and silently leak mass. So the distributor
must:

1. first pass: assign each cell its equal share;
2. clamp each deduction to the flow available at that cell;
3. redistribute the un-deducted remainder across cells that still carry flow
   (or carry it forward as unmet NBS storage).

This is the same rule the culvert/tributary router already obeys
(`AGENTS.md` → Mass conservation): draw the *available* amount, never the
*requested* amount.

- Deduction at an **outflow-boundary** cell reduces the flow leaving the
  footprint there.
- Deduction at an **internal-accumulation** cell reduces the inflow to the
  **trap** that contains that cell (which may not be downstream of the NBS at
  all — it simply shares footprint with it).

## 5. Emission (re-emission at outlets)

Outlets are added to the network as **seeds**. A path grown from an outlet cell
is an ordinary flow path; the NBS discharge is injected at the outlet position
and routed downstream to whatever trap the path leads into — no NBS-specific
routing machinery beyond the injection.

**Outlet record shape.** A path records its NBS outlet as
`DynFlowPath.nbs_outlets`. Default resolution for now: **all outflowing
(non-runoff) layers discharge at the same single outlet point**, so the current
2-tuple `(placement_ix, position)` suffices — the consumer sums all piped-layer
discharge and injects it at that one position. When per-layer distinct outlets
are needed later, extend to `(placement_ix, layer_ix, position)`. (Marked as a
`# @@@` extension point at the field.)

## 6. NBS inflow — static baseline, dynamic correction

Capture `C` sources from precipitation on the footprint plus inflow-cell
contributions. An inflow cell is handled like any other terrain cell:

- computed **statically** as a baseline (precip + static inflow-cell flow);
- **corrected dynamically** when the inflow cell lies on an active flow path
  (then its contribution is whatever the network delivers there).

Whether a given inflow cell is "static" or "dynamic" is only knowable once the
paths are built, so it is treated as static and corrected where a path passes
through it — exactly the existing flow-path convention.

## 7. Division of labour

**`build_network.jl` (topology):**
- identify the four cell-sets per NBS;
- add seeds for outlets, outflow-boundary (external `ds` cell), and
  internal-accumulation cells;
- grow the monolithic network, then split into connected components.

**Distributor / rate layer (flux — re-implementation in progress):**
- compute capture `C`, storage, and emission `E` from the layer model;
- apportion the negative correction over (outflow ∪ accumulation) cells with the
  clamp+redistribute of §4;
- land each share on its network element: an outflow-boundary share on the
  seeded path, an accumulation share on the **containing trap**.

### Resolving an accumulation cell to its trap (solver-side)

No persisted cell→trap map is needed. The trap is already a `DynTrap` in the
`DynNetwork` (the accumulation seed — a trap bottom, hence a flow-graph sink —
traces to and registers it during build). At solve time the distributor:

- deducts from the `DynTrap` whose `footprint` **contains** the accumulation cell
  (robust against sub/supertrap nesting — always the node actually in the net;
  do **not** re-derive `regions[cell] → trap_ix` independently, it may pick a
  different trap than the one build registered);
- builds the cell→local-trap lookup **once per solve**, not `findfirst` per cell
  per rate evaluation (this runs inside the ODE hot path — per-event overhead
  matters).

## 8. Cycle safety

NBS piped outlets, like culverts, are not constrained to follow terrain flow
direction, so an outlet placed upstream of the NBS's own inflow could form a
loop and break the "water always flows downstream" invariant. Handled by
**user contract**: outlets must not be placed at a higher elevation than the
footprint. (No automatic detection for now — matches the deferred reverse-culvert
handling.)

## 9. Splitting into connected components

`split_network_into_connected_components(net, tstruct)` (public, in
`src/dynamics/network_utils.jl`) builds an undirected graph over path nodes
(`1:np`) and trap nodes (`np+1:np+nt`), takes `Graphs.connected_components`, and
rebuilds each as a standalone `DynNetwork` with all cross-references remapped to
local 1-based indices.

Connectivity must follow more than terrain flow — a culvert or an NBS bridges
otherwise-disjoint regions and must land them in the *same* component (shared
mass balance). Edges come from:

- `path → target_trap`, `trap → spill_path`, and `path` merges (terrain + tribs);
- **culverts**: unite the path/trap owning the inlet with the one owning the
  outlet;
- **NBS**: unite each NBS's coupled elements —
  - *emission*: paths carrying one of its outlets (`nbs_outlets`);
  - *outflow deduction*: paths whose **first cell** is one of the NBS's
    `footprint_outflow_cells` (the seed stores the external `ds` cell, so this is
    a `net`-only test);
  - *accumulation deduction*: the network trap holding its
    `internal_accumulation_cells`. Accumulation seeds register a trap (a sink)
    with no flow path, so they can't be matched cell-on-path; instead map each
    cell to its lowest-level region (`tstruct.regions[c]`, deduped — absorbs
    flat many-celled bottoms) and then to the one network trap in that region's
    `supertraps_of` hierarchy (asserted unique). This is the only NBS link that
    needs `tstruct`.

### Deferred

- Per-layer distinct outlets (§5) — single shared outlet for now.
- The distributor / rate-layer re-implementation itself.

---

## Relationship to the earlier implementation

The earlier approach modelled the NBS as a `DynNBS` **overlay element** carried
inside the `DynNetwork`, with its own layered state appended to the solver's
state vector and re-emission derived at terrain-exit boundaries plus piped
outlets. The redesign:

- **drops `DynNBS`**; the `DynNetwork` carries `NBSPlacement` directly (NOTES:
  possibly rename to `DynNBSPlacement` for consistency with `DynCulvert`), and
  `NBSPlacement` is now mutable so its build-time cell-sets can be filled in;
- reframes the NBS as a **pure correction to the fixed flow** rather than an
  overlay node — capture as negative flow at boundary/accumulation cells,
  emission as ordinary seeded flow paths at outlets;
- builds the network **monolithically** (grow from all seeds into one
  `DynNetwork`, dedup via a shared `pathmap`), then splits into connected
  components, instead of tracing per-seed subnetworks and merging with offset
  remapping.
