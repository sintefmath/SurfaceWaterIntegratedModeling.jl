# NBS as an interior region with a natural spillpoint (Option 2) — REJECTED ALTERNATIVE

**Status: not chosen.** The selected approach is the overlay model in
`NBS_OPTION1_OVERLAY_PLAN.md`. This file is kept for reference — it lost mainly on
static risk, the submerge/emerge state handoff, and forced ponding on slopes.
Motivation and option analysis: `agent/prompts/nbs_integration_reconsideration.org`.

## Why the change

The current implementation digs each NBS footprint into a flat artificial trap
and bolts outlets onto a computed `Spillpoint`. Three defects:

1. **Green roofs on buildings** — carving the footprint destroys the building
   obstacle it should sit on.
2. **Non-level NBS** (sloped permeable pavement) — flattening the footprint
   discards the slope.
3. **Outlet leak-points in the static trap structure** — a dynamic-only concept
   smuggled into static analysis (`_resolve_nbs!` admits it doesn't even rewrite
   the spillpoint it validates against).

## Chosen model: interior region with a natural spillpoint

An NBS footprint becomes **one ordinary positive spill region on real terrain**
(no carving). It seeds a region via the `-4` flag (§3), but is *not* an exit
node, so `spillpoints` computes its spillpoint = **lowest boundary saddle** and
it participates in `sshierarchy` like any trap. The dynamic solver tags the trap
as NBS and overrides simple volume-filling with the layered model.

### Why interior, not exterior (design history)

An earlier draft made the footprint an *exterior* (negative, sink-like) region
to avoid giving it a spillpoint — motivated by "a spillpoint must be the lowest
boundary point, else water flows in and out the same saddle → cycle." That was
an over-correction: `spillpoints` **always** picks the lowest saddle by
construction, so a natural spillpoint is cycle-safe automatically — the same
guarantee every ordinary trap relies on. There was no problem to dodge.

Going exterior also *threw away submergence*: an NBS genuinely can be submerged
(a neighbor trap fills to the shared saddle and the two merge; e.g. a piped
outlet discharges into a trap that spills back over terrain into the NBS). That
merge is a sub/supertrap relationship, which requires the NBS to *have* a
spillpoint. Interior gives it one, so submergence + merging come from the
existing hierarchy + dynamic-network merge for free (§7) — no bespoke geometric
check.

Inflow stays trivial either way: the footprint's whole upstream catchment is one
region, so NBS inflow = rain-over-region + upstream.

## Core unifying idea

**Every layer's overflow is the same operation: deliver flow to a cell, from
which it flows downstream as normal.** There is no "runoff vs piped" distinction
in the model — a layer's outlet is just a cell (like a culvert outlet), routed
through the position-exact delivery built in PR-3.

The only difference is a *resolution-time default*: an **unspecified** outlet
backfills to the trap's **natural spillpoint downstream cell** (the surface
re-emergence point — the lowest saddle, already computed by `spillpoints`). An
**explicit** outlet is used as given. The natural spillpoint is cycle-safe by
construction; an explicit (piped) outlet is validated to be cycle-safe (§4).

---

## 1. Data model (`NBSPlacement`)

```julia
struct NBSPlacement <: DynObject
    system::NBSSystem
    footprint::Vector{Int}                # linear indices, unchanged
    outlets::Vector{CartesianIndex{2}}    # one per outflowing layer (Kout>0)
end
```

Same shape as culvert outlets. An unspecified entry (`CartesianIndex(0,0)`)
backfills to the footprint's lowest-saddle cell in `_resolve_nbs!`; an explicit
entry is validated and used as given. The solver only ever sees "layer L
delivers to cell C."

## 2. Static pipeline — no carve, no dig, no spillpoint bolt-on

New order in `spillanalysis`:

```julia
field, slope = spillfield(gridcpy, ...)          # unchanged — real slopes kept
_flag_nbs_cells!(field, nbs)                      # NEW: footprint cells → -4
regions, flowgraph, bottoms = spillregions(field, ...)
spoints, regbnd = spillpoints(gridcpy, regions, ...)   # skips negatives already
...
_resolve_nbs!(nbs, regions, spoints, gridcpy)     # runoff cell + piped outlets
```

Delete: `_dig_nbs_traps!`, the `Spillpoint` bolt-on, the supertrap cycle guard
(the natural spillpoint is cycle-safe by construction; the hierarchy handles
merging).

## 3. New spillfield flag `-4` (NBS capture cell) — two touch points

Behaves like `-1` (trap bottom), so the footprint forms an ordinary *positive*
region on real terrain:

- **`_spillfield_flow_edges!`** — like `-1`: no outward edge, push a self-edge.
  Streamlines terminate on the footprint (all surface water hitting it is
  captured — this is what forces a region even on a sloped footprint that isn't
  a natural depression).
- **`_flat_zone_connecting_edges!`** — connect neighboring `-4 ↔ -4` (like `-1`)
  so the footprint merges into **one** region. Assumes footprints are
  disjoint/non-adjacent (already assumed today).

Do **not** add `-4` to `_find_exit_nodes`. That keeps the region *positive*, so
`spillpoints` computes its spillpoint (lowest saddle) and `sshierarchy` includes
it. Also: **exempt `-4` regions from the zero-volume flat-trap elimination** in
`spillregions._process_domain!`, so a flat footprint isn't merged away.

## 4. Outlets: natural spillpoint (default) + validated piped outlets

The runoff outlet is just the trap's **natural spillpoint downstream cell** —
`spoints[nbs_region].downstream_region_cell`, already computed, already the
lowest saddle. `_resolve_nbs!` backfills any unspecified layer outlet to it. No
bespoke saddle helper is needed.

Explicit (piped) outlets keep the existing validation: in-bounds, region ≠
footprint region, strictly below the spillpoint elevation (so a piped drain
can't route straight back over the saddle into the footprint).

## 5. Static flow-graph wiring (network inclusion)

The runoff (natural spillpoint) already links the NBS region to its downstream
region through the ordinary spillgraph — no extra wiring. For **piped** outlets,
reuse the culvert / `_nbs_links` connectivity already built: add a connector
link footprint → piped cell so `_build_dyn_networks` pulls the outlet's
downstream component into the same dynamic network.

## 6. Dynamic rate function (`dynNetworkRateFunction!`, NBS branch)

State per placement = per-layer storage `S` (mm). Each ODE step:

1. **Inflow claim.** The negative NBS region's accumulated routed inflow (+
   direct rain on footprint) feeds the **top** layer. One new hook: a map
   `nbs_region_id → element`; in inflow accumulation, redirect that region's
   inflow into the element instead of discarding it as "left domain." *This is
   the only genuinely new dynamic plumbing.*
2. **Cascade** (existing): layer `i` infiltrates to `i+1` at
   `Kinf·max(S−Smin,ε)^ninf`; bottom layer's infiltration → lost to ground.
3. **Overflow delivery** (reuse PR-3): each layer with `Kout>0` computes
   `Qout = compute_outflow(...)·1e-3` and delivers it position-exact to its
   resolved cell — runoff_cell or piped cell, same code path. Runoff water then
   flows downhill through the normal flow graph like any inflow.
4. `dS/dt = inflow − infiltration_out − overflow_out − evap` per layer.

## 7. Submergence — handled by the hierarchy

Because the NBS is now an interior trap with a spillpoint, submergence is a
sub/supertrap event, captured by the existing machinery rather than a bespoke
check:

- A neighbor trap filling to the shared saddle merges with the NBS into a
  supertrap; the sub/supertrap hierarchy already records this, and the dynamic
  network merge (`_merge_networks`) already couples them into one solved
  component.
- **When the NBS trap is a filled subtrap of an active supertrap** (its
  spillpoint elevation is submerged by the supertrap's water surface), the NBS
  is drowned: **stop infiltration** (no head differential to ground) and **stop
  any outlet whose cell is likewise submerged in that supertrap** (no gradient).
  This is a state read off the same fill state the solver already tracks
  (`cur_amounts` / filled-trap set), not a new geometric scan.
- The merged pool then fills and spills per ordinary supertrap dynamics toward
  the supertrap's own spillpoint.

Residual detail to confirm during implementation: the exact predicate for
"NBS submerged" in terms of the supertrap water level vs the NBS spillpoint
elevation, and how the layered state freezes/thaws across submerge/emerge — see
§7a.

## 7a. Submergence: complexity & handling

Complexity assessment: **medium**, and *not* dominated by the rerouting.

**Inherited for free (the rerouting worry).** "Water that spilled to a sibling
now routes up to the parent once submerged" is already implemented for ordinary
traps — `_compute_changetime_estimate` (`fill_sequence.jl:263-300`, the
`parentof` / `filled_traps[parent]` branch), `_reconcile_spillgraph!`
(`network_context.jl:274`), `_touch_networks!`. An interior NBS trap plugs
straight in; no new rerouting code.

**Key design choice — runoff outlet = natural spillpoint (collapses the hard
knot).** Route top-layer overflow to the natural spillpoint (the default). Then
the NBS's surface overflow spills over the *same terrain saddle the hierarchy
tracks*, so it is **indistinguishable from an ordinary trap spill** for routing —
sibling-spill and route-up-to-parent apply automatically. Any *other* runoff
target would create two competing overflow mechanisms (layer `Kout` to a cell
vs geometric spill over the saddle) needing a storage↔pond-volume reconciliation
to decide which fires first. Spillpoint runoff makes them the **same event**.
(A modeler-specified non-spillpoint / weighted runoff — open decision #2 — must
therefore accept this reconciliation cost; the single-spillpoint default avoids
it.)

**Genuinely new work (bounded, localized):**

1. **Submerge/emerge state handoff.** On drowning, the NBS's stored layer water
   folds into the parent pool volume *exactly once*; on emergence it returns.
   This rides the same trap-absorption seam that already produced two documented
   stale-state bugs (`staletime-bug-fix`, `stale-state0-network-absorption-bug`)
   — the highest-risk seam in the code. No new architecture, but demands careful
   projection to `cur_time` and targeted mass-conservation tests. **Main cost /
   risk.**
2. **Scheduler fill-time consistency.** `_compute_changetime_estimate` reads
   `trapvolumes − subvolumes` (terrain volume). The terrain volume below the
   spillpoint is correct (no carving), but the layer model governs the *rate*, so
   `getinflow` / `getsmax` / `getsmin` for an NBS trap must be fed a
   storage↔volume-consistent value or the scheduler mis-times the fill.
3. **Partial submergence is binary.** A lumped model has no elevation profile, so
   there is no consistent spatial fraction. Declare submergence binary
   (threshold at spillpoint elevation). Do not invent a partial-fraction physics.

**De-risking first cut.** When the NBS becomes a filled subtrap of an active
supertrap, **freeze** the layer dynamics (`dS/dt = 0`), fold current storage into
the parent pool, gate outlets by the submerged-cell rule; **thaw** on emergence.
The layer cascade then runs only while the NBS is surface-active (never
underwater), keeping the handoff a clean binary transition.

## 8. Mass ledger (must hold every step)

```
inflow_to_top = Σ_layers dS/dt·A/1000 + Σ overflow_delivered
                + bottom_infil_to_ground + evap
```

Every overflow term is delivered to exactly one cell; nothing created or
destroyed. Same invariant the culvert path already satisfies.

## 9. Delta vs current PR (nbs-dynamic-core)

**Keep:** `NBSLayer` / `NBSSystem` / `compute_outflow`, factories, mm↔m³
conversion, the layer cascade, position-exact delivery, `_nbs_links`
connectivity.

**Replace:**
- dig/carve → `-4` spillfield flag (interior region on real terrain);
- bolted-on spillpoint → the *natural* spillpoint (lowest saddle);
- supertrap cycle guard → not needed (natural spillpoint is cycle-safe;
  hierarchy handles merging/submergence);
- single-outlet semantics → per-layer outlet cells, unspecified backfilling to
  the natural spillpoint downstream cell.

**Net delta from current impl:** minus carving, minus spillpoint bolt-on; plus
the `-4` flag on real terrain and using the natural spillpoint as the runoff
outlet.

---

## Open decisions

1. **Default outlet policy for ≥3-layer systems** (mantilla:
   surface/soil/drainage). Proposal: topmost outflowing layer → runoff (natural
   spillpoint), bottommost outflowing layer → piped, middle layers percolate only
   (`Kout=0`). Matches the physical models. Alternative: fully explicit per-layer,
   no default. **TBD.**

2. **Runoff concentration vs terrain spread.** The interior model delivers a
   layer's overflow to one cell by default (the natural spillpoint). That's a
   *default*, not a constraint: a layer overflow may instead be delivered to a
   **weighted list of runoff cells** (`Vector{Pair{CartesianIndex{2},Float64}}`,
   weights summing to 1), each flowing downhill naturally; default = the single
   natural spillpoint at weight 1.
   - The single *static* spillpoint stays single regardless — it drives
     merge/submergence topology, which must be terrain-driven (correctness, not
     a limitation).
   - `NBSSystem` is **lumped** (one bucket per layer), so overflow is a single
     aggregate flux. Option 1's "spread per terrain" would distribute that one
     number by heuristic, not physics — the spatial info was already collapsed
     at intake. So terrain-spread is largely illusory on a lumped core.
   - Genuine limitation: a **large footprint spanning multiple downhill
     directions** (ridge pavement) funnels all overflow to one saddle. Fix by
     the weighted runoff list above, or by splitting the NBS into several
     placements (one per sub-catchment). Document the split guidance.
   **TBD:** ship the single-cell default now and add the weighted list only when
   a use-case needs it, or build the weighted list up front?

3. (space for further questions/comments during iteration)
