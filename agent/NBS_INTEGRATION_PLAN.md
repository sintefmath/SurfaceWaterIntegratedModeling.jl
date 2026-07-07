# NBS integration plan

Integrate Nature-Based-Solution (NBS) systems into SWIM's dynamic-network
subsystem, so an `NBSPlacement` behaves as a rate-limited network element
alongside traps and rate-limited culverts.

Authoritative requirements + settled decisions live in
`agent/prompts/nbs_instructions.org` (questions Q1–Q10, all answered). This
document is the implementation plan derived from them. Where a decision is
referenced it is tagged `[Qn]`.

## 0. Settled decisions (recap)

- **Dedicated element type** `DynNBS` in the network — *not* a reused `DynTrap`
  `[Q1a]`. No topology-change (`:fill`/`:empty`/`:unspill`) triggers on its
  layer states `[Q1a]`.
- Each NBS with `n` layers contributes `n` state variables to the solver state
  vector; inflow enters the **topmost** layer only, sourced from the NBS
  "trap"'s spill region `[Q1]`.
- **All** per-layer outlets are routed dynamically; the downstream structure of
  every outlet cell must be pulled into the DynNetwork (like culvert endpoints)
  `[Q2]`.
- Lowest-layer infiltration leaves the modelled system (deep-ground loss),
  exactly like ordinary terrain infiltration `[Q2a]`.
- All solver arithmetic is in grid volume units = (cell area)×(height unit),
  assume metres if ambiguous, comment it `[Q3]`. NBS layer parameters stay in
  **mm** (SUrbArea shares the types); convert at the NBS/terrain boundary with a
  clear comment `[Q3a]`.
- The flattened "trap"'s geometric capacity is used **only** to attribute a
  spill region + spillpoint in static analysis; it is otherwise ignored by the
  solver `[Q3b]`.
- Terrain infiltration is zeroed over NBS footprint cells; the NBS layer model
  is the sole water authority there. Ordinary `DynTrap` `wetted_infiltration`
  is **unchanged** `[Q8]`,`[Q8a]`.
- Layer area `A` is overwritten from the footprint (footprint = source of
  truth), all layers share it `[Q6]`,`[Q9]`.
- Unspecified outlet sentinel is `CartesianIndex(0,0)`; an unspecified
  non-lowermost outlet defaults to the bottom outlet + warning; an unspecified
  lowermost outlet defaults to the natural spillpoint + warning `[Q7]`,`[Q7a]`.
- Static sequencing: flatten → provisional spillpoint (bottom outlet) → run
  regions → validate/backfill outlets → warn if a spillpoint was moved
  `[Q4]`,`[Q10]`. Dig depth: one common fixed depth for all NBS footprints, deep
  enough to guarantee a nonzero-volume trap `[Q10]`.

## 1. Data structures

### 1.1 `DynNBS` (new, `src/dynamics/elements.jl`, next to `DynCulvert`)

```julia
struct DynNBS <: DynObject
    trap_ix::Int                        # artificial NBS trap index in the spillanalysis structure
    placement_ix::Int                   # index into tstruct.nbs (the source NBSPlacement)
    outlets::Vector{CartesianIndex{2}}  # resolved outlet cell per layer (bottom→top order matches system.layers)
    # solver bookkeeping filled at rate-param build time (not stored here; see §3.2)
end
```

Rationale for the fields: `trap_ix` ties the NBS to its artificial trap (state
carryover, covered-set bookkeeping, spill-region lookup `[Q5]`);
`placement_ix` recovers the `NBSSystem` (layer parameters) and footprint from
`tstruct.nbs`; `outlets` are the *resolved* cells (post-backfill), one per layer.
The per-layer→state-index map is derived per solve (§3.2), not stored on the
element, mirroring how culvert routing data is rebuilt per solve rather than
cached on `DynCulvert`.

### 1.2 `TrapStructure` new field `[Q1 preamble / §22 of the org]`

Add one field:

```julia
nbs::Vector{NBSPlacement}   # NBS installations attributed to this terrain (empty if none)
```

**Include-order problem.** `NBSPlacement` is defined in
`src/dynamics/nbs_elements.jl`, included at line 24 of
`SurfaceWaterIntegratedModeling.jl`, but `TrapStructure.jl` is included at line
12 — so the struct field cannot name `NBSPlacement` as written. Resolution
(primary): move the two *type-definition* includes earlier:

```
include("spillpoints.jl")
include("trapvolumes.jl")
include("dynamics/elements.jl")       # DynObject + Dyn* types   (moved up)
include("dynamics/nbs_elements.jl")   # NBSLayer/NBSSystem/NBSPlacement (moved up)
include("spillanalysis.jl")
include("TrapStructure.jl")           # now sees NBSPlacement
...
```

`elements.jl`/`nbs_elements.jl` reference `TrapStructure` only inside function
bodies (untyped `tstruct` args), so moving them ahead of `TrapStructure.jl` is
safe *provided* a full test run confirms no include-time dependency on
later-defined names. **Fallback** if the reorder surfaces a cycle: type the
field `Vector{NBSPlacement}` still requires the type early, so instead declare
an empty abstract supertype `abstract type AbstractNBSPlacement end` in an
early file and make `NBSPlacement <: AbstractNBSPlacement`, then type the field
`Vector{<:AbstractNBSPlacement}` via a parametric field. The reorder is
cleaner; prefer it.

The `TrapStructure` constructor call in `spillanalysis.jl:146` gains the new
positional argument (append `nbs` last; update the docstring in
`TrapStructure.jl` too — remember `no-new-trapstructure-fields` is explicitly
waived here per the org §22).

## 2. Static analysis (`spillanalysis` path)

New file `src/nbs_placement.jl` (included before `spillanalysis.jl`) holds the
static helpers; `spillanalysis` calls them. Signature change:

```julia
function spillanalysis(grid; ...,
                       nbs = nothing)      # untyped default, resolved in body (see below)
```

`nbs` is untyped in the signature and resolved to `NBSPlacement[]` in the body —
same pattern `fill_sequence` uses for `culverts` (`fill_sequence.jl:40-43`),
because `NBSPlacement` is defined in a later include *if* the §1.2 reorder is
not adopted; keep the body-resolution regardless for symmetry.

Order of operations inside `spillanalysis`:

1. **Validate + normalise placements** (`_prepare_nbs!`):
   - Overwrite every layer's `A` = `length(footprint) * cell_area` `[Q6,Q9]`
     (cell_area = 1.0 for now, `@@@` TODO to read real resolution — same
     hardcode the `DynCulvert` constructor carries).
   - Sanity-check the footprint (non-empty, in-bounds, contiguous — warn if not
     contiguous).
2. **Flatten terrain** (`_dig_nbs_traps!` on `gridcpy`, before `spillfield`):
   set every footprint cell to a common fixed low elevation `NBS_DIG_LEVEL`
   `[Q10]`. Choose it relative to the global grid minimum, e.g.
   `minimum(gridcpy) - NBS_DIG_DROP` with `NBS_DIG_DROP` a fixed constant (e.g.
   1 height-unit; note that for visualization purposes it should not be very far
   from the regular grid elevation values) so a nonzero-volume trap is
   guaranteed and the cells are visually recognisable. Comment the unit
   assumption `[Q3]`.
   - Set `NBS_DIG_LEVEL = minimum(gridcpy) - NBS_DIG_DROP` so the footprint
     becomes the strictly-lowest terrain in the grid; a *uniform* flat bottom is
     fine. The flat-trap elimination in `spillregions`
     (`spillanalysis.jl:103-105`) is gated on **zero volume**, not on flatness —
     digging below the global minimum guarantees the trap's volume is positive
     (depth ≥ `NBS_DIG_DROP` up to the surrounding saddle), so the NBS trap
     always survives. `NBS_DIG_DROP` need only be a small positive value (≈1
     height-unit) for both a positive volume and visual recognisability.
     Validate in the §6.1 static test (assert the footprint is exactly one
     surviving trap).
3. Set a **provisional spillpoint**: nothing to do pre-regions except record
   which footprint each NBS owns; the actual spillpoint is a property of the
   trap that `spillpoints`/`sshierarchy` will compute. The provisional target is
   the **lowermost layer's** outlet cell `[Q2, §24]` — but since spillpoints are
   derived, we instead *post-process* in step 5.
4. Run the existing pipeline unchanged (`spillfield`, `spillregions`,
   `spillpoints`, `sshierarchy!`, `trapvolumes`, footprints, …) on `gridcpy`.
5. **Resolve NBS traps + outlets** (`_resolve_nbs!`, after regions/spillpoints
   exist):
   - For each placement, look up its trap index from `regions[footprint[1]]`
     `[Q5]` (all footprint cells share one region after digging; assert this).
   - For each **specified** outlet: error if it lies inside the NBS trap's spill
     region (`regions[outlet] == nbs_trap_region`) `[§25]`; else confirm it sits
     at a **strictly lower** elevation than the NBS trap's natural spillpoint —
     error otherwise `[§15]`. Rationale (user): the elevation test is not a
     literal "graph-downstream" check but a **cycle guard** — an outlet strictly
     below the natural spillpoint cannot route water back up into the trap, so
     no feedback loop can form. (Strict `<`, not `≤`: an equal-elevation outlet
     gives no such guarantee.)
   - For each **unspecified** outlet `CartesianIndex(0,0)` `[Q7]`: the lowermost
     defaults to the trap's natural spillpoint (`spoints[trap].current_region_cell`
     / spillpoint cell); non-lowermost default to the bottom outlet; **warn** in
     both cases `[Q7a]`.
   - **Supertrap cycle guard**: no outlet may share a supertrap with the NBS
     trap (`intersect(supertraps_of[nbs_trap], supertraps_of[regions[outlet]])`
     non-empty → error).  If they shared an enclosing supertrap, filling it
     would later merge the two and route discharge back into the NBS — a cycle
     the strict-below-spillpoint elevation test alone does not catch (that test
     only compares to the NBS's own spillpoint, not to a shared parent's fill
     level).
   - Backfill of an unspecified outlet uses the spillpoint's **downstream** cell
     (`downstream_region_cell`, where the trap overflows *to*), never
     `current_region_cell` (which is inside the NBS region — using it would make
     discharge re-enter the NBS itself).
   - Write the resolved outlets back into the `NBSPlacement.outlets` vector
     (mutating its contents — the field is a `Vector`, contents are mutable).
6. Store the (mutated) `nbs` vector in the returned `TrapStructure`.

No terrain-infiltration change happens here — infiltration is a `fill_sequence`
input, handled in §3.

## 3. Dynamic solver (`fill_sequence` + `dynamics/`)

This is the substantial part and can be split into 3a/3b/3c sub-PRs.

### 3.1 Auto-trigger + network seeding

- `fill_sequence` already accepts `dyn_traps`/`culverts`. Add: NBS trap indices
  (from `tstruct.nbs` via §2.5 lookup) are unioned into the effective dynamic
  seed set automatically `[Q5]`, so the presence of any NBS triggers DynNetwork
  construction even with `dyn_traps === nothing`.
- `_dyn_seeds` (`network_context.jl:137`) additionally emits, per NBS: a
  representative footprint cell **and** every resolved outlet cell (so the
  downstream structure of each outlet is traced into the network, exactly as
  culvert inlet/outlet cells are) `[Q2]`.
- `setup_network`/`_subnetwork`/`_expand_with_culverts` (`elements.jl:187+`):
  when a traced trap index is an NBS trap, emit a `DynNBS` node instead of a
  `DynTrap`, and expand the network from each of its outlet cells the same way
  `_expand_with_culverts` expands from culvert endpoints. NBS nodes may merge
  otherwise-disjoint components (same as culverts). The incremental `SubnetCache`
  occupied-cell logic (`_subnet_cells`, `elements.jl:221`) must count NBS
  footprint cells as occupied.

### 3.2 State-vector layout + rate params

- The per-context `state` / solver `V` keeps one entry per `net.traps` (trap
  volumes, m³) and **appends** NBS layer states after them: for each `DynNBS` in
  the net (stable order), a contiguous block of `nlayers` entries, stored in the
  same volume unit (m³) as the rest of `V` for solver-uniform tolerances.
- The NBS layer block is initialised (at every context build/rebuild) from the
  **persistent NBS-state store** (§3.5), *not* zeroed — that store is the single
  source of truth for NBS layer volumes across events and weather periods, the
  same role `cur_amounts` plays for trap volumes.
- `_build_rate_params` (`networksolver.jl:611`) precomputes, per `DynNBS`:
  - `state_base` (offset of its layer block in `V`) and `nlayers`;
  - per-layer `(Kout, nout, Smax_vol, Kinf, ninf, Smin_vol, A)` with the mm→vol
    conversions applied once (`Smax_vol = Smax_mm * 1e-3 * A`, etc.) `[Q3a]`;
  - per-layer outlet → downstream injection target (which flow-path index +
    cell position, or trap), resolved like `_build_culvert_plan`
    (`networksolver.jl:250`);
  - the NBS node's external inflow (its spill-region leaves) feeding layer 1.
- `DynNetworkRateParams` gains an optional `nbsplan` field (nothing when the net
  has no NBS), mirroring the optional `cvplan`.

### 3.3 Rate function

Extend `dynNetworkRateFunction!` / `_routed_inflow` (`networksolver.jl:667-714`)
so NBS layer outflow is injected into routing and NBS layer states evolve:

- Before/within `_route_flow`, compute each NBS layer's **outflow** from its
  current state `V[state_base+l]`:
  `out_mm_m2 = compute_outflow(Kout, nout, Smax_mm, S_mm)` where
  `S_mm = V_layer_m3 * 1e3 / A`, then `out_m3 = out_mm_m2 * 1e-3` `[Q3a]`.
  Inject `out_m3` into the downstream structure of that layer's outlet cell —
  the same delivery mechanism as `culvert_actual_delivered`
  (`_route_flow`, mass-conservation rule in AGENTS.md: delivered == drawn).
- Compute each layer's **infiltration**
  `inf_mm_m2 = compute_outflow(Kinf, ninf, Smin_mm, S_mm)`, `inf_m3` likewise.
  Route `inf_m3` into layer `l+1`'s state; the **bottom** layer's infiltration
  leaves the system (deep-ground loss, not injected anywhere) `[Q2a]`.
- Layer-state derivative (m³/time):
  `dV[state_base+l] = inflow_l - out_m3(l) - inf_m3(l)` where
  `inflow_l = (l==1 ? nbs_external_inflow : inf_m3(l-1))` `[Q1]`. No `/A` because
  we work in volumes, not depths — the `A` division that `NBSNetworkRateFunction!`
  does is folded into the mm↔m³ conversion. Evapotranspiration term omitted
  (deferred) — assert/ignore `EVCoeff`.
- Ordinary trap entries are unchanged. Terrain infiltration over NBS footprint
  cells is zeroed (§3.4), so no double count `[Q8]`.

### 3.4 Zeroing terrain infiltration under NBS footprints `[Q8]`

In `fill_sequence`, after building the infiltration map (`fill_sequence.jl:49`),
zero it on every NBS footprint cell:
`infiltration[cell] = 0` for `cell in placement.footprint` for each NBS. This is
NBS-specific; ordinary trap footprints keep their `wetted_infiltration` `[Q8a]`.
Do this on a copy so a caller-supplied map isn't mutated.

### 3.5 Events + steady state `[Q1a]`

- NBS layer states get **no** `:fill`/`:empty`/`:unspill` conditions — skip them
  in `_event_conditions` (`networksolver.jl:774`); the event system indexes only
  the `nt` trap entries. Downstream traps fed by NBS outflow keep their normal
  events.
- The `evolving` set and `all(abs(du) <= abstol)` steady-state checks
  (`networksolver.jl:1090-1098`) operate over the whole `V`, so NBS layer states
  are naturally included — a net with only NBS still evolving will not be
  declared steady prematurely. Confirm `_validate_network` and the FULL/EMPTY
  fast paths iterate `1:nt` (trap entries) only, not the NBS block.
### 3.6 NBS-state persistence across events and weather periods `[required]`

NBS layer volumes **must persist across weather events** (user requirement) — an
NBS half-saturated at the end of one rain period must start the next period from
that same saturation, exactly as trap volumes carry via `cur_amounts`.

`cur_amounts` cannot hold this: it is one `FilledAmount` per trap, whereas an
NBS has `nlayers` states. So add a dedicated store:

```julia
nbs_state::Dict{Int, Vector{Float64}}   # placement_ix -> per-layer volume (m³)
```

- **Owned by `fill_sequence`**, threaded through the weather-event loop
  (`fill_sequence.jl:66`) alongside `cur_amounts` — *not* a `TrapStructure`
  field (only the placements themselves go on the struct; per
  `no-new-trapstructure-fields`, evolving state stays in `fill_sequence`).
- **Initialised** to `zeros(nlayers)` per placement at the start of
  `fill_sequence` (all layers empty at t=0). *(Future option: let
  `NBSPlacement` carry an initial saturation; zero for now — note it.)*
- **Read** at every context build/rebuild: `_make_context`
  (`network_context.jl:175`) fills each `DynNBS`'s layer block from
  `nbs_state[placement_ix]` (§3.2).
- **Written back** whenever a context's evolved state is committed — the same
  points that update `cur_amounts` from network state
  (`_finalize_networks!`, `fill_sequence.jl:226`, and any intra-period commit
  before a `_build_dyn_networks` rebuild). Extract each `DynNBS`'s layer block
  from `ctx.state` and store it. This makes `nbs_state` the single source of
  truth, so intra-period rebuilds and cross-period restarts both read consistent
  values.
- Because the network solve runs to the weather-period boundary (`end_time`),
  NBS states evolve continuously up to each boundary and resume from there — no
  gap, no reset.

Persisting NBS state into per-event `SpillEvent` *snapshots* (for replay/plotting
of intermediate NBS saturation) is a **separate, later** concern; the
requirement here is dynamical continuity of the simulation, which `nbs_state`
delivers. Extending `SpillEvent` can follow if the animation work needs it.

## 4. Acyclicity assumption

The DynNetwork requires water to flow strictly downstream (no cycles); culverts
already assume this. For NBS this is **enforced, not merely assumed**: the §2.5
outlet check requires every outlet to sit *strictly below* its NBS trap's
natural spillpoint (the elevation at which water leaves that trap). Each routing
hop is therefore strictly decreasing in spillpoint elevation, so no chain of NBS
outlets — including two NBS discharging toward each other — can close into a
cycle. No separate cycle-detection pass is needed; the elevation validation is
the guard. (Should a future feature relax the elevation rule, reinstate an
explicit acyclicity check, mirroring the deferred reverse-culvert case in
AGENTS.md.)

## 5. Files touched

| File | Change |
|------|--------|
| `src/SurfaceWaterIntegratedModeling.jl` | reorder type-definition includes (§1.2); add `nbs_placement.jl` |
| `src/TrapStructure.jl` | add `nbs` field + docstring |
| `src/nbs_placement.jl` (new) | `_prepare_nbs!`, `_dig_nbs_traps!`, `_resolve_nbs!`, constants |
| `src/spillanalysis.jl` | `nbs` kwarg, call the helpers, pass field to constructor |
| `src/dynamics/elements.jl` | `DynNBS` struct + docstring; network build recognises NBS traps |
| `src/dynamics/nbs_rate.jl` (new) | NBS per-layer outflow/infiltration + injection (analogue of `culvert_rate.jl`) |
| `src/dynamics/networksolver.jl` | `nbsplan` in rate params; NBS branch in rate fn + steady-state; skip in events |
| `src/dynamics/network_context.jl` | `_dyn_seeds` emits NBS seeds/outlets; state layout; `_make_context` reads `nbs_state` (§3.6) |
| `src/fill_sequence/fill_sequence.jl` | auto-union NBS traps; zero footprint infiltration; own + thread `nbs_state`, write it back on commit/finalize (§3.6) |
| `test/` | new `NBS` test set (see §6) |

## 6. Tests

Deterministic, real-terrain, per AGENTS.md (no mocked artifacts):

1. **Static**: `spillanalysis` with a `puddle` placement on a small real grid —
   assert the footprint became one trap, spillpoint resolved to the outlet,
   `tstruct.nbs` populated, `A` overwritten, footprint infiltration zeroing hook
   fires.
2. **Outlet validation**: outlet inside spill region → error; unspecified outlet
   → backfilled + warning.
3. **Mass conservation** (the load-bearing test): a single `mantillaRRmodel`
   NBS in a small catchment under a constant-rain weather event; assert
   `inflow_total ≈ Δstorage(traps) + Δstorage(NBS layers) + infiltration_lost +
   out_of_domain` to tight tolerance. This is the primary guard `[AGENTS.md
   mass conservation]`.
4. **Degenerate parity**: a `puddle` with `Kinf=0`, single layer, `Smax=0`,
   large `Kout` should behave ≈ a pass-through (spill ≈ inflow), matching the
   no-NBS routing within tolerance.
5. Run the existing `Sequencing` + `Trapping structure` sets (fill_sequence
   touched) and `SpillfieldTests`/`SpillregionTests` (spillanalysis touched),
   per AGENTS.md §Making changes.
6. **Cross-period persistence** (§3.6): two back-to-back weather events; assert
   the NBS layer volumes at the start of event 2 equal those at the end of event
   1 (no reset to zero), and that a single 2-event run gives the same NBS
   trajectory as one continuous event of equal total duration + rain.

## 7. Suggested PR sequencing

- **PR-1 (static)**: §1.2 field + include reorder, §2 (`nbs_placement.jl`,
  `spillanalysis` kwarg), test §6.1–6.2. Self-contained, no solver changes.
- **PR-2 (dynamic core)**: §1.1 `DynNBS`, §3.1–3.4 (seeding, state, rate,
  infiltration zeroing) and §3.6 (`nbs_state` persistence — required for
  correct multi-period runs, and coupled to the §3.2 state layout, so it lands
  here not later), test §6.3–6.4, §6.6.
- **PR-3 (events/steady + hardening)**: §3.5, §4 acyclicity guard, full test
  §6.5.

Each PR branches from `main`, keeps public keyword interfaces intact (the only
new public surface is the `nbs` kwarg on `spillanalysis` and the `DynNBS` /
`NBSPlacement` exports), and runs the mandated test sets before merge.
