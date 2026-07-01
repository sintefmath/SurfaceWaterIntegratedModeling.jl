# Review of `prompt_modification.org` — solveDynNetwork redesign

## What is clearly right

- **Exact equality for FULL/EMPTY** (V == C, V == 0) makes state classification
  unambiguous and moves the numerical responsibility to the caller, where it
  belongs — not inside the integrator.
- **Validating the contract at entry** is much better than trying to recover
  from bad states mid-solve. The current code works around bad states with
  tolerance tricks; rejecting them upfront is cleaner.
- **FULL→TRANSITORY check before integration** (the `:unspill` replacement):
  if any FULL trap already has negative net inflow at t=0, report it immediately
  without running the ODE. This is a clean replacement for the current
  `prevfloat(cap)` sentinel approach.
- **EMPTY→TRANSITORY is a non-event for leaf traps**: correct, no topology
  change occurs.
- **Reduced-dimensional ODE** is the natural consequence of the design: if only
  TRANSITORY traps go into the ODE state vector, non-evolving traps (FULL or
  EMPTY) simply are not there, and the drift problem disappears entirely without
  any freeze or clamping logic.

---

## Clarifications needed

### 1. Simultaneous fills — who detects them?

The ODE callback fires for one trap at a time. If traps A and B both hit
capacity in the same ODE step, `solveDynNetwork` returns `:fill` for A, with
`state[B] ≈ C_B` (within ODE error, but not exactly C_B). The caller receives
a state where B is not yet FULL by the exact-equality definition — but B will
fill "instantaneously" on the very next call.

Two options:

**(a)** The caller, after receiving any `:fill`, scans `state` for any trap
where `V >= C` and clamps those to exactly C, treating them as FULL and
updating topology accordingly. This makes the caller responsible for
simultaneous detection. `solveDynNetwork` stays simple.

**(b)** `solveDynNetwork` scans after the callback fires and returns ALL
simultaneously-filled traps in a single response (e.g. a vector of trap
indices rather than one). More self-contained, but changes the return type.

Which do you intend?

### 2. What "has a spill path" means for validation

A FULL trap must have a designated spill path out of it. In the current
`DynFlowPath` structure, paths have a `target_trap` field (where they lead),
but not a `source_trap`. The routing code sweeps in topological order and
routes spill from full traps, but the path structure does not explicitly record
which trap "owns" a given path as its outlet.

How should the validator determine that a specific FULL trap has a spill path?
Does the network structure need a new field (e.g. `source_trap` on
`DynFlowPath`), or is there an existing way to query this?

### 3. TRANSITORY→EMPTY — which traps get the empty callback?

The spec says TRANSITORY→EMPTY is only relevant for traps **with subtraps**
(because exposing subtraps changes topology). Two options:

**(a)** Only register the `:empty` callback for evolving traps that have
subtraps in the `TrapStructure`. Leaf traps get no `:empty` callback.

**(b)** Keep the callback for all evolving traps but filter in the `affect!`
handler — ignore the event for leaf traps.

Option (a) is semantically cleaner. Option (b) is easier to implement
incrementally. Which do you prefer?

### 4. Exactly what the caller does after each event

The prompt says: *"the caller should keep track of which traps are currently
full/empty and reflect this when sending a network along with its state."*

Does this mean:

**(a)** After a `:fill` event for trap T, the caller:
  1. Sets `state[T] = C[T]` exactly.
  2. Calls `setup_network` (or equivalent) to rebuild the network topology
     (opening the downstream spill path from T).
  3. Passes the rebuilt network and updated state to the next `solveDynNetwork`
     call.

**(b)** The same network object is reused with an updated state vector, without
rebuilding the topology. The network is only rebuilt at coarser events.

The current `_run_cascade` in the verification script reuses the same network
across all events. If the intent is (a) — topology rebuild after every fill —
the verification script and any other outer loop would need updating. This also
determines what "valid network" means: option (a) requires the network to
always reflect the current full/transitory/empty partition, while option (b)
allows stale topology.

### 5. Reduced-dimensional ODE — internal mapping

If the ODE state vector contains only TRANSITORY traps (not FULL or EMPTY), the
solver needs to:

- Build a mapping from the full `nt`-length external state to the smaller ODE
  state at the start of each solve.
- Reconstruct the full state (inserting FULL traps at C and EMPTY traps at 0)
  before returning.

Is this the intended architecture, or should the ODE continue to use the full
`nt`-length vector with FULL/EMPTY traps frozen at their boundary values?

---

## Potential inconsistency with `fill_sequence`

The main structural difference: `fill_sequence` tracks the sub/supertrap
hierarchy explicitly. When a trap fills, its subtraps are exposed and become
the new active leaf nodes. `solveDynNetwork` currently operates only at the
leaf-trap level (it never models subtrap exposure).

The new spec mentions TRANSITORY→EMPTY for traps **with subtraps**, which
implies `solveDynNetwork` will eventually need to detect when a supertrap
empties and expose its subtraps. Is this in scope for the current redesign, or
deferred to a later increment?

The EMPTY→TRANSITORY direction has the same subtrap implication: when a
supertrap empties and exposes its subtraps, the subtraps become the new leaf
nodes, and the outer loop must rebuild the network. This is consistent with the
"topology rebuild after each event" approach — just worth being explicit that
the outer loop (not `solveDynNetwork`) is responsible for this rebuild.

---

## Summary of open questions

| # | Question | Options |
|---|----------|---------|
| 1 | Who detects simultaneous fills? | (a) caller scans after each event, (b) solver returns all fills at once |
| 2 | How does validation check "has a spill path"? | Needs a query mechanism or new field |
| 3 | Which traps get the `:empty` callback? | (a) only subtrapped traps, (b) all but filter in affect! |
| 4 | Does each fill event trigger a topology rebuild? | (a) yes, always, (b) no, reuse network |
| 5 | ODE dimensionality — full nt or reduced? | Reduced is cleaner; full with freeze is easier |
| 6 | Is subtrap exposure in scope for this redesign? | In scope or deferred? |

---

## Response to answers in `prompt_modification.org` (2026-06-25)

### Q1 — Simultaneous fills: agreed, with one precision

The "let it fire instantly on the next call" approach is clean and avoids
changing the return type.  One subtlety the answer doesn't address: after the
ODE fires for trap A at time t_A, trap B's state could be `V_B > C_B`
(ODE overshoot past capacity, not just stopped just short).  `V_B > C_B` is
neither FULL (not exactly `C_B`) nor a valid TRANSITORY state (`V > C` is
unphysical).  The caller will need to handle this by clamping: **if `V >= C`,
clamp to exactly `C` and treat as FULL**.  This should be documented as part
of the caller contract.  Concretely: after every `solveDynNetwork` return, scan
`state[i] >= C[i]` (not just `== C[i]`) and clamp those, then process them as
fills before the next call.

### Q2 — Spill path validation via `DynTrap.spill_path`: agreed

`DynTrap.spill_path > 0` already exists as a field (`elements.jl:60`).  The
validation check is simply: for every FULL trap in the network, assert
`trap.spill_path > 0`.  No new fields needed.  The converse check (no
TRANSITORY trap has `spill_path > 0`) is already structurally guaranteed by
how `_build_network` constructs paths — the last trap in the chain gets
`spill_path = 0` — but it is worth asserting explicitly anyway.

### Q3 — TRANSITORY→EMPTY for parent traps only: agreed, with a cleaner test

`trap_ix > numregions(tstruct)` works, but `!isempty(subtrapsof(tstruct,
trap.trap_ix))` is more self-documenting and uses the exported API directly.
These are equivalent (`numregions` counts lowest-level traps, and only
lowest-level traps have no subtraps), but the `subtrapsof` form makes the
intent clear.  Either is fine; preference goes to `subtrapsof`.

### Q4 — Topology rebuild after every event: agreed

Rebuild on every event is correct.  Computational cost: `setup_network` is
O(path length from start cell to domain edge), which is small compared to the
ODE solve.  No concern there.

One clarification request: after a `:fill` event for trap T, the caller sets
`state[T] = C[T]` exactly and adds T to `full_traps` before calling
`setup_network`.  After an `:unspill` (FULL→TRANSITORY) event for trap T,
the caller removes T from `full_traps` and rebuilds.  Is the reverse also
correct — the caller does NOT need to adjust `state[T]` for the unspill case,
since by definition `state[T]` returned from `solveDynNetwork` will already
have been updated by the FULL→TRANSITORY check and will be exactly `C[T]` (we
detected it at t=0)?  Or will `state[T]` in the returned result be unchanged
(still `C[T]`) and the caller is expected to leave it as-is for the next call,
now as a TRANSITORY trap (since `C[T]` is no longer exactly `C` for the
new network geometry)?

Actually — this is fine: the old network had T as FULL with `V = C`.  The new
network has T as TRANSITORY.  The capacity in the new network is the same C.
So passing `V = C` for T into the new call makes T FULL again according to the
new network's classification, which would immediately trigger the FULL→TRANSITORY
check at t=0 again (infinite loop).  **The caller must set `state[T] < C[T]`
after an unspill event** — presumably to the value returned by
`solveDynNetwork` minus a small decrement, or the rate function's prediction
over a tiny dt.  This needs to be specified.

### Q5 — Reduced-dimensional ODE: FULL traps can be safely excluded

The key insight: **with constant inflows and no culverts, FULL traps cannot
transition to TRANSITORY during integration**.  Here is why:

A FULL trap's net rate is:

    dV/dt = inflow + upstream_spill − infiltration_full − culvert_drain − spill_out = 0

`spill_out` is the adjustable term that makes this exactly zero.  The routing
pattern is fixed for the duration of one `solveDynNetwork` call (topology does
not change mid-call), so `inflow + upstream_spill` is constant.  The trap is
fully submerged, so `infiltration_full` is constant.  Therefore `dV/dt = 0`
throughout, and the trap stays full.

**With culverts**: a culvert's drain rate depends on the outlet trap's water
level.  If the outlet trap is TRANSITORY (filling up), its level rises, which
*reduces* the culvert head differential, which *reduces* `culvert_drain`.  A
decreasing `culvert_drain` makes the net rate less negative — the FULL trap
is *more* likely to stay full during integration, not less.  FULL→TRANSITORY
can therefore only be triggered at t=0 (when culvert drain exceeds inflow at
the start), never during integration.

**Conclusion**: FULL traps are safe to exclude from the ODE state vector in
all current and planned configurations.  The t=0 FULL→TRANSITORY check is
sufficient.  This should be stated clearly in a code comment when implemented.

### Q6 — fill_sequence wrinkle: I believe the current behavior is correct

The scenario: parent trap is EMPTY, both children are FULL.

`_unfilled_trap_at` looks up `supertraps_of[cur_reg]`, filters out entries in
`full_traps`, and returns the minimum.  Since children are in `full_traps` and
the parent is not, it returns the parent.  The network is built with the EMPTY
parent as the leaf node — which is **the correct behavior**:

- Children are full (water is at their spillpoints).
- The parent represents the additional volume *above* both children's
  spillpoints.  With V=0, the water level is exactly at that merger elevation.
- As inflow continues, water accumulates in the parent (EMPTY → TRANSITORY →
  FULL), not in the children (they are already at capacity).

Removing the children from `full_traps` at this point would be wrong: it would
expose them as leaf nodes and try to fill them again, ignoring that they are
already physically full.

**However**, there is one genuine edge case to check: if `flow_path_from` is
called with a start cell that lies *within* one of the full children's
footprints, will `_unfilled_trap_at` correctly skip the child and return the
parent?  Looking at the code (`elements.jl:263−268`):

```julia
unfilled = filter(x -> x ∉ full_traps, tstruct.supertraps_of[cur_reg])
isempty(unfilled) ? 0 : minimum(unfilled)
```

`supertraps_of[cur_reg]` lists the lowest-level trap for this region plus all
its supertraps in order.  The child trap *is* the lowest-level trap for its
region, so it appears first.  Since the child is in `full_traps`, it is
filtered out.  The parent appears second and is not in `full_traps`, so it is
returned.  This is exactly right.

**I believe no change to `flow_path_from` or `_subnetwork` is needed for this
case.**  The "wrinkle" resolves correctly in the existing code.  Worth adding a
comment to `_unfilled_trap_at` explaining this invariant, but no logic change
required.  Please verify this reading against a known test case before
implementing.

---

## Remaining open question

After reviewing all answers, one new question has emerged (see Q4 above):

**After an `:unspill` (FULL→TRANSITORY) event**, what value should the caller
pass for `state[T]` on the next call?  The returned state has `state[T] = C[T]`
(the trap was full at t=0 when the event was detected), but passing `C[T]`
again would classify T as FULL, triggering the same event infinitely.

Options:
- **(a)** `solveDynNetwork` returns `state[T] = C[T] - ε` for the unspill case
  (a tiny decrement to clearly mark T as TRANSITORY), so the caller can use the
  returned state directly.
- **(b)** The caller is responsible for decrementing `state[T]` after an
  unspill event before the next call — but by how much?  An arbitrary ε is
  fragile.
- **(c)** The unspill event causes the caller to rebuild the network WITHOUT T
  in `full_traps`, and the new network passes `state[T] = C[T]` — but now C[T]
  in the new network is different (the new geometry has a different capacity for
  the active trap), so T is no longer exactly at capacity in the new geometry.

Option (c) seems the most principled: after an unspill, the network is rebuilt
with T removed from `full_traps`.  In the new network, the relevant trap for
T's water might be represented differently (a subtrap, or T itself but with a
different capacity).  The caller passes the old volume for T, which is `C[T]`
— but in the new network, T's capacity may differ, so T is naturally
TRANSITORY.  This needs confirmation.
