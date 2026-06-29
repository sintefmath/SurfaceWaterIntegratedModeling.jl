# Review of `integrate_networks_plan.md` against the code

The architecture is solid.  The spillgraph-exclusion approach (§4) and the
changetimeest ordering (§5/§6) hold up well against the code.

Two kinds of finding are tracked below:

- **Confirmed bugs** — defects in the existing code (constructor/solver) that
  must be fixed for the integration to be correct.  Numbered `B1`, `B2`, … so
  the numbering stays stable as items are added.
- **Observations** — gaps, risks, or clean-ups in the *plan*, grouped by
  severity.

---

## Confirmed bugs

### B1. Parent-trap infiltration double-counted (RESOLVED — Design A)

**Status: constructor fix implemented & tested; inflow side folded into the plan
(Design A — always subsume).**

**Implemented (this branch).**  The constructor half of Design A is in code:
`_subnetwork` now calls `_subsume_terminal_parent` (with `_is_descendant`) in
`src/dynamics/elements.jl`, collapsing a terminal *unfilled* parent and its
recorded full descendants into the single parent node and dropping the path
segments internal to its footprint.  Regression test:
`@testset "setup_network subsumes terminal parent (B1)"` in `test/dynamics_test.jl`
(lone-trap and external-approach branches, plus a scan over every multi-child
parent).  All dynamics tests pass; the only suite failure is the pre-existing
`MaskingTests`.  The inflow side (`_external_inflow`, gross composite =
Σ leaf-descendant `trap_inflow`) lives with the *not-yet-built* fill_sequence
integration and is documented in the plan, not implemented here.

**Note (found while implementing).**  A parent's footprint is *not* exactly the
union of its children's footprints — it has exclusive cells above the child
spillpoints (e.g. cell 2366 of parent 414, region 18, wetted only when 414
fills).  Design A is unaffected: every footprint cell's region is a leaf
descendant, so the gross-inflow sum still captures all composite runoff, and the
solver charges the whole composite footprint.  But this is why the inflow uses
the *gross* leaf-descendant sum rather than any "parent has no own infiltration"
shortcut.

**Corrected diagnosis.**  The original writeup said a subsumed parent's external
inflow is *undercounted* (≈ 0).  That was wrong on two counts:

- Subsumed children are *non-network* (not in `net_trap_set`), so they are not
  masked, and their `C→P` spillgraph edges *do* populate `trap_inflow[P]`
  (`flow.jl:53-58`; sibling→parent rule `spillgraph.jl:180-185`).
- The real defect is **infiltration double-counting**.  The spillgraph delivers
  `trap_inflow[P]` *net* of child infiltration (`signed_flow = inflow − Smax`),
  but the solver's parent geometry charges the *whole composite* footprint
  infiltration again (`wetted_infiltration` / `_footprint_infiltration`), with no
  `−Smin` term.  `fill_trap_until` gets this right
  (`infilfun = sum(fprint_infil) − Smin`, `fill_sequence.jl:470`); the
  DynNetwork solver omits it.

**Compounding: constructor inconsistency.**  `flow_path_from` subsumes children
when the parent is full, but records them as separate nodes when the parent is
the terminal *unfilled* trap.  A full child appearing as a separate node also
hits a latent hazard: when its runoff is below its infiltration it computes
`dV < 0` and would spuriously `:unspill`, though it cannot drain while submerged.

**Chosen fix — Design A (always subsume).**  A parent is always one node
subsuming all descendants; children are never separate nodes.  The node is fed
the *gross* composite inflow
`external_inflow[node] = Σ trap_inflow[leaf descendants]` (= `trap_inflow[node]`
for a leaf), which matches the solver's existing composite-footprint geometry —
so **no solver change** is required.  The only constructor change is to subsume a
terminal-unfilled parent's full descendants (the full-parent case is already
subsumed by `flow_path_from`).  Folded into plan §3 (Build), §4, §6; subsumed
descendants are also excluded from the standard changetime machinery (plan §7,
`net_covered_set`).

---

## Agreed & folded into the plan

These are resolved and now live in `integrate_networks_plan.md`; kept here only
as a record.

- **Relative vs. absolute time** — `solveDynNetwork!` returns *elapsed* time from
  the solve start; the plan now adds the absolute solve-start time
  (`+ cur_time`) before storing into `changetimeest` / `next_event` (plan §5, §6;
  `next_event` documented as absolute in §2).
- **Bounded integration / state over-advance** — `solveDynNetwork!` gains an
  optional `tmax`, and the plan adopts a predict-on-copy / commit-to-`T_commit`
  protocol so a context's `state` is never over-advanced past the committed event
  time (plan §5, §8).
- **Network-trap amounts come only from the solver** — network traps are excluded
  from the two `fill_trap_until` call sites (`_update_affected_amounts` and the
  end-of-period finalization loop); their volumes are read exclusively from the
  integrated `ctx.state` (plan §9).
- **Gate re-solves on actual inflow change** — instead of re-solving every network
  every event, the plan adopts a *touch* model: a context is advanced/re-predicted
  only when one of its traps fires or its external inflow changes.  A cached
  `extern_inflow` field lets an untouched network be advanced later in one bounded
  solve (plan §2, §3 touch criterion, §8).
- **State remap across a rebuild** — the plan spells out the remap: traps present
  before carry their committed volume by global index; newly-introduced traps are
  seeded from `cur_amounts[trap_ix]` (or the boundary value if just filled) (plan
  §3 "State remap across a rebuild").
- **Post-event boundary values** — already covered: plan §3 names all three
  (`:fill` → C, `:empty` → 0, `:unspill` → `prevfloat(C)`); confirm during
  implementation.

---

## Withdrawn

- **Rebuild trigger too narrow** — *incorrect; withdrawn.*  The claim was that a
  *non-network* fill could extend network membership.  But `_subnetwork`
  (`elements.jl:250-256`) appends the first unfilled trap downstream
  (`_unfilled_trap_at`) to `net.traps`, so the network's downstream frontier is
  itself a member.  Terrain flow is static, so the only traps whose fill-state can
  alter the trace are those *on* the trace — all members (full traps in `ftraps`,
  the single unfilled stopping point added explicitly).  A network whose frontier
  is *not* a member spills out of the domain and has nothing to extend into.  So
  membership only ever changes on a **member** fill, which the trigger
  `fill_updates ∩ net_trap_set` already catches.  Non-member fills change
  *inflow*, not membership — handled by re-solving, not rebuilding.

---

## Observations

None currently open — all have been folded into the plan, resolved, or withdrawn
(see the sections above).

---

## Bottom line

`B1` is resolved: corrected to an infiltration double-count (plus constructor
inconsistency), with Design A (always subsume) chosen and folded into the plan.
Everything raised in this review has now been agreed and folded into the plan,
resolved, or withdrawn — no outstanding items.
