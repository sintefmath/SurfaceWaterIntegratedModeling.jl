# Culvert support — work in progress (note to self)

## Status

**Task 1 (network construction) — partially done.** `setup_network` now takes an
optional `culverts` kwarg, expands the network around touched culverts, dedups
traps, groups paths+traps into components (with culvert links), and assigns
culverts to their owning path/trap. All existing dynamics tests pass plus new
unit tests for `_components` (cross-component culvert link, lone trap survival).

**Task 1 — DONE (construction + tests + dead-code cleanup).** Added the testset
`"setup_network culverts on mini.txt"` in `test/dynamics_test.jl` (46 assertions,
all passing).  It covers the five scenarios below plus three edge cases (multiple
culverts in one network; fix-point/chained inclusion where one culvert's expansion
pulls in another; terrain-*inlet* expansion mirroring the terrain-outlet case).
`valid_network` was extended to bounds-check culvert references in paths/traps, a
`cvlt(inlet, outlet)` helper builds a `DynCulvert` with throwaway hydraulic
parameters, and `_culvert_owners` got a unit test for the trap-over-path tie-break.

**Dead-code cleanup.** Removed vestigial culvert handling left over from the earlier
"culverts ride on paths" design (superseded by end-stage cell-ownership in
`_culvert_owners`/`_build_component`): `_resolve_cell_overlaps!` lost its
`all_culverts` parameter and truncation-drop filter, and `_dedup_traps` lost its
culvert-list union.  Both ran only on always-empty lists.  Pure no-op refactor —
full dynamics suite unchanged.  (`_combine_networks` was likewise simplified and
renamed `_combine_subnets`: it now flattens paths/traps only, since subnets are
always culvert-free and culvert ownership is resolved later by cell.)

1. **Inclusion / assignment** — `:long` network (start `(7,119)`); culvert inlet
   `(7,119)` → trap 233, outlet `(199,4)` → trap 13.  1 network, 1 culvert,
   registered on the inlet trap's `culvert_inlets` and the outlet trap's
   `culvert_outlets`.
2. **Cross-network merge** — starts `[(195,7),(179,37)]` (2 disjoint nets); culvert
   `(179,37)`→`(196,6)` merges them into 1 network containing the culvert.
3. **Non-inclusion** — culvert `(1,1)`→`(10,10)` (both off the `(179,37)` net) is
   not added; network unchanged, `culverts` empty.
4. **Terrain outlet → expansion** — culvert inlet `(179,37)` (in net), outlet
   `(8,119)` (bare slope cell of the long chain); the outlet traces a fresh
   downstream chain, growing trap/path counts (3→60 traps in practice).
5. **Inlet/outlet on flow paths** — culvert `(182,34)`→`(190,31)` (both path cells)
   is registered on the owning paths' `culvert_inlets`/`culvert_outlets`, on no
   trap.

## Task 2 (solver) — IN PROGRESS

**Rate function — DONE.** `culvert_rate` (`src/dynamics/culvert_rate.jl`) computes
capacity from a simplified HDS-5 model (inlet control weir/orifice vs. outlet
control full-barrel/free-outfall, take the more restrictive).  A convenience
`DynCulvert(tstruct, inlet, outlet; r, n, …)` constructor derives `Kf`/length and
fills SI defaults.  Settled simplifications are under "Clarifications" in
`agent/prompts/culvert_rate_implementation.org`; unit tests cover both controls,
the constructor, and the free-outfall branch (`test/dynamics_test.jl`), plus a
visual `examples/verification/culvert_rate.jl`.  Not yet called by the solver.

**Routing integration — NOT STARTED.** `networksolver.jl`/`_route_flow` still has
no culvert handling. Plan (see conversation):
- Precompute per-culvert routing data (inlet/outlet owner: path+cell-pos or trap;
  inlet elevation for submersion events).
- Handle culvert inlets/outlets inside `_route_flow`, in the unified per-path event
  traversal (inlet = like high infiltration at its cell; outlet = like a merge at
  its cell). **Mass conservation is paramount**: the actual flow drawn at a culvert
  inlet (possibly capped by available flow) must equal the amount delivered at the
  outlet — track `culvert_actual_delivered[ci]` and use only that on the outlet
  side. See AGENTS.md (Mass conservation).
- Trap-inlet culverts drain the trap; trap-outlet culverts add to it — all computed
  inside `_route_flow` using topological order (inlet processed before outlet).
- Event detection: a trap-inlet culvert's submersion status flips (water level
  crosses the inlet elevation) → termination event (`:culvert_inlet_change`).
  Outlet submersion does NOT trigger an event (only affects rate).
- **Solver ordering must respect culvert direction.** Construction now does this:
  `_topological_order` (in `_build_component`) adds an inlet-owner→outlet-owner
  edge per culvert, so each built network is ordered inlet-before-outlet and an
  uphill/reverse culvert (which would make the network cyclic) is rejected at
  construction. But the solver re-derives its own order via `_network_order`
  (`networksolver.jl`), which still ignores culverts — when routing lands it must
  add the same culvert edges (or trust the construction order) so `_route_flow`
  processes each culvert's inlet before its outlet. Uphill/reverse culverts stay
  deferred (currently they fail loud at construction).

## Design decisions already settled (with the user)

- `DynFlowPath.merges` is `Vector{Tuple{Int,Int}}` (tributary idx, junction cell idx).
- Inclusion rule (final, simple): a culvert is in the network iff an endpoint cell
  is in an in-network trap footprint or on an in-network flow path.
- Expansion traces a fresh `_subnetwork` from BOTH endpoints of each included culvert.
- `culvert_rate` is a stub for now; real hydraulics come in task 2's second phase.
- One routing code path only (no separate approximate path).
