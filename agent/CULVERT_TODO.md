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

**Routing integration — IN PROGRESS.**

*Increment 1 — trap↔trap culverts — DONE.* `_route_flow` now routes culverts whose
both endpoints own traps: per call it evaluates `culvert_rate` from the two traps'
live water levels (`CulvertPlan` precomputes owners/diameters/topography handle;
`dynNetworkRateFunction!` passes `trap_level = water_level(geom, V)`), then **drains
the inlet trap and delivers the same amount to the outlet trap** (`trap_inflow[in]
-= q; trap_inflow[out] += q`), so `dV = inflow − loss − spill` stays mass-exact.
Downhill-only (`allow_reverse=false`); a would-be reverse yields 0.  `_network_order`
now adds the inlet→outlet culvert edge (the solver-side ordering fix).  Mass-
conservation unit test on a two-trap network; end-to-end solve verified on mini.txt.

*Remaining:*
- **Flow-path endpoints (increment 2).** A culvert inlet/outlet on a flow path:
  inlet = abstract at its cell position, **capped by the flow passing there**;
  outlet = add the delivered amount at its position (like a `merges` junction).
  Path endpoints are treated as not-submerged with head = diameter `D` (so capacity
  is `weir(D)`).  Needs the unified per-path traversal (tributaries + culvert
  inlets/outlets merged by cell position).  Currently such culverts carry 0.
- **Unspill event for culvert-drained full traps (increment 3).** A full trap drained
  by a culvert until its net crosses zero should raise `:unspill`, but that condition
  is still dormant (`condition` outputs `1.0`).  The `:unspill` scaffold should also
  be enabled for traps with culvert inlets/outlets (state-dependent net).  No new
  *submersion* event is needed (settled: events are for topology changes only; a
  culvert turning on/off or switching weir↔orifice doesn't change connectivity).
- Uphill/reverse culverts stay deferred (rejected at construction; reverse → 0).

## Design decisions already settled (with the user)

- `DynFlowPath.merges` is `Vector{Tuple{Int,Int}}` (tributary idx, junction cell idx).
- `DynFlowPath.culvert_inlets`/`culvert_outlets` are `Vector{Tuple{Int,Int}}`
  `(culvert idx, cell position in this path's cells)` — the position lets routing
  charge infiltration up to the abstraction/addition point, like a `merges` junction.
  `DynTrap.culvert_inlets`/`culvert_outlets` stay `Vector{Int}` (a trap is a reservoir
  with no along-path position).  The inlet/outlet *elevation* needed for trap-inlet
  submersion events is NOT stored in `DynTrap`: the cell is already in `DynCulvert`
  and elevations are derived from `tstruct`, so the solver precomputes it per culvert
  (see the routing-integration bullet on `inlet elevation for submersion events`).
- Inclusion rule (final, simple): a culvert is in the network iff an endpoint cell
  is in an in-network trap footprint or on an in-network flow path.
- Expansion traces a fresh `_subnetwork` from BOTH endpoints of each included culvert.
- One routing code path only (no separate approximate path).
