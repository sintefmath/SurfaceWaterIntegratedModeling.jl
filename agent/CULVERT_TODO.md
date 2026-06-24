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

**Routing integration — DONE** (downhill, no-reverse scope).

*Increment 1 — trap↔trap culverts.* `_route_flow` evaluates `culvert_rate` from the
two traps' live surface levels (`CulvertPlan` precomputes owners/diameters/topography
handle; `_routed_inflow` passes `trap_level`), then **drains the inlet trap and
delivers the same amount to the outlet trap**, so `dV = inflow − loss − spill` stays
mass-exact.  `_network_order` adds the inlet→outlet edge (solver-side ordering).

*Increment 2 — flow-path endpoints.* A culvert inlet/outlet on a flow path is handled
inside the path's segmented traversal: tributary junctions and culvert inlet/outlet
positions are merged into one in-order stream.  A path inlet abstracts **capped by
the flow passing its cell**; a path outlet adds the source's drawn amount.  Path
endpoints are not-submerged with head = diameter `D` (capacity `weir(D)`).  All four
endpoint combinations now carry flow with drawn == delivered.

*Increment 3 — general `:unspill` event.* A full trap begins draining when its net
inflow (`inflow − footprint loss`) crosses zero; this is now monitored for **every**
full trap (not culvert-specific), since the net can vary with a feeding/draining
culvert and, later, with NBS inflows.  `_routed_inflow` exposes the raw inflow the
clamped rate hides; the `condition` feeds the net into the `:unspill` root.  Inert
for constant-inflow networks (never crosses).  No *submersion* event (settled:
events are for topology changes only).  Verified: `culvert_rate` heads use the real
surface level (`_surface_level`, not the `Inf` full-trap sentinel — `Inf − Inf` was
giving `NaN`); end-to-end a culvert drain triggers `:unspill` on mini.txt.

*Still deferred:* runtime reverse flow (uphill culverts rejected at construction;
reverse currently yields 0).

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
