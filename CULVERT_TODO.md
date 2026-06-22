# Culvert support — work in progress (note to self)

## Status

**Task 1 (network construction) — partially done.** `setup_network` now takes an
optional `culverts` kwarg, expands the network around touched culverts, dedups
traps, groups paths+traps into components (with culvert links), and assigns
culverts to their owning path/trap. All existing dynamics tests pass plus new
unit tests for `_components` (cross-component culvert link, lone trap survival).

**Task 1 integration tests on mini.txt — NOT YET WRITTEN.** This is the main
outstanding item before moving on. Need integration tests that:

1. **Inclusion / assignment**: build the `:long` network (start `(7,119)`), add a
   `DynCulvert` whose inlet is a footprint cell of one in-network trap and whose
   outlet is a footprint cell of another. Assert the result is still 1 network,
   `length(net.culverts) == 1`, the inlet trap lists it in `culvert_inlets`, and
   the outlet trap in `culvert_outlets`.
2. **Cross-network merge**: starts `[(195,7), (179,37)]` give 2 disjoint networks.
   A culvert from a trap in net 1 to a trap in net 2 should yield a single merged
   network containing the culvert.
3. **Non-inclusion**: a culvert whose inlet and outlet are both on cells not in the
   network (scan the grid for a cell outside `_occupied_cells`) must NOT be added;
   the network is unchanged and `net.culverts` is empty.
4. **Outlet in terrain → downstream expansion**: a culvert whose outlet lands on a
   bare-terrain cell should trace a new downstream flow path / trap that becomes
   part of the network (verify trap/path count grows).
5. **Inlet/outlet on a flow path** (not a trap): assert the culvert is registered
   on the path's `culvert_inlets` / `culvert_outlets` and the junction handling is
   correct once the solver side lands.

Helper needed in tests: convert `ts.footprints[trap_ix]` linear indices to
`CartesianIndex` via `CartesianIndices(ts.topography)` to pick endpoint cells.

## Task 2 (solver) — NOT STARTED

`networksolver.jl` still has no culvert handling. Plan (see conversation):
- `culvert_rate(culvert, ...)` **mock stub** returning a fixed constant for now.
- Precompute per-culvert routing data (inlet/outlet owner: path+cell-pos or trap;
  inlet elevation for submersion events).
- Handle culvert inlets/outlets inside `_route_flow`, in the unified per-path event
  traversal (inlet = like high infiltration at its cell; outlet = like a merge at
  its cell). **Mass conservation is paramount**: the actual flow drawn at a culvert
  inlet (possibly capped by available flow) must equal the amount delivered at the
  outlet — track `culvert_actual_delivered[ci]` and use only that on the outlet
  side. See CLAUDE.md.
- Trap-inlet culverts drain the trap; trap-outlet culverts add to it — all computed
  inside `_route_flow` using topological order (inlet processed before outlet).
- Event detection: a trap-inlet culvert's submersion status flips (water level
  crosses the inlet elevation) → termination event (`:culvert_inlet_change`).
  Outlet submersion does NOT trigger an event (only affects rate).

## Design decisions already settled (with the user)

- `DynFlowPath.merges` is `Vector{Tuple{Int,Int}}` (tributary idx, junction cell idx).
- Inclusion rule (final, simple): a culvert is in the network iff an endpoint cell
  is in an in-network trap footprint or on an in-network flow path.
- Expansion traces a fresh `_subnetwork` from BOTH endpoints of each included culvert.
- `culvert_rate` is a stub for now; real hydraulics come in task 2's second phase.
- One routing code path only (no separate approximate path).
