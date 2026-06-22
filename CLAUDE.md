# CLAUDE.md — project guidance for AI agents

## Mass conservation

**Mass conservation is paramount.** Water that leaves one element of the network
must arrive at exactly one other element; no flow may be created or destroyed.

This applies everywhere: trap-to-trap routing, path infiltration, tributary
merges, culvert flows, and ODE rate functions.  Concretely:

- In `_route_flow`, the *actual* flow drawn at a culvert/tributary inlet must
  equal the flow delivered at its outlet — never the *requested* rate when
  available flow is insufficient.
- `culvert_actual_delivered[ci]` tracks the real delivered amount and is the
  only quantity used on the outlet side.
- The spill emitted by a full trap must equal `inflow - infiltration - all culvert drains`
  so that `dV ≈ 0` at steady state.

If a design choice would break mass conservation, find a different design.

## Culverts

Culverts are not constrained to follow terrain flow direction.  However the
initial implementation assumes downhill culverts (inlet processed before outlet
in topological order).  Uphill / reverse culverts are deferred; add a TODO if
a case arises.

The `culvert_rate` function is currently a stub returning a fixed constant.
It will be replaced with a proper hydraulic formula (orifice/weir, inlet/outlet
control) once network topology is validated.  Do not rely on its return value
for anything other than topology testing.

## Local paths

`examples/Project.toml` has a machine-local `[sources]` path and is protected
with `git update-index --skip-worktree`.  Do not commit it.  `examples/Manifest.toml`
is in `.gitignore` for the same reason.  See memory file `local_path_overrides.md`.
