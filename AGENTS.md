# AGENTS.md — SurfaceWaterIntegratedModeling.jl

## Project overview

SWIM is a Julia package for static and dynamic surface water and flood modelling
based on terrain topography. It performs spill-point analysis (originally
developed for CO₂ storage) adapted to model water accumulation, watershed
delineation, and urban flooding — without numerical time-stepping.

The central data structure is `TrapStructure` (`src/TrapStructure.jl`), which
holds the raster terrain grid, flow graph, trap regions, spillpoints, volumes,
and the sub/supertrap hierarchy. Nearly every algorithm reads from or writes
into one of these objects.

Key subsystems:

| Area                        | Files                              |
|-----------------------------|------------------------------------|
| Flow field & regions        | `spillfield.jl`, `spillregions.jl` |
| Spill-point geometry        | `spillpoints.jl`                   |
| Trap volumes & hierarchy    | `trapvolumes.jl`, `sshierarchy.jl` |
| Full static analysis        | `spillanalysis.jl`                 |
| Dynamic fill sequencing     | `fill_sequence/` (six files)       |
| Flow intensity over terrain | `watercourses.jl`                  |
| IO & visualisation          | `IOandplot.jl`                     |

## Agent context folder

Prompts, task specs, and background/working documents that support agent work
live under `agent/`:

| Path                | Contents                                             |
|---------------------|------------------------------------------------------|
| `agent/prompts/`    | Source prompts and original task instructions        |
| `agent/reports/`    | Generated analysis/review reports                    |
| `agent/CULVERT_TODO.md` | Active culvert work tracker                       |

Read these for background, but this file (`AGENTS.md`) remains the authoritative
standing guidance.

## Language and conventions

- Julia 1.9+ (CI tests on 1.9; compatible with current stable).
- No linter is enforced, but the codebase uses 4-space indentation throughout.
- Exported public API is documented with docstrings; internal helpers are
  prefixed with `_` and are not exported.
- Type parameters use `{T<:Real}` / `{<:Real}` consistently; avoid
  introducing concrete-type constraints in new function signatures.
- Parallelism is not used yet (commented-out tiling tests exist in
  `test/runtests.jl`); do not add threading or task-based concurrency without
  discussion.

## Running tests

Tests require the `swim_testdata` artifact (downloaded automatically via
`LazyArtifacts`).

```julia
using Pkg
Pkg.test()           # from the package root
```

Or in the REPL:

```
] activate .
] test
```

Tests are deterministic — they compare checksums/counts, not floating-point
tolerances — so a failing test always means something changed in the algorithm.
Do not update expected values to make a test pass; fix the underlying logic.

## Making changes

1. Branch from `main`. Branch names follow a short, descriptive kebab-case
   pattern (e.g. `fix-compute-spillfield-graph-missing-vertices`).
2. All public functions must retain their existing keyword-argument interface
   unless the change is intentionally breaking (increment minor version and
   note in PR description).
3. After any change to `fill_sequence/`, run the `Sequencing` and
   `Trapping structure` test sets at minimum.
4. After any change to `spillfield.jl` or `spillregions.jl`, run
   `SpillfieldTests` and `SpillregionTests`.
5. Submit a PR against `main`. CI runs tests and builds docs automatically.

## Dynamic fill sequencing (`fill_sequence/`)

This is the most algorithmically complex subsystem. Key invariants to preserve:

- `fill_sequence` returns a `Vector{SpillEvent}` where events are strictly
  ordered by timestamp.
- `SpillEvent` entries can be *full* (storing complete state snapshots) or
  *incremental* (storing only `IncrementalUpdate` diffs). The first event in
  each weather period is always a full snapshot.
- `_identify_next_status_change!` mutates `changetimeest`, `filled_traps` (via
  a temporary flip), and `cur_amounts`, but must leave them consistent on
  return.
- `fill_trap_until` integrates water-level ODEs using
  `DifferentialEquations.jl`. The ODE terminates via `VectorContinuousCallback`
  with three conditions: trap full, trap empty, or rate sign-change. Do not
  remove any of these three conditions.

## Dos and don'ts for agents

**Do:**
- Read `TrapStructure` field definitions in `TrapStructure.jl` before touching
  any algorithm — the struct has grown over time and some fields are optional
  or context-dependent.
- Check `spillanalysis.jl` for how `TrapStructure` objects are constructed
  end-to-end before adding new fields to the struct.
- Prefer explicit loops over broadcast/matrix operations in inner loops on large
  grids (see comment at the top of `spillfield.jl` — profiling confirmed this).
- Keep `watercourses.jl` aware of `waterbodies` (see commit `b044c45`); any
  flow-routing change must account for pre-filled water cells.

**Don't:**
- Don't mock or stub the test data artifacts — tests must use real terrain grids
  to catch numerical regressions.
- Don't change expected values in `runtests.jl` without a clear algorithmic
  justification in the PR.
- Don't add `Makie`-dependent code outside `IOandplot.jl`; plotting is
  intentionally isolated so the package can be used headlessly.
- Don't add new top-level exports without updating the docstring and adding at
  least a basic test.

## Documentation

Docs live in `docs/src/` and are built with Documenter.jl. Run locally with:

```julia
julia --project=docs docs/make.jl
```

Examples in `examples/` have their own `Project.toml` and must be run with
that environment activated (see README for the step-by-step guide).

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

The `culvert_rate` function (`src/dynamics/culvert_rate.jl`) computes capacity
from a simplified HDS-5 model: inlet control (weir below submergence, orifice
above) vs. outlet control (full-barrel pressure flow, free-outfall tailwater
when the outlet is dry), returning the more restrictive of the two.  Heads are
water level above the cell invert; the settled simplifications are recorded under
"Clarifications" in `agent/prompts/culvert_rate_implementation.org`.  It accepts
`allow_reverse` (default `false`): when `true` a drowned culvert returns a
*negative* (outlet→inlet) rate instead of `0`, but the routing layer does not yet
consume negative rates.  The outlet-control friction term uses the HDS-5
hydraulic-radius form `Kf = Ku·n²·L / R^(4/3)` (`Ku = 19.63` SI, `R = D/4` for a
full circular pipe); the 1 m/cell grid-resolution assumption is still taken on
trust and flagged in-code with `@@@`.  The function is wired into the solver via `_culvert_flow` → `_path_delivered!` → `_route_flow`.

## Local paths

`examples/Project.toml` has a machine-local `[sources]` path and is protected
with `git update-index --skip-worktree`.  Do not commit it.  `examples/Manifest.toml`
is in `.gitignore` for the same reason.  See memory file `local_path_overrides.md`.
