# Dynamic-network membership — plan

Status: design settled, not yet implemented.
Related: `agent/NBS_ROUTING_REDESIGN.md` (the network this runs on),
`src/dynamics/build_network.jl`, `src/dynamics/network_utils.jl`.

## Context

A `DynNetwork` is grown from dynamic **seeds** (culvert outlets, NBS outlets /
outflow cells, `dyn_coords` — the last mostly for testing). As the solve
proceeds, a trap can **stop overflowing**; its spill path then disappears, which
can leave a downstream chunk no longer connected to any seed. Such a chunk is no
longer dynamic and must revert to the static `fill_sequence` handling.

We need to detect that detachment exactly and incrementally. The chosen scheme:
an in-degree **counter on traps**, plus a **re-root-or-die** operation on the
spill path when a trap stops spilling, cascading downstream.

## Design decisions

- **Counter on `DynTrap` only.** Traps are the persistent nodes; flow paths are
  edges. A path survives iff some injection still feeds it (its head trap, a
  tributary, or a culvert/NBS outlet). Keeping the topology exact already forces
  us to *re-root* a surviving path to its live source, so the re-root/die
  decision subsumes any path counter — a path counter would be redundant.
  `DynFlowPath` stays immutable (re-root replaces the object, never mutates it).
- **Sources are the seed-tied elements** (culvert outlets, NBS outlets/outflow
  cells, `dyn_coords`). They are permanent within a solve window, so they never
  detach. Marked `is_source`; not counted.
- **Alive iff `is_source || in_count > 0`.**
- **Incremental during evolution.** Growth increments one trap; shrink decrements
  + cascades. A full count pass is used only at **bulk build** (the split), never
  while the net merely grows or shrinks.

## 1. Data model — `src/dynamics/elements.jl`

Add two fields to `DynTrap` only; keep the existing 4-arg constructor working via
an inner constructor that defaults them.

```julia
mutable struct DynTrap <: DynObject
    trap_ix::Int
    spill_path::Int
    culvert_inlets::Vector{Int}
    culvert_outlets::Vector{Int}
    in_count::Int      # live incoming flow paths
    is_source::Bool    # seed / culvert-outlet / NBS-fed -> never detaches
    DynTrap(ix, sp, ci, co) = new(ix, sp, ci, co, 0, false)
end
```

`DynFlowPath` is unchanged (no counter, stays immutable).

## 2. Primitives — new `src/dynamics/network_reachability.jl` (included after `network_utils.jl`)

**Bulk init (build/split only):** count incoming paths per trap. Run once when a
component's `DynNetwork` is minted (`_build_subnetwork`), where trap indices are
freshly remapped. **Not** run during growth/shrink.
```julia
function init_in_counts!(net::DynNetwork)
    for t in net.traps; t.in_count = 0; end
    for p in net.flow_paths
        p.target_trap > 0 && (net.traps[p.target_trap].in_count += 1)
    end
    return net
end
```
(A path's `target_trap` is its single downstream trap — the only edge feeding a
trap. Culvert/NBS feeds are captured by `is_source` on the fed trap, set at build
time.) Reusing this to `@assert` the incremental value in tests is a test-only
nicety, not a production path.

**Grow (trap starts spilling → new path targets `B`):** local increment, no
cascade — adding a feed can't detach anything.
```julia
B.in_count += 1
```

**Shrink — detach on a trap ceasing to spill:** clear the spill edge, re-root or
kill the path, and cascade. Returns the detached traps for the caller to hand off.
```julia
# trap_id stopped overflowing. Returns detached trap ids (to migrate to static handling).
function detach_spill!(net::DynNetwork, trap_id::Int) -> Vector{Int}
    sp = net.traps[trap_id].spill_path
    net.traps[trap_id].spill_path = 0
    sp == 0 && return Int[]
    _reroot_or_kill_path!(net, sp)   # see §3; decrements/cascades trap in_counts
end
```

## 3. Re-root-or-kill — the centerpiece

When a spill path `P` loses its head (its trap stopped spilling):

1. Scan `P`'s injections top-down (`merges` junctions + culvert/NBS outlet
   positions) for the **uppermost still-present** one.
2. **Found** → promote it: new surviving path = that injection's source path +
   `P` from the junction down to `P.target_trap`; drop the dead head stub; swap
   the new `DynFlowPath` into `net.flow_paths`; keep lower injections as
   tributaries (re-based positions). `is_source` set if promoted from a
   culvert/NBS outlet. `P.target_trap`'s `in_count` is **unchanged** (still fed).
3. **None** → `P` fully dies: `t = P.target_trap`; `t.in_count -= 1`; if
   `!t.is_source && t.in_count == 0` → `t` detaches → recurse via
   `detach_spill!(net, t)`.

Invariant that keeps "uppermost live" cheap: the net holds **only live paths**
(re-root/remove eagerly), so uppermost-live == uppermost-present — no recursive
aliveness check on the tributary. Use a worklist; mind ordering when a tributary
and its host die in one sweep.

## 4. Build-time init — `setup_network`

After the net is built: set `is_source = true` on traps fed by a seed / culvert
outlet / NBS emit, then call `init_in_counts!` so every component starts
consistent. The seed set is known here (it is the input to the growth).

## 5. First-cell-truncation is already safe

A culvert/NBS/`dyn_coord` seed is the first cell of its own path. A terrain path
that passes through such a cell is terrain-downstream of *its* seed, so the
seed-cell has a strictly higher topological rank and is grown first (see
`_seeds_downstream_first`); the passing path then truncates at it (`ix > 1`) and
merges. So promoting seeds to sources introduces no new truncation hazard.

## 6. Phase 2 — live-lifecycle wiring (deferred)

The "trap stopped overflowing" event lives in the dynamic solve / `fill_sequence`
layer, which for the new `build_network` representation is not wired into the
module yet. Phase 1 lands data + primitives + isolated tests only. Phase 2, once
the new net is driven by `fill_sequence`: at the full→not-full transition call
`detach_spill!` and route the returned traps back to static handling, reusing the
existing absorption path — carefully w.r.t. state transfer (cf. the prior
stale-state absorption bugs).

Build-time split (`network_utils.jl`) is unaffected.

## 7. Verification

Isolated harness (as for the split tests): chain and diamond `DynNetwork`s.
- `init_in_counts!` gives expected per-trap counts.
- `detach_spill!` on a mid trap detaches exactly the traps that lose all feeds;
  survivors keep `in_count > 0`; `is_source` roots never detach.
- Re-root case: a path with a live tributary keeps its target trap fed (no
  detach), and the surviving path emanates from the tributary's source.
- (Test-only) incremental counts equal a fresh `init_in_counts!` after a sequence
  of grows/detaches.

## Open questions

- Confirm culverts/NBS are structurally permanent within a window (so
  `is_source` never flips). If a culvert can be removed, it needs the same
  cascade.
- Phase-2 hook location + how a detached subtree's water state migrates to
  static handling.
- Re-root of a promoted tributary whose own source later dies — handled by the
  cascade, but nail the worklist ordering.
