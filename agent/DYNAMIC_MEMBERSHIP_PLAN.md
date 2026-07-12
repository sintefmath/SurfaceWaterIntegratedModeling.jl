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
- **No source flag.** A seed-anchored trap (dyn_coord / culvert or NBS outlet / NBS
  accumulation) is counted by the persistent source-headed connector path its seed
  grows — no floor needed. The only floor is a culvert *inlet* in a trap (a
  persistent coupling that is neither a seed/connector nor a feed). Floor and
  source-headed connectors are never decremented (orphaning fires only when a
  *trap* stops spilling, and neither has a trap head).
- **Alive iff `in_count > 0`.**
- **Incremental during evolution.** Growth increments one trap; shrink decrements
  + cascades. A full count pass is used only at **bulk build** (the split), never
  while the net merely grows or shrinks.

## 1. Data model — `src/dynamics/elements.jl`

Add one field to `DynTrap`; keep the existing 4-arg constructor working via an
inner constructor that defaults it.

```julia
mutable struct DynTrap <: DynObject
    trap_ix::Int
    spill_path::Int
    culvert_inlets::Vector{Int}
    culvert_outlets::Vector{Int}
    in_count::Int      # live incoming flow paths (+1 seed if dyn_coord-anchored)
    DynTrap(ix, sp, ci, co) = new(ix, sp, ci, co, 0)
end
```

`DynFlowPath` is unchanged (no counter, stays immutable).

## 2. Primitives — new `src/dynamics/network_reachability.jl` (included after `network_utils.jl`)

**Build-time seed pass** (once, on the monolithic net). Counts are
component-invariant (a trap's feeders are all in its component), so the split just
**copies** each `in_count` through `_localize_trap` — no per-component recompute.

`in_count` = live incoming paths (`target_trap`) **plus a floor of 1 if a culvert
inlet sits in the trap**. Every seed anchor (dyn_coord / culvert or NBS outlet /
NBS accumulation) already grows a persistent source-headed connector path that
targets its trap, so it is counted via `target_trap` — a floor for those would
double-count. A culvert *inlet* is the only coupling with neither a connector (not
a seed) nor a feed (it draws).
```julia
function init_in_counts!(net::DynNetwork)
    for t in net.traps; t.in_count = 0; end
    for p in net.flow_paths
        p.target_trap > 0 && (net.traps[p.target_trap].in_count += 1)
    end
    for t in net.traps
        isempty(t.culvert_inlets) || (t.in_count += 1)
    end
    return net
end
```
Reusing this to `@assert` the incremental value in tests is a test-only nicety,
not a production path.

**Grow (trap starts spilling → new path targets `B`):** local increment, no
cascade — adding a feed can't detach anything. (Topology growth — trace + attach
the new path/trap — is §4; the increment is the counter side of it.)
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
   tributaries (re-based positions). `P.target_trap`'s `in_count` is **unchanged**
   (still fed). A path promoted from a culvert/NBS outlet is source-headed, so it
   is never orphaned thereafter.
3. **None** → `P` fully dies: `t = P.target_trap`; `t.in_count -= 1`; if
   `t.in_count == 0` → `t` detaches → recurse via `detach_spill!(net, t)`.

Invariant that keeps "uppermost live" cheap: the net holds **only live paths**
(re-root/remove eagerly), so uppermost-live == uppermost-present — no recursive
aliveness check on the tributary. Use a worklist; mind ordering when a tributary
and its host die in one sweep.

## 4. Grow — trace and attach

When a trap `A` already in the dynamic net becomes full and starts spilling, its
spill path and whatever it reaches must be added — mirror
`_grow_network_from_seed!` from `A`'s spillpoint:

1. Trace `A`'s spill (`_trace_to_next_trap`): the connecting path (possibly
   zero-length), any culverts/NBS on it, and the target it reaches.
2. Attach the target:
   - existing dynamic trap `B` → link + `B.in_count += 1`;
   - new trap → add the `DynTrap`, wire its culvert/NBS couplings + floor,
     `in_count = 1` via the new path;
   - path merges into an existing path (pathmap) → register the merge.
3. Set `A.spill_path` to the new connector.

No cascade — a new feed can only keep/revive downstream, never detach it.

**State absorption.** A newly attached trap enters the dynamic net carrying its
current (static) water level; seed its volume from that, not from stale state
(cf. the prior stale-state absorption bugs).

## 5. Component re-partitioning

Grow and detach can change the connected-component partition — each component is
one independent ODE solve, so this matters:

- **Fusion (grow) — required.** If `A`'s new spill path lands in a *different*
  component (`B` has its own seed/net), the two are now coupled; solving them
  apart drops the coupling flow and breaks mass conservation. Detected for free
  at grow time (the target lookup already identifies `B`'s component); merge the
  two (edge insertion = union), or rebuild the pair.
- **Fission (detach) — deferred.** When a chain detaches, the surviving net may
  fall into two pieces that no longer exchange flow. Solving them as one ODE is
  still *correct* (the dead link carries zero flow), just larger. Fission is the
  hard direction incrementally (edge deletion / dynamic connectivity), so keep
  the over-large component and re-split lazily on the next rebuild if perf needs.

## 6. Build-time init — `setup_network`  *(done, `106a897`)*

`setup_network` calls `init_in_counts!(net)` after growth (paths + the culvert-inlet
floor; seed anchors are counted via their connector paths); the split then copies each
`in_count` through `_localize_trap`, since counts are component-invariant. Self-contained
(no `fill_sequence`); verified on the real `mini.txt` terrain.

## 7. Build prerequisites (in place, commit `d131a30`)

The counter relies on the build already guaranteeing:
- **Every trap-to-trap link is a path.** Adjacent traps (spillpoint on the
  downstream footprint, no cell between) are linked by a zero-length connector
  (empty `cells`, `target_trap` set), so `init_in_counts!` counts them like any
  other edge.
- **A path's source is `departure_point`**, valid even when `cells` is empty or
  truncated, so source-ID never depends on `cells[1]`. The old "first cell never
  truncated" assert is gone (a first-cell intersection truncates to a zero-length
  connector); `_seeds_downstream_first` still grows a seed cell before any path
  passing through it, so a seed's own path is never mis-rooted.

## 8. Phase 2 — live-lifecycle wiring (deferred)

The "trap stopped overflowing" event lives in the dynamic solve / `fill_sequence`
layer, which for the new `build_network` representation is not wired into the
module yet. Phase 1 lands data + primitives + isolated tests only. Phase 2, once
the new net is driven by `fill_sequence`: at the full→not-full transition call
`detach_spill!` and route the returned traps back to static handling, reusing the
existing absorption path — carefully w.r.t. state transfer (cf. the prior
stale-state absorption bugs).

Build-time split (`network_utils.jl`) is unaffected.

## 9. Verification

Isolated harness (as for the split tests): chain and diamond `DynNetwork`s.
- `init_in_counts!` gives expected per-trap counts.
- `detach_spill!` on a mid trap detaches exactly the traps that lose all feeds;
  survivors keep `in_count > 0`; a dyn_coord-anchored trap (seeded `1`) never detaches.
- Re-root case: a path with a live tributary keeps its target trap fed (no
  detach), and the surviving path emanates from the tributary's source.
- (Test-only) incremental counts equal a fresh `init_in_counts!` after a sequence
  of grows/detaches.

## Open questions

- Confirm culverts/NBS are structurally permanent within a window (so their
  feeding source-path is never orphaned). If a culvert can be removed, it needs
  the same cascade.
- Phase-2 hook location + how a detached subtree's water state migrates to
  static handling.
- Re-root of a promoted tributary whose own source later dies — handled by the
  cascade, but nail the worklist ordering.
