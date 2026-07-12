# Dynamic-network membership — plan

Status: Phase 1 (data + primitives), build-time init (§6), and the full set-level
mechanism — `grow_spill!`/`detach_spill!` (§3–4), fusion (§5), and the `apply_fill!` /
`apply_unfill!` entry points — implemented + tested standalone on real terrain. What
remains is the **gate**: driving these from `fill_sequence`/`network_context` (state
handoff at the transitions) plus the distributor re-impl / retiring `DynNBS`.  Commits
`2cea3e1`, `75fcbc3`, `282234f`, `106a897`, `96530e4`, `a16a842`, `5ea40a6`, `52d1f71`,
`03a5828`, `3fd1319`.
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
    in_count::Int      # live incoming flow paths (+1 floor if a culvert inlet is in the trap)
    DynTrap(ix, sp, ci, co) = new(ix, sp, ci, co, 0)
end
```

`DynFlowPath` is unchanged (no counter, stays immutable).

## 2. Primitives — `src/dynamics/network_updating.jl` (included after `network_utils.jl`)

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
kill the path, cascade, then **compact in place** (`_compact!`) — one event fires at a
time, so there is no batching to gain by separating the passes. Returns the detached
traps' **`trap_ix`** (stable spillanalysis indices — local indices are renumbered by the
compaction) for the caller to hand off.
```julia
# trap_id stopped overflowing. Returns detached trap_ix (to migrate to static handling).
function detach_spill!(net::DynNetwork, trap_id::Int) -> Vector{Int}
    detached = Int[]                       # local indices during the cascade
    sp = net.traps[trap_id].spill_path
    net.traps[trap_id].spill_path = 0
    sp > 0 && _reroot_or_kill_path!(net, sp, detached)   # §3; decrements/cascades in_counts
    detached_ix = [net.traps[t].trap_ix for t in detached]  # before compaction drops them
    _compact!(net)                         # drop tombstones + in_count==0 traps, reindex
    return detached_ix
end
```
`_compact!` rebuilds each `DynFlowPath` slot with remapped `target_trap`/`merges` and
patches `spill_path` in place (so `DynFlowPath` stays **immutable** — no boxing in the ODE's
hot `flow_paths` vector; `DynNetwork`/`DynTrap` are already mutable). Culverts/NBS are
untouched (a detach removes none).

## 3. Re-root-or-kill — the centerpiece

When a spill path `P` loses its head (its trap stopped spilling):

1. Scan `P`'s `merges` for the **uppermost still-present** tributary junction.
   (Merges only: a culvert/NBS outlet is always seeded as its own connector, which
   enters an intersected path as a co-located merge — or heads its own source path —
   so an outlet never needs promoting on its own.)
2. **Found** → promote it: new surviving path = that tributary `Q` + `P` from the
   junction down to `P.target_trap`; drop the dead head stub; write the survivor into
   **Q's slot** (its head trap already points there — no redirect) and tombstone the
   orphaned head path; keep lower injections as tributaries (re-based positions).
   `P.target_trap`'s `in_count` is **unchanged** (still fed).
3. **None** → `P` fully dies: `t = P.target_trap`; `t.in_count -= 1`; if
   `t.in_count == 0` → `t` detaches → recurse via `detach_spill!(net, t)`.

Invariant that keeps "uppermost live" cheap: the net holds **only live paths**
(re-root/remove eagerly), so uppermost-live == uppermost-present — no recursive
aliveness check on the tributary. Use a worklist; mind ordering when a tributary
and its host die in one sweep.

## 4. Grow — trace and attach  *(done: `grow_spill!` in `network_updating.jl`)*

When a trap `A` already in the dynamic net becomes full and starts spilling, its
spill path and whatever it reaches must be added — `grow_spill!(net, tstruct,
full_traps, trap_id)` mirrors `detach_spill!`, reusing the `_grow_network_from_seed!`
tracer from `A`'s spillpoint:

1. Trace `A`'s spill: the connecting path(s) (possibly zero-length), any culverts/NBS
   on them, and the target(s) reached.  The tracer runs with `departing_trap_ix = A`
   (so `A.spill_path` is set to the first connector) and `stop_at_present = true` (stop
   at the first already-present trap — its downstream is already represented; only the
   connector is attached).  `A` spilling out of the domain → `A.spill_path = -1`, nothing
   attached.
2. `in_count` is then stamped **incrementally over the added elements only** (never a
   full `init_in_counts!` pass): each new path `+1`s its `target_trap`; each new trap
   additionally gets the culvert-inlet floor.  An existing trap the chain terminates into
   is `+1`ed via the connector's `target_trap`; a connector that merges into an existing
   path (target 0) adds nothing.  This is exactly `init_in_counts!` restricted to the new
   elements.
3. Returns the newly-added traps' `trap_ix` (for the caller to migrate from static and
   seed their state).

No cascade — a new feed can only keep/revive downstream, never detach it.  Verified on
real `mini.txt` terrain (inject `build_network`): across ~90 grow cases the incremental
counts equal a fresh `init_in_counts!`, refs stay in range, and no `trap_ix` duplicates.
Cross-component landing (fusion) is **not** handled here — see §5.

**State absorption.** A newly attached trap enters the dynamic net carrying its
current (static) water level; seed its volume from that, not from stale state
(cf. the prior stale-state absorption bugs).

## 5. Component re-partitioning

Grow and detach can change the connected-component partition — each component is
one independent ODE solve, so this matters:

- **Fusion (grow) — done (`apply_fill!` + `_fuse_components!`).** If a grow lands in a
  *different* component the two are now coupled; solving them apart drops the coupling flow
  and breaks mass conservation. Detected against the pre-grow index (`_index_components`):
  a newly-added trap that already lives elsewhere, or a new path cell that falls on another
  component (a shared corridor **or** an NBS footprint — NBS coupling is
  decomposition-relevant even though the rate effect is geometric). The coupled components
  are then **regrown from their seeds** (`_fuse_components!` → `_regrow`): the build tracer
  redoes overlaps, duplicate-trap collapsing, and culvert/NBS assignment for free, so no
  hand-merging is needed. (NBS-crossing paths within the *same* component need nothing —
  the rate function derives capture/emit from footprint geometry at eval time.)
- **Fission (detach) — deferred.** When a chain detaches, the surviving net may
  fall into two pieces that no longer exchange flow. Solving them as one ODE is
  still *correct* (the dead link carries zero flow), just larger. Fission is the
  hard direction incrementally (edge deletion / dynamic connectivity), so keep
  the over-large component and re-split lazily on the next rebuild if perf needs.
  `detach_spill!` still `_compact!`s the surviving component (drops the tombstones /
  dead traps) — that is index cleanup, independent of re-partitioning.

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
  connector). Because of this, **seed order no longer matters** — the downstream-first
  `topological_sort` (an O(terrain) cost paid on every build *and* regrow) was removed;
  builds under any seed order give identical coverage/partition (verified).

## 8. Phase 2 — live-lifecycle wiring

Done and tested standalone (`network_updating.jl`): the per-network mechanisms
`detach_spill!` (+ `_compact!`, §3) and `grow_spill!` (§4), and the set-level entry points
`apply_unfill!` (locate → detach) and `apply_fill!` (locate → grow → detect coupling →
fuse, §5).  Each returns the migrated `trap_ix` (newly-static on unfill, newly-dynamic on
fill) for the caller's state handoff.

What remains is the **gate** — driving them from `fill_sequence`/`network_context`, which
does not yet use the new `build_network` representation:
- at full→not-full call `apply_unfill!` and route the returned `trap_ix` back to static
  handling — carefully w.r.t. state transfer (cf. the prior stale-state absorption bugs);
- at not-full→full call `apply_fill!` and seed the newly-added traps' state from their
  current static level;
- plus the distributor / rate-layer re-impl and retiring `DynNBS`.

Build-time split (`network_utils.jl`) is unaffected.

## 9. Verification *(done for Phase 1 + §6)*

Isolated harness (as for the split tests): chain and diamond `DynNetwork`s.
- `init_in_counts!` gives expected per-trap counts.
- `detach_spill!` on a mid trap detaches exactly the traps that lose all feeds;
  survivors keep `in_count > 0`; a seed-anchored trap (counted via its connector)
  never detaches.
- Re-root case: a path with a live tributary keeps its target trap fed (no
  detach), and the surviving path emanates from the tributary's source.
- (Test-only) incremental counts equal a fresh `init_in_counts!`.

Also integration tests on real `mini.txt` terrain (inject build_network):
- `setup_network` runs `init_in_counts!`, the split copies counts (stored == fresh
  recompute), counts equal incoming-path counts, and `detach_spill!` cascades on a 20+-trap
  chain;
- `grow_spill!` over ~90 fills: incremental counts == fresh recompute, refs in range, no
  duplicate `trap_ix`;
- `_fuse_components!` regrow round-trips single/pair components to the exact trap coverage;
- `apply_fill!` sweep (fills incl. 4 real fusions) and `apply_unfill!` keep counts
  consistent, no duplicate `trap_ix`, no trap loss, migrated traps genuinely new.

Still to add: a composition harness scripting random fill/unfill sequences and asserting the
incremental component set matches a fresh `setup_network` at each step.

## Open questions

- Confirm culverts/NBS are structurally permanent within a window (so their
  feeding source-path is never orphaned). If a culvert can be removed, it needs
  the same cascade.
- Phase-2 hook location + how a detached subtree's water state migrates to
  static handling.
- Re-root of a promoted tributary whose own source later dies — handled by the
  cascade, but nail the worklist ordering.
