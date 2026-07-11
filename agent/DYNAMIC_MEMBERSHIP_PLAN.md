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
- **No source flag.** A `dyn_coord`-anchored trap (the trap *is* the seed, no
  feeding path) is seeded with `in_count = 1` at build; that `+1` is never
  decremented (decrements fire only per feeding-path death), so it floors at 1.
  Culvert/NBS-fed traps need nothing special — their feeding source-path is never
  orphaned (orphaning fires only when a *trap* stops spilling), so it persists and
  keeps `in_count ≥ 1`.
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

`in_count` = live incoming spill paths (transient) **plus a permanent floor** of
one per persistent coupling to a dynamic element that is *not* a path: a culvert
with an endpoint in the trap (intrinsic, `trap.culvert_*`), an NBS accumulation
coupling (`nbs_accum_traps`, from the split's accumulation resolution), and a
dyn_coord anchor (`dyn_coord_traps`, a setup input). Floor terms are never
decremented.
```julia
function init_in_counts!(net::DynNetwork, dyn_coord_traps, nbs_accum_traps)
    for t in net.traps; t.in_count = 0; end
    for p in net.flow_paths
        p.target_trap > 0 && (net.traps[p.target_trap].in_count += 1)   # spill / source paths
    end
    for t in net.traps                                                  # direct culvert coupling
        (!isempty(t.culvert_inlets) || !isempty(t.culvert_outlets)) && (t.in_count += 1)
    end
    for i in dyn_coord_traps; net.traps[i].in_count += 1; end           # dyn_coord seed
    for i in nbs_accum_traps; net.traps[i].in_count += 1; end           # NBS accumulation
    return net
end
```
(A culvert/NBS feed arriving *via a path* is already counted through `target_trap`
and its source-headed path is never orphaned; the floor is for the *direct*
couplings that are not paths.) Reusing this to `@assert` the incremental value in
tests is a test-only nicety, not a production path.

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
   tributaries (re-based positions). `P.target_trap`'s `in_count` is **unchanged**
   (still fed). A path promoted from a culvert/NBS outlet is source-headed, so it
   is never orphaned thereafter.
3. **None** → `P` fully dies: `t = P.target_trap`; `t.in_count -= 1`; if
   `t.in_count == 0` → `t` detaches → recurse via `detach_spill!(net, t)`.

Invariant that keeps "uppermost live" cheap: the net holds **only live paths**
(re-root/remove eagerly), so uppermost-live == uppermost-present — no recursive
aliveness check on the tributary. Use a worklist; mind ordering when a tributary
and its host die in one sweep.

## 4. Build-time init — `setup_network`

After the monolithic net is built, call
`init_in_counts!(net, dyn_coord_traps, nbs_accum_traps)` once (paths + the
permanent floor); the split then copies each `in_count` through `_localize_trap`,
since counts are component-invariant. `dyn_coord_traps` is a setup input;
`nbs_accum_traps` comes from the split's accumulation resolution; culvert coupling
is read off the trap.

## 5. Build prerequisites (in place, commit `d131a30`)

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
