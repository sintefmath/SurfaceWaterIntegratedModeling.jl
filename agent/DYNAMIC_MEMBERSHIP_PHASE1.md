# Dynamic-network membership — Phase 1 (intended changes)

Concrete change set for Phase 1 of `agent/DYNAMIC_MEMBERSHIP_PLAN.md`: data +
primitives + isolated tests. Build-time wiring and the live-lifecycle hook are
Phase 2. Nothing here is written yet.

Prerequisite already in place (commit `d131a30`): every trap-to-trap link is a
path — adjacent traps get a zero-length connector with `target_trap` set — so
`init_in_counts!`'s `p.target_trap` count picks them up. See `PLAN.md` §5.

## No source flag — seed anchors counted via connectors

There is no source flag. A trap anchored by a seed (dyn_coord / culvert or NBS
outlet / NBS accumulation) is fed by the persistent source-headed **connector path**
its seed grows (zero-length if the seed sits in the trap), which targets the trap and
is therefore counted via `target_trap`. That path is never orphaned (orphaning fires
only when a *trap* stops spilling, and a seed-headed path has no trap head), so the
count never drops. **Aliveness is uniformly `in_count > 0`.**

The only coupling needing a floor is a culvert **inlet** in a trap — not a seed (no
connector) and it draws rather than feeds (no `target_trap`).

## 1. `src/dynamics/elements.jl` — one field on `DynTrap`

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

The 4-arg inner constructor keeps every existing `DynTrap(ix, sp, ci, co)` call
and the 2-arg convenience working. `DynFlowPath` is untouched (no path counter,
stays immutable — re-root replaces the object).

## 2. New `src/dynamics/network_reachability.jl` (+ `include` after `network_utils.jl`)

**Build-time seed pass** (once, on the monolithic net). Counts are
component-invariant (a trap's feeders are all in its component), so the split just
**copies** each `in_count` through `_localize_trap` — no per-component recompute.

Implemented as `init_in_counts!(net)` — floor simplified to culvert **inlets** only
(commit `106a897`): every seed anchor (dyn_coord / culvert or NBS outlet / NBS
accumulation) is already counted via the persistent connector path it grows, so only a
culvert inlet (not a seed, and it draws) still needs a floor. Counts are
component-invariant, so the split copies them through `_localize_trap`.
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
(A culvert/NBS feed arriving *via a path* is already counted through `target_trap`
and its source-headed path is never orphaned; the floor is for the *direct*
couplings that are not paths.)

**Grow:** no function — the call site does `net.traps[B].in_count += 1` when a new
spill path targeting `B` is added. No cascade (a new feed can't detach anything).

**Shrink / detach:**
```julia
# trap_id stopped overflowing. Returns detached trap ids (to hand off to static handling).
function detach_spill!(net::DynNetwork, trap_id::Int)
    detached = Int[]
    sp = net.traps[trap_id].spill_path
    net.traps[trap_id].spill_path = 0
    sp != 0 && _reroot_or_kill_path!(net, sp, detached)
    return detached
end

# The orphaned path lost its head: promote its uppermost live injection, else kill it
# and cascade to its target trap.
function _reroot_or_kill_path!(net, path_id, detached)
    fp = net.flow_paths[path_id]
    j  = _uppermost_injection(fp)           # min pos over merges + culvert/NBS outlets
    if j !== nothing
        net.flow_paths[path_id] = _promote(net, fp, j)     # target trap unchanged
    else
        t = fp.target_trap
        _remove_path!(net, path_id)                        # tombstone, do not shift indices
        if t > 0
            net.traps[t].in_count -= 1
            if net.traps[t].in_count == 0
                push!(detached, t)
                sp = net.traps[t].spill_path
                net.traps[t].spill_path = 0
                sp != 0 && _reroot_or_kill_path!(net, sp, detached)
            end
        end
    end
end
```

Careful bits (to solve during implementation):

- `_uppermost_injection(fp)` = smallest `pos` over `fp.merges` ∪ `fp.culvert_outlets`
  ∪ `fp.nbs_outlets`. `culvert_inlets` are excluded — they draw, not feed.
  Tributaries in `merges` are live by the "only live paths present" invariant, so
  no recursive aliveness check.
- `_promote(net, fp, j)` rebuilds a `DynFlowPath`: culvert/NBS injection at `j` →
  new cells `fp.cells[j:end]`; tributary `Q` → `Q.cells ++ fp.cells` after the
  junction. Same `target_trap`; **re-base** surviving lower injections by the cut;
  tombstone the folded `Q`. Target trap `in_count` unchanged.
- `_remove_path!` = **tombstone, not `deleteat!`** — deleting would shift indices
  and break `target_trap` / `spill_path` / `merges`. Mark inactive; compact after
  the cascade settles (and hand detached traps off then).

## 3. Deferred to Phase 2

The build-time `init_in_counts!(net)` call already lives in `setup_network` and the
split copies counts through `_localize_trap` (done, `106a897`; validated on real
terrain without a solver). What remains for Phase 2: the live "trap stopped
overflowing" hook calling `detach_spill!` (the transition lives in the solve /
`fill_sequence` layer, which does not yet drive the new `build_network`
representation), plus grow / state-handoff / re-partitioning.

## 4. Tests (scratchpad harness first, like the split tests)

Chain and diamond `DynNetwork`s:
- `init_in_counts!` gives expected per-trap counts;
- `detach_spill!` on a mid trap returns exactly the traps that lose all feeds;
  survivors keep `in_count > 0`; a seed-anchored trap (counted via its connector)
  never detaches;
- re-root: a path with a live tributary keeps its target fed and the survivor
  emanates from the tributary's source;
- (test-only) incremental counts == a fresh `init_in_counts!` after a sequence of
  grows/detaches.

Highest-risk piece is `_promote`'s position re-base + tombstone bookkeeping; build
it behind the diamond test first.
