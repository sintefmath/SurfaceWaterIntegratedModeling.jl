# Dynamic-network membership — Phase 1 (intended changes)

Concrete change set for Phase 1 of `agent/DYNAMIC_MEMBERSHIP_PLAN.md`: data +
primitives + isolated tests. Build-time wiring and the live-lifecycle hook are
Phase 2. Nothing here is written yet.

Prerequisite already in place (commit `d131a30`): every trap-to-trap link is a
path — adjacent traps get a zero-length connector with `target_trap` set — so
`init_in_counts!`'s `p.target_trap` count picks them up. See `PLAN.md` §5.

## No `is_source` flag — dyn_coord seed = starting `in_count` of 1

There is no source flag. A trap anchored directly by a `dyn_coord` (the trap *is*
the seed, no feeding path — dyn_coords are mostly a test/debug input and can sit
inside a trap) is seeded with `in_count = 1` at build. That `+1` is never
decremented — a decrement fires only per real feeding-path death, and there are
only as many of those as there were paths — so `in_count` floors at 1 and the
trap never detaches. **Aliveness is uniformly `in_count > 0`.**

Culvert/NBS-fed traps need nothing special: their feeding source-path is never
orphaned (orphaning fires only when a *trap* stops spilling, and a
culvert/NBS/dyn_coord-headed path has no trap head), so it persists and keeps
`in_count ≥ 1`. (Assumes culvert/NBS outlets discharge onto terrain via a path;
an outlet directly inside a trap would need its feed represented — revisit then.)

## 1. `src/dynamics/elements.jl` — one field on `DynTrap`

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

The 4-arg inner constructor keeps every existing `DynTrap(ix, sp, ci, co)` call
and the 2-arg convenience working. `DynFlowPath` is untouched (no path counter,
stays immutable — re-root replaces the object).

## 2. New `src/dynamics/network_reachability.jl` (+ `include` after `network_utils.jl`)

**Build-time seed pass** (once, on the monolithic net). Counts are
component-invariant (a trap's feeders are all in its component), so the split just
**copies** each `in_count` through `_localize_trap` — no per-component recompute.

`in_count` = live incoming spill paths (transient) **plus a permanent floor** of
one per persistent coupling to a dynamic element that is *not* carried by a path:
a culvert with an endpoint in the trap (intrinsic, from `trap.culvert_*`), an NBS
accumulation coupling (`nbs_accum_traps`, from the split's accumulation
resolution), and a dyn_coord anchor (`dyn_coord_traps`, a setup input). Floor
terms are never decremented, so path-deaths can't drop below them.
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

The build-time `init_in_counts!` call (with the dyn_coord and NBS-accumulation
trap sets) lives in `setup_network`, and the split copies counts through
`_localize_trap`. The "trap stopped overflowing" hook lives in the solve /
`fill_sequence` layer, which does not yet drive the new `build_network`
representation. Phase-1 tests build `DynNetwork`s by hand and call the primitives
directly.

## 4. Tests (scratchpad harness first, like the split tests)

Chain and diamond `DynNetwork`s:
- `init_in_counts!` gives expected per-trap counts;
- `detach_spill!` on a mid trap returns exactly the traps that lose all feeds;
  survivors keep `in_count > 0`; a dyn_coord-anchored trap (seeded `in_count = 1`)
  never detaches;
- re-root: a path with a live tributary keeps its target fed and the survivor
  emanates from the tributary's source;
- (test-only) incremental counts == a fresh `init_in_counts!` after a sequence of
  grows/detaches.

Highest-risk piece is `_promote`'s position re-base + tombstone bookkeeping; build
it behind the diamond test first.
