# Dynamic-network membership: a reachability counter on traps.  A trap is dynamic while
# `in_count > 0`.  See agent/DYNAMIC_MEMBERSHIP_PLAN.md.

export init_in_counts!, detach_spill!

# ----------------------------------------------------------------------------
"""
    init_in_counts!(net, dyn_coord_traps, nbs_accum_traps) -> net

Set every trap's `in_count`: live incoming spill paths (`target_trap`) plus a permanent
floor of one per direct coupling to a dynamic element that is not a path — a culvert with
an endpoint in the trap, an NBS accumulation coupling (`nbs_accum_traps`), or a dyn_coord
anchor (`dyn_coord_traps`).  Floor terms are never decremented.  Run once at build; the
split copies counts through `_localize_trap`.
"""
function init_in_counts!(net::DynNetwork,
                         dyn_coord_traps = Int[],
                         nbs_accum_traps = Int[])
    for t in net.traps; t.in_count = 0; end
    for p in net.flow_paths
        p.target_trap > 0 && (net.traps[p.target_trap].in_count += 1)
    end
    for t in net.traps
        (!isempty(t.culvert_inlets) || !isempty(t.culvert_outlets)) && (t.in_count += 1)
    end
    for i in dyn_coord_traps; net.traps[i].in_count += 1; end
    for i in nbs_accum_traps; net.traps[i].in_count += 1; end
    return net
end

# ----------------------------------------------------------------------------
"""
    detach_spill!(net, trap_id) -> Vector{Int}

Trap `trap_id` stopped overflowing: clear its spill edge and re-root-or-kill the orphaned
spill path, cascading downstream.  Returns the local indices of traps that thereby left
the dynamic network (for the caller to hand back to static handling).
"""
function detach_spill!(net::DynNetwork, trap_id::Int)
    detached = Int[]
    sp = net.traps[trap_id].spill_path
    net.traps[trap_id].spill_path = 0
    sp > 0 && _reroot_or_kill_path!(net, sp, detached)
    return detached
end

# The spill path lost its head: promote its uppermost live injection, else kill it.
function _reroot_or_kill_path!(net::DynNetwork, path_id::Int, detached::Vector{Int})
    j = _uppermost_injection(net.flow_paths[path_id])
    j === nothing ? _kill_path!(net, path_id, detached) : _promote!(net, path_id, j)
end

# Uppermost injection position on a path — a merge junction or a culvert/NBS outlet — or
# nothing.  culvert_inlets draw rather than feed, so are excluded.
function _uppermost_injection(fp::DynFlowPath)
    j = nothing
    upd(pos) = (j === nothing || pos < j) && (j = pos)
    for (_, pos) in fp.merges;          upd(pos); end
    for (_, pos) in fp.culvert_outlets; upd(pos); end
    for (_, pos) in fp.nbs_outlets;     upd(pos); end
    return j
end

# No live injection: the path dies.  If it targets a trap, decrement it and cascade when
# that trap loses its last feed; if it is a tributary/domain-exit path, just remove it.
function _kill_path!(net::DynNetwork, path_id::Int, detached::Vector{Int})
    t = net.flow_paths[path_id].target_trap
    if t > 0
        _tombstone_path!(net, path_id)
        net.traps[t].in_count -= 1
        if net.traps[t].in_count == 0
            push!(detached, t)
            sp = net.traps[t].spill_path
            net.traps[t].spill_path = 0
            sp > 0 && _reroot_or_kill_path!(net, sp, detached)
        end
    else
        _detach_from_host_merges!(net, path_id)
        _tombstone_path!(net, path_id)
    end
end

# Drop this (dying tributary) path from whatever path lists it in `merges`.  No trap is
# decremented — a tributary feeds a path, not a trap.
function _detach_from_host_merges!(net::DynNetwork, path_id::Int)
    for host in net.flow_paths
        filter!(m -> m[1] != path_id, host.merges)
    end
end

# Replace a path slot with a dead marker (departure_point == (0,0)) without shifting
# indices; live references are updated by the caller before this is reached.
_tombstone_path!(net::DynNetwork, path_id::Int) =
    net.flow_paths[path_id] = DynFlowPath(CartesianIndex{2}[], CartesianIndex(0, 0), 0,
        Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[])

# @@@ Phase 1c: re-root (promote a live tributary/culvert/NBS injection to take over the
#     orphaned path's downstream role).  Until then a path with an injection cannot be
#     detached.
_promote!(net::DynNetwork, path_id::Int, j::Int) =
    error("_promote! (re-root) not yet implemented — Phase 1c")
