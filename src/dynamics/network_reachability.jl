# Dynamic-network membership: a reachability counter on traps.  A trap is dynamic while
# `in_count > 0`.  See agent/DYNAMIC_MEMBERSHIP_PLAN.md.

export init_in_counts!, detach_spill!

# ----------------------------------------------------------------------------
"""
    init_in_counts!(net) -> net

Set every trap's `in_count` = its live incoming paths (`target_trap`) plus a floor of 1
if a culvert *inlet* sits in the trap.  A culvert inlet is a persistent coupling that no
path carries — it draws, and (unlike an outlet) is not a seed.  Every seed anchor
(dyn_coord / culvert or NBS outlet / NBS accumulation) is already counted: each seeds a
persistent source-headed connector path that targets its trap.  Run once at build; the
split copies counts through `_localize_trap`.
"""
function init_in_counts!(net::DynNetwork)
    # initialize all trap counts to zero
    for t in net.traps; t.in_count = 0; end

    # each flow path that targets a trap counts as one feed to that trap
    for p in net.flow_paths
        p.target_trap > 0 && (net.traps[p.target_trap].in_count += 1)
    end

    # each culvert inlet counts as an additional feed (it's not actually feeding
    # the trap, but the additional count keeps the trap alive at all times,
    # ensuring it stays dynamic.
    for t in net.traps
        isempty(t.culvert_inlets) || (t.in_count += 1)
    end
    return net
end

# ----------------------------------------------------------------------------
"""
    detach_spill!(net, trap_id) -> Vector{Int}

Trap `trap_id` stopped overflowing: clear its spill edge and re-root-or-kill the orphaned
spill path, cascading downstream.  Then compact the net in place (drop the tombstoned paths
and detached traps, reindex the survivors).  Returns the `trap_ix` (stable spillanalysis
indices) of the traps that left the dynamic network, for the caller to hand back to static
handling — stable across the compaction, unlike local indices.
"""
function detach_spill!(net::DynNetwork, trap_id::Int)
    detached = Int[]                       # local trap indices, collected during the cascade
    sp = net.traps[trap_id].spill_path
    net.traps[trap_id].spill_path = 0
    sp > 0 && _reroot_or_kill_path!(net, sp, detached)
    detached_ix = [net.traps[t].trap_ix for t in detached]  # grab before compaction drops them
    _compact!(net)
    return detached_ix
end

# Remove the dead nodes a detach cascade leaves — tombstoned paths (departure_point == (0,0))
# and traps with in_count == 0 — and reindex the survivors' cross-references, in place.
# Culverts/NBS are untouched (a detach removes none).  DynNetwork/DynTrap are mutable, so the
# vectors are reassigned and spill_path is patched directly; each DynFlowPath slot is rebuilt.
function _compact!(net::DynNetwork)
    live_p = [p for p in eachindex(net.flow_paths)
              if net.flow_paths[p].departure_point != CartesianIndex(0, 0)]
    live_t = [t for t in eachindex(net.traps) if net.traps[t].in_count > 0]
    pmap = zeros(Int, length(net.flow_paths)); for (n, o) in enumerate(live_p); pmap[o] = n; end
    tmap = zeros(Int, length(net.traps));      for (n, o) in enumerate(live_t); tmap[o] = n; end
    rmp(ix, m) = ix <= 0 ? ix : m[ix]          # 0 (none) / -1 (out-of-domain) pass through

    net.flow_paths = [let p = net.flow_paths[o]
        DynFlowPath(p.cells, p.departure_point, rmp(p.target_trap, tmap),
                    p.culvert_inlets, p.culvert_outlets, p.nbs_outlets,
                    [(pmap[m], j) for (m, j) in p.merges if pmap[m] > 0])
    end for o in live_p]
    net.traps = [let t = net.traps[o]; t.spill_path = rmp(t.spill_path, pmap); t end
                 for o in live_t]
    return net
end

# The spill path lost its head: promote its uppermost live tributary, else kill it.
function _reroot_or_kill_path!(net::DynNetwork, path_id::Int, detached::Vector{Int})
    m = _uppermost_merge(net.flow_paths[path_id])
    m === nothing ? _kill_path!(net, path_id, detached) : _promote_tributary!(net, path_id, m[2], m[1])
end

# Uppermost tributary merging into this path, as (tributary_path, junction_pos), or nothing.
# Only merges feed a re-root: a culvert/NBS outlet is always seeded as its own connector
# path, which enters an intersected path as a co-located merge (or heads its own source
# path), so an outlet never needs promoting on its own — the merge already covers it.
# culvert_inlets draw rather than feed, so they are irrelevant here too.
function _uppermost_merge(fp::DynFlowPath)
    best = nothing
    for (q, pos) in fp.merges
        (best === nothing || pos < best[2]) && (best = (q, pos))
    end
    return best
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

# Re-root: the tributary Q at junction j takes over the orphaned path's downstream role.
# The survivor is Q followed by P from the junction down, kept in Q's slot (it is Q's
# continuation): Q's head trap already spills into q, and path_id's head trap was cleared
# before re-root, so no spill_path redirect is needed — tombstone the orphaned head path
# (path_id) instead.  The dead head stub (P.cells 1:j-1) is dropped; the target trap stays
# fed (in_count unchanged); injections all sit at positions >= j, so none is lost.
# @@@ a culvert *inlet* in the dropped stub (position < j) is dropped with it.
function _promote_tributary!(net::DynNetwork, path_id::Int, j::Int, q::Int)
    P, Q = net.flow_paths[path_id], net.flow_paths[q]
    lq = length(Q.cells)
    reb(p) = lq + p - j + 1
    ci = vcat(Q.culvert_inlets,  [(c, reb(p)) for (c, p) in P.culvert_inlets  if p >= j])
    co = vcat(Q.culvert_outlets, [(c, reb(p)) for (c, p) in P.culvert_outlets if p >= j])
    no = vcat(Q.nbs_outlets,     [(n, reb(p)) for (n, p) in P.nbs_outlets     if p >= j])
    mg = vcat(Q.merges,          [(m, reb(p)) for (m, p) in P.merges          if p != j])
    net.flow_paths[q] = DynFlowPath(vcat(Q.cells, P.cells[j:end]),
        Q.departure_point, P.target_trap, ci, co, no, mg)
    _tombstone_path!(net, path_id)
end
