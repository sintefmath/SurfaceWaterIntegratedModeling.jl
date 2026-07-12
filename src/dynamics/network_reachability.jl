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
    for t in net.traps; t.in_count = 0; end
    for p in net.flow_paths
        p.target_trap > 0 && (net.traps[p.target_trap].in_count += 1)
    end
    for t in net.traps
        isempty(t.culvert_inlets) || (t.in_count += 1)
    end
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

# Re-root: the injection at position j takes over the orphaned path's downstream role.
# The dead head stub (cells 1:j-1) is dropped; the target trap stays fed (in_count
# unchanged).  Injections all sit at positions >= j, so none is lost.
# @@@ a culvert *inlet* in the dropped stub (position < j) is dropped with it.
function _promote!(net::DynNetwork, path_id::Int, j::Int)
    P  = net.flow_paths[path_id]
    mi = findfirst(m -> m[2] == j, P.merges)
    mi === nothing ? _promote_source!(net, path_id, j) :
                     _promote_tributary!(net, path_id, j, P.merges[mi][1])
end

# Injection at j is a tributary Q: the survivor is Q followed by P from the junction down.
function _promote_tributary!(net::DynNetwork, path_id::Int, j::Int, q::Int)
    P, Q = net.flow_paths[path_id], net.flow_paths[q]
    lq = length(Q.cells)
    reb(p) = lq + p - j + 1
    ci = vcat(Q.culvert_inlets,  [(c, reb(p)) for (c, p) in P.culvert_inlets  if p >= j])
    co = vcat(Q.culvert_outlets, [(c, reb(p)) for (c, p) in P.culvert_outlets if p >= j])
    no = vcat(Q.nbs_outlets,     [(n, reb(p)) for (n, p) in P.nbs_outlets     if p >= j])
    mg = vcat(Q.merges,          [(m, reb(p)) for (m, p) in P.merges          if p != j])
    net.flow_paths[path_id] = DynFlowPath(vcat(Q.cells, P.cells[j:end]),
        Q.departure_point, P.target_trap, ci, co, no, mg)
    for t in net.traps; t.spill_path == q && (t.spill_path = path_id); end  # redirect Q's head
    _tombstone_path!(net, q)
end

# Injection at j is a culvert/NBS outlet: the survivor is P from j down, now source-headed.
function _promote_source!(net::DynNetwork, path_id::Int, j::Int)
    P = net.flow_paths[path_id]
    reb(p) = p - j + 1
    ci = [(c, reb(p)) for (c, p) in P.culvert_inlets  if p >= j]
    co = [(c, reb(p)) for (c, p) in P.culvert_outlets if p >= j]
    no = [(n, reb(p)) for (n, p) in P.nbs_outlets     if p >= j]
    mg = [(m, reb(p)) for (m, p) in P.merges          if p >= j]
    net.flow_paths[path_id] = DynFlowPath(P.cells[j:end], P.cells[j], P.target_trap, ci, co, no, mg)
end
