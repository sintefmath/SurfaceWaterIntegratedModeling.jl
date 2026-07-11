# Dynamic-network membership: a reachability counter on traps.  A trap is dynamic while
# `in_count > 0`.  See agent/DYNAMIC_MEMBERSHIP_PLAN.md.

export init_in_counts!

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
