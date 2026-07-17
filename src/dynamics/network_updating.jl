# Dynamic-network updating: incremental grow/shrink of a DynNetwork as traps fill and
# empty during a solve, tracked by a reachability counter on traps (a trap is dynamic while
# `in_count > 0`).  `grow_spill!` and `detach_spill!` are the two symmetric per-network
# mutations; `apply_fill!` / `apply_unfill!` drive them across the live component set
# (fusing components a grow couples together).

export init_in_counts!, detach_spill!, grow_spill!, apply_fill!, apply_unfill!, apply_empty!

# ----------------------------------------------------------------------------
"""
    init_in_counts!(net) -> net

Recompute every trap's reachability `in_count` from scratch: the live flow paths targeting
it, plus a floor of 1 when a culvert *inlet* sits in the trap (a persistent coupling that
draws but is carried by no path; seed anchors are already counted via the connector paths
they grow).  Full pass — used at build only, never during incremental grow/shrink.

# Arguments
- `net::DynNetwork`: network whose trap counts are (re)set in place.

# Returns
`net`, mutated.
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

Handle trap `trap_id` ceasing to overflow: clear its spill edge, re-root or kill the
orphaned spill path and cascade the reachability decrement downstream, then compact the
network in place (drop tombstoned paths and `in_count==0` traps, reindex survivors).

# Arguments
- `net::DynNetwork`: network, mutated in place.
- `trap_id::Int`: local index (in `net.traps`) of the trap that stopped overflowing.

# Returns
The `trap_ix` (stable spillanalysis indices) of the traps that thereby left the dynamic
network, for the caller to hand back to static handling.  Stable across the compaction,
unlike local indices.
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

"""
    _compact!(net) -> net

Drop the dead elements a detach left behind — tombstoned paths (`departure_point == (0,0)`) and
traps with `in_count == 0` — and reindex the survivors' cross-references.

Mutates and returns `net`: `flow_paths` and `traps` are reassigned to the compacted vectors and
`spill_path` / `target_trap` / `merges` are relabelled.  Culverts and NBS are untouched, since a
detach removes none.  Local indices therefore change; callers hand back stable `trap_ix` values
instead.
"""
function _compact!(net::DynNetwork)
    live_p = [p for p in eachindex(net.flow_paths)
              if net.flow_paths[p].departure_point != CartesianIndex(0, 0)]
    live_t = [t for t in eachindex(net.traps) if net.traps[t].in_count > 0]
    pmap = zeros(Int, length(net.flow_paths)); for (n, o) in enumerate(live_p); pmap[o] = n; end
    tmap = zeros(Int, length(net.traps));      for (n, o) in enumerate(live_t); tmap[o] = n; end
    rmp(ix, m) = ix <= 0 ? ix : m[ix]          # 0 (none) / -1 (out-of-domain) pass through

    net.flow_paths = [let p = net.flow_paths[o]
        DynFlowPath(p.cells, p.departure_point, rmp(p.target_trap, tmap),
                    p.culvert_inlets, p.culvert_outlets, p.nbs_inlets, p.nbs_outlets,
                    [(pmap[m], j) for (m, j) in p.merges if pmap[m] > 0])
    end for o in live_p]
    net.traps = [let t = net.traps[o]; t.spill_path = rmp(t.spill_path, pmap); t end
                 for o in live_t]
    return net
end

"""
    _reroot_or_kill_path!(net, path_id, detached) -> nothing

The orphaned spill path `path_id` lost its head trap: promote its uppermost live tributary to
take over the downstream role, or kill the path outright if it has none.

Mutates `net` and appends to `detached`, the accumulator collecting the local indices of traps
that lose all their feeds as a kill cascades downstream.
"""
function _reroot_or_kill_path!(net::DynNetwork, path_id::Int, detached::Vector{Int})
    m = _uppermost_merge(net.flow_paths[path_id])
    m === nothing ? _kill_path!(net, path_id, detached) : _promote_tributary!(net, path_id, m[2], m[1])
end

"""
    _uppermost_merge(fp) -> (tributary_path, junction_pos) or nothing

The tributary merging into `fp` highest up its length, or `nothing` if it has none.  Nothing is
mutated.

Only merges can feed a re-root.  A culvert or NBS outlet is always seeded as its own connector
path, which enters an intersected path as a co-located merge (or heads its own source path), so
an outlet never needs promoting on its own — the merge already covers it.  Culvert inlets draw
rather than feed, so they are irrelevant here too.
"""
function _uppermost_merge(fp::DynFlowPath)
    best = nothing
    for (q, pos) in fp.merges
        (best === nothing || pos < best[2]) && (best = (q, pos))
    end
    return best
end

"""
    _kill_path!(net, path_id, detached) -> nothing

Path `path_id` has no live tributary, so it dies.  If it targets a trap, that trap's `in_count`
is decremented, and when it loses its last feed it leaves the network too: it is appended to
`detached`, its spill path cleared, and the kill cascades down that path.  A tributary or
domain-exit path is simply removed.

Mutates `net` (tombstones the path, patches `in_count` / `spill_path`) and appends to
`detached`.
"""
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

"""
    _detach_from_host_merges!(net, path_id) -> nothing

Drop the dying tributary `path_id` from the `merges` list of whatever path hosts it.  Mutates
those paths' `merges`.  No trap `in_count` is touched: a tributary feeds a path, not a trap.
"""
function _detach_from_host_merges!(net::DynNetwork, path_id::Int)
    for host in net.flow_paths
        filter!(m -> m[1] != path_id, host.merges)
    end
end

"""
    _tombstone_path!(net, path_id) -> DynFlowPath

Mark path `path_id` dead by replacing its slot with an empty marker
(`departure_point == (0,0)`), leaving every other index intact so references stay valid until
`_compact!` reindexes.  Mutates `net.flow_paths`; live references to the slot are the caller's
to clear beforehand.
"""
_tombstone_path!(net::DynNetwork, path_id::Int) =
    net.flow_paths[path_id] = DynFlowPath(CartesianIndex{2}[], CartesianIndex(0, 0), 0,
        Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[])

"""
    _promote_tributary!(net, path_id, j, q) -> nothing

Re-root: the tributary `Q` in slot `q`, joining at junction position `j`, takes over the
downstream role of the orphaned path `P` in slot `path_id`.

The survivor is `Q` followed by `P` from the junction down, kept in `Q`'s slot since it is `Q`'s
continuation — `Q`'s head trap already spills into `q`, and `P`'s head trap was cleared before
the re-root, so no `spill_path` redirect is needed and `P`'s slot is tombstoned.  `P`'s dead
head stub (cells `1:j-1`) is dropped, the target trap stays fed (`in_count` unchanged), and
injections all sit at positions `>= j`, so none is lost.  Mutates `net`.

!!! note
    @@@ a culvert *inlet* in the dropped stub (position < j) is dropped with it.
"""
function _promote_tributary!(net::DynNetwork, path_id::Int, j::Int, q::Int)
    P, Q = net.flow_paths[path_id], net.flow_paths[q]
    lq = length(Q.cells)
    reb(p) = lq + p - j + 1
    ci = vcat(Q.culvert_inlets,  [(c, reb(p)) for (c, p) in P.culvert_inlets  if p >= j])
    co = vcat(Q.culvert_outlets, [(c, reb(p)) for (c, p) in P.culvert_outlets if p >= j])
    ni = vcat(Q.nbs_inlets,      [(n, reb(p)) for (n, p) in P.nbs_inlets      if p >= j])
    no = vcat(Q.nbs_outlets,     [(n, reb(p)) for (n, p) in P.nbs_outlets     if p >= j])
    mg = vcat(Q.merges,          [(m, reb(p)) for (m, p) in P.merges          if p != j])
    net.flow_paths[q] = DynFlowPath(vcat(Q.cells, P.cells[j:end]),
        Q.departure_point, P.target_trap, ci, co, ni, no, mg)
    _tombstone_path!(net, path_id)
end

# ----------------------------------------------------------------------------
"""
    grow_spill!(net, tstruct, full_traps, trap_id) -> Vector{Int}

Handle trap `trap_id` becoming full and starting to spill: trace its spill downstream and
attach the new path(s)/trap(s), maintaining `in_count` incrementally over just the added
elements (never a full `init_in_counts!` pass).  The mirror of [`detach_spill!`].

# Arguments
- `net::DynNetwork`: network, mutated in place.
- `tstruct`: spillanalysis structure (terrain, spillpoints, footprints).
- `full_traps`: set of currently-full trap indices; must include `trap_id`'s `trap_ix`.
- `trap_id::Int`: local index (in `net.traps`) of the trap that started spilling.

# Returns
The `trap_ix` of the traps newly brought into the dynamic network, for the caller to
migrate from static handling and seed their state.

# Notes
The trace stops at the first already-present trap (its downstream is already represented;
only the connector is attached).  Growth crossing into a *different* component (fusion) is
**not** handled here — the caller (`apply_fill!`) detects the coupling and fuses.
"""
function grow_spill!(net::DynNetwork, tstruct, full_traps, trap_id::Int)
    spoint = tstruct.spillpoints[net.traps[trap_id].trap_ix]
    if spoint.downstream_region_cell == spoint.current_region_cell
        net.traps[trap_id].spill_path = -1        # spills straight out of the domain
        return Int[]
    end

    LI = LinearIndices(tstruct.topography)
    CI = CartesianIndices(tstruct.topography)
    pathmap = Dict{Int,Int}()                     # existing net: cell (linear) -> path index
    for (pi, p) in enumerate(net.flow_paths), c in p.cells
        pathmap[LI[c]] = pi
    end

    np0, nt0 = length(net.flow_paths), length(net.traps)
    _grow_network_from_seed!(net, pathmap, CI[spoint.downstream_region_cell], tstruct, full_traps;
                             departing_trap_ix=trap_id)

    # Incremental in_count: each new path feeds its target trap; each new trap additionally
    # gets the culvert-inlet floor.  (A new path merging into an existing one has target 0 and
    # adds nothing — the merge feeds a path, not a trap.)  Same rule as init_in_counts!, over
    # the added elements only.
    for pi in (np0 + 1):length(net.flow_paths)
        tt = net.flow_paths[pi].target_trap
        tt > 0 && (net.traps[tt].in_count += 1)
    end
    added = Int[]
    for ti in (nt0 + 1):length(net.traps)
        isempty(net.traps[ti].culvert_inlets) || (net.traps[ti].in_count += 1)
        push!(added, net.traps[ti].trap_ix)
    end
    return added
end

# ----------------------------------------------------------------------------
"""
    _index_components(comps, tstruct) -> (trapmap, cellmap)

Global lookup over the live components, rebuilt on demand each event — cheap, cannot go stale,
and doubles as the fusion detector.  Nothing is mutated.

# Returns
- `trapmap`: `trap_ix -> (component index, local trap index)`, to locate an event's trap.
- `cellmap`: occupied cell (linear) `-> component index`, covering flow-path cells, trap
  footprints and NBS footprints.  A grown path cell landing on a *different* component here is
  what triggers a fusion.

Components are disjoint by construction, so every cell and trap resolves to exactly one.
"""
function _index_components(comps::Vector{DynNetwork}, tstruct)
    LI = LinearIndices(tstruct.topography)
    trapmap = Dict{Int, Tuple{Int,Int}}()
    cellmap = Dict{Int, Int}()
    for (ci, net) in enumerate(comps)
        for (li, t) in enumerate(net.traps)
            trapmap[t.trap_ix] = (ci, li)
            for k in tstruct.footprints[t.trap_ix]; cellmap[k] = ci; end
        end
        for p in net.flow_paths, c in p.cells;  cellmap[LI[c]] = ci; end
        for nb in net.nbs, k in nb.footprint;    cellmap[k]     = ci; end
    end
    return trapmap, cellmap
end

"""
    _locate_trap(comps, trap_ix) -> (component index, local trap index) or nothing

Locate `trap_ix` in the component set, or `nothing` if it is not currently dynamic.  An
`O(traps)` scan that builds no cell map — for when only the trap is needed and the full
[`_index_components`](@ref) lookup would be waste.  Nothing is mutated.
"""
function _locate_trap(comps::Vector{DynNetwork}, trap_ix::Int)
    for (ci, net) in enumerate(comps), (li, t) in enumerate(net.traps)
        t.trap_ix == trap_ix && return (ci, li)
    end
    return nothing
end

# All trap_ix currently in the dynamic network (as a set).
_live_trap_ix(comps::Vector{DynNetwork}) = Set(t.trap_ix for net in comps for t in net.traps)

"""
    _fill_subsumes(tstruct, trap_ix, full_traps) -> Bool

True when filling `trap_ix` completes its parent's subtraps, so the siblings must collapse into
the parent supertrap.  That is a hierarchy change the incremental grow cannot express, so the
caller regrows the component instead.  Nothing is mutated.

Testing one level up suffices to trigger: the regrow's tracer cascades any further collapse.
"""
function _fill_subsumes(tstruct, trap_ix, full_traps)
    p = parentof(tstruct, trap_ix)
    p === nothing && return false
    return all(s in full_traps for s in subtrapsof(tstruct, p))
end

# ----------------------------------------------------------------------------
"""
    apply_fill!(comps, tstruct, full_traps, trap_ix) -> Vector{Int}

Apply a not-full → full transition for trap `trap_ix` across the live component set: grow
the owning component, then fuse in any component the growth coupled into.

# Arguments
- `comps::Vector{DynNetwork}`: live components, mutated in place.
- `tstruct`: spillanalysis structure.
- `full_traps`: currently-full trap indices (must already include `trap_ix`).
- `trap_ix::Int`: spillanalysis index of the trap that filled (must be dynamic already).

# Returns
The `trap_ix` of traps newly brought into the dynamic network, for the caller's state handoff.
"""
function apply_fill!(comps::Vector{DynNetwork}, tstruct, full_traps, trap_ix::Int)
    # index the pre-grow components (fusion-detector snapshot); locate the event trap
    trapmap, cellmap = _index_components(comps, tstruct)
    haskey(trapmap, trap_ix) || return Int[]        # not dynamic — nothing to grow
    ci, li = trapmap[trap_ix]
    before = _live_trap_ix(comps)

    if _fill_subsumes(tstruct, trap_ix, full_traps)
        # hierarchy boundary: incremental grow can't collapse the sibling group into the parent,
        # so regrow the owner (the build tracer picks the right supertrap node).
        _fuse_components!(comps, [ci], tstruct, full_traps)
    else
        # incremental grow of the owning component
        owner = comps[ci]
        np0 = length(owner.flow_paths)
        grow_spill!(owner, tstruct, full_traps, li)

        # fusion: a new trap already living elsewhere, or a new path cell on another component
        # (a shared corridor or an NBS footprint), couples that component — regrow the union
        LI = LinearIndices(tstruct.topography)
        coupled = Set{Int}()
        for t in owner.traps
            haskey(trapmap, t.trap_ix) && trapmap[t.trap_ix][1] != ci && push!(coupled, trapmap[t.trap_ix][1])
        end
        for k in (np0 + 1):length(owner.flow_paths), c in owner.flow_paths[k].cells
            cj = get(cellmap, LI[c], ci)
            cj == ci || push!(coupled, cj)
        end
        isempty(coupled) || _fuse_components!(comps, push!(collect(coupled), ci), tstruct, full_traps)
    end

    # newly-dynamic traps (set difference handles both the incremental and regrow paths)
    return collect(setdiff(_live_trap_ix(comps), before))
end

# ----------------------------------------------------------------------------
"""
    apply_unfill!(comps, tstruct, trap_ix) -> Vector{Int}

Apply a full → not-full transition for trap `trap_ix`: detach its spill in the owning
component (compaction runs inside `detach_spill!`; fission is deferred).

# Arguments
- `comps::Vector{DynNetwork}`: live components, mutated in place.
- `tstruct`: spillanalysis structure.
- `trap_ix::Int`: spillanalysis index of the trap that stopped overflowing.

# Returns
The `trap_ix` of traps that left the dynamic network, for handoff back to static handling.
"""
function apply_unfill!(comps::Vector{DynNetwork}, tstruct, trap_ix::Int)
    # locate the trap; if it is not dynamic there is nothing to detach
    loc = _locate_trap(comps, trap_ix)
    loc === nothing && return Int[]
    ci, li = loc
    # detach its spill: cascades downstream and compacts that component in place.  Fission is
    # deferred — the survivor may split into disconnected pieces but stays one correct solve.
    # (A trap draining further, below the level that merged its subtraps, is de-subsumption —
    # a separate drain event handled by `apply_empty!` on the parent node, not here.)
    return detach_spill!(comps[ci], li)
end

# ----------------------------------------------------------------------------
"""
    apply_empty!(comps, tstruct, full_traps, trap_ix) -> Vector{Int}

Apply a de-subsumption: supertrap `trap_ix` has drained below the level that merged its
subtraps into one pool, so the basin must split back into the now separately-tracked children.
Regrow the owning component with the current `full_traps` — the build tracer emits the child
nodes (and `split` re-partitions them if they are no longer connected).

# Arguments
- `comps::Vector{DynNetwork}`: live components, mutated in place.
- `tstruct`: spillanalysis structure.
- `full_traps`: currently-full trap indices (the emptied children already removed).
- `trap_ix::Int`: the supertrap node that de-subsumed.

# Returns
The `trap_ix` of the nodes newly split out (for the caller to distribute the parent's water
across them).
"""
function apply_empty!(comps::Vector{DynNetwork}, tstruct, full_traps, trap_ix::Int)
    loc = _locate_trap(comps, trap_ix)
    loc === nothing && return Int[]        # not a live node
    before = _live_trap_ix(comps)
    _fuse_components!(comps, [loc[1]], tstruct, full_traps)   # regrow → de-subsume (+ re-split)
    return collect(setdiff(_live_trap_ix(comps), before))
end

# ----------------------------------------------------------------------------
"""
    _component_seeds(net) -> Vector{CartesianIndex{2}}

The root seed cells `net` was grown from: the departure points of its source-headed paths, i.e.
those no trap spills into.  Regrowing from these reproduces the whole component.

Culvert outlets and NBS seed cells are among them, since each heads its own path.  The
trap-spill connectors are deliberately omitted — the tracer regrows those from `full_traps`.
Nothing is mutated.
"""
function _component_seeds(net::DynNetwork)
    spill = Set(t.spill_path for t in net.traps if t.spill_path > 0)   # paths a trap spills into
    return [net.flow_paths[k].departure_point
            for k in eachindex(net.flow_paths) if !(k in spill)]
end

# ----------------------------------------------------------------------------
"""
    _regrow(seeds, culverts, nbs, tstruct, full_traps) -> Vector{DynNetwork}

Grow a fresh network from `seeds` carrying the given culverts and NBS, and split it into
connected components.

The core of [`setup_network`](@ref) minus the input validation and the NBS cell precompute —
the placements passed here already carry theirs, and their `id`s are already stamped.  Used to
regenerate fused or hierarchy-collapsed components.  `culverts` and `nbs` are held by reference
on the result rather than copied.
"""
function _regrow(seeds, culverts::Vector{DynCulvert}, nbs::Vector{DynNBSPlacement}, tstruct, full_traps)
    net = DynNetwork(culverts, nbs)
    pathmap = Dict{Int,Int}()
    foreach(s -> _grow_network_from_seed!(net, pathmap, s, tstruct, full_traps), seeds)
    init_in_counts!(net)
    return split_network_into_connected_components(net, tstruct)
end

# ----------------------------------------------------------------------------
"""
    _fuse_components!(comps, idxs, tstruct, full_traps) -> comps

Fuse the components at indices `idxs` — coupled by a grow — into fresh, clean component(s) by
regrowing them from their combined seeds.  The build tracer redoes cell overlaps, duplicate-trap
collapsing and culvert/NBS assignment for free, so no bespoke surgery is needed.

**Mutates and returns `comps`**, deleting the fused slots and appending the regrown result
(normally one component, but the regrow may legitimately split it into several).  Passing a
single index is the idiom for regrowing one component in place.
"""
function _fuse_components!(comps::Vector{DynNetwork}, idxs, tstruct, full_traps)
    seeds    = CartesianIndex{2}[]
    culverts = DynCulvert[]
    nbs      = DynNBSPlacement[]
    for c in idxs
        append!(seeds,    _component_seeds(comps[c]))
        append!(culverts, comps[c].culverts)
        append!(nbs,      comps[c].nbs)
    end
    unique!(culverts); unique!(nbs)              # a coupling element could sit in two of the fused
    regrown = _regrow(seeds, culverts, nbs, tstruct, full_traps)
    deleteat!(comps, sort(collect(idxs)))
    append!(comps, regrown)
    return comps
end
