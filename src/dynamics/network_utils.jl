import Graphs

export split_network_into_connected_components

# ----------------------------------------------------------------------------
"""
    split_network_into_connected_components(net::DynNetwork, tstruct) -> Vector{DynNetwork}

Split a monolithic [`DynNetwork`](@ref) into independently-solvable connected
components.  Two elements share a component iff a flow link, a merge, a culvert, or
an NBS couples them, so each component can be solved without reference to the
others.  Every returned network has its cross-references (target traps, spill
paths, merges, culverts, NBS) remapped to local 1-based indices.
"""
function split_network_into_connected_components(net::DynNetwork, tstruct)
    np = length(net.flow_paths)
    nbs_coupled = _nbs_coupled_nodes(net, tstruct, np)
    g  = _coupling_graph(net, np, nbs_coupled)
    return [_build_subnetwork(net, nodes, np, nbs_coupled)
            for nodes in Graphs.connected_components(g)]
end

"""
    _coupling_graph(net, np, nbs_coupled) -> Graphs.SimpleGraph

Undirected graph whose connected components are the independently-solvable subnetworks.  Nodes
are the flow paths (`1:np`) then the traps (`np+1 : np+nt`); an edge means the two elements are
hydraulically coupled and must be solved together.

Edges come from flow links (path to its target trap, trap to its spill path), merges, and the
culverts and NBS placements that tie otherwise-separate elements together.  Nothing is mutated.
"""
function _coupling_graph(net::DynNetwork, np::Int, nbs_coupled)
    g = Graphs.SimpleGraph(np + length(net.traps))
    tnode(t) = np + t
    for (p, fp) in enumerate(net.flow_paths)
        fp.target_trap > 0 && Graphs.add_edge!(g, p, tnode(fp.target_trap))
        for (m, _) in fp.merges; Graphs.add_edge!(g, p, m); end
    end
    for (t, tr) in enumerate(net.traps)
        tr.spill_path > 0 && Graphs.add_edge!(g, tnode(t), tr.spill_path)
    end
    for nodes in _culvert_owner_nodes(net, np); _link_all!(g, nodes); end
    for nodes in nbs_coupled;                   _link_all!(g, nodes); end
    return g
end

# Connect a set of nodes into one component (a chain suffices for connectivity).
_link_all!(g, nodes) = for i in 2:length(nodes); Graphs.add_edge!(g, nodes[1], nodes[i]); end

"""
    _culvert_owner_nodes(net, np) -> Vector{Vector{Int}}

For each culvert, the coupling-graph nodes hosting its inlet or outlet.  A culvert moves water
between its ends, so both sides must share a component — the caller links each list into one.
Nothing is mutated.
"""
function _culvert_owner_nodes(net::DynNetwork, np::Int)
    owners = [Int[] for _ in net.culverts]
    for (p, fp) in enumerate(net.flow_paths), (c, _) in vcat(fp.culvert_inlets, fp.culvert_outlets)
        push!(owners[c], p)
    end
    for (t, tr) in enumerate(net.traps), c in vcat(tr.culvert_inlets, tr.culvert_outlets)
        push!(owners[c], np + t)
    end
    return owners
end

"""
    _nbs_coupled_nodes(net, tstruct, np) -> Vector{Vector{Int}}

For each NBS placement, the coupling-graph nodes causally tied to it and so required to share
its component: its outlet paths, the paths leaving its footprint, the paths entering it, and
its internal-accumulation trap.

An outflow cell seeds a path, so it is matched against the path's `departure_point`; an inflow
cell is crossed mid-path, so it is matched anywhere in `cells`.  Nothing is mutated.
"""
function _nbs_coupled_nodes(net::DynNetwork, tstruct, np::Int)
    coupled = [Int[] for _ in net.nbs]
    outflow = [Set(n.footprint_outflow_cells) for n in net.nbs]
    inflow_of = Dict{CartesianIndex{2},Int}()          # inflow cell -> the NBS it feeds
    for (n, nb) in enumerate(net.nbs), c in nb.footprint_inflow_cells
        inflow_of[c] = n
    end
    for (p, fp) in enumerate(net.flow_paths)
        for (n, _) in fp.nbs_outlets; push!(coupled[n], p); end
        for n in eachindex(net.nbs)
            fp.departure_point in outflow[n] && push!(coupled[n], p)
        end
        for cell in fp.cells
            n = get(inflow_of, cell, 0); n > 0 && push!(coupled[n], p)
        end
    end
    _add_accumulation_traps!(coupled, net, tstruct, np)
    return [unique!(v) for v in coupled]
end

"""
    _add_accumulation_traps!(coupled, net, tstruct, np) -> nothing

Unite each NBS with the network trap holding its internal-accumulation cells — water ponding
inside the footprint is that trap's, so the two must solve together.

Cells map to lowest-level regions by direct `regions` lookup (dedup absorbs flat, many-celled
bottoms), and each region to the single network trap in its supertrap hierarchy; that
uniqueness is asserted.  **Mutates `coupled`**, appending the trap node to each NBS's list.
"""
function _add_accumulation_traps!(coupled, net::DynNetwork, tstruct, np::Int)
    tnode = Dict(t.trap_ix => np + l for (l, t) in enumerate(net.traps))
    for (n, nb) in enumerate(net.nbs)
        for r in unique(tstruct.regions[c] for c in nb.internal_accumulation_cells)
            covering = [tnode[st] for st in tstruct.supertraps_of[r] if haskey(tnode, st)]
            @assert length(covering) == 1 "NBS accumulation region $r maps to $(length(covering)) network traps (expected 1)"
            push!(coupled[n], covering[1])
        end
    end
end

"""
    _build_subnetwork(net, nodes, np, nbs_coupled) -> DynNetwork

Rebuild the component spanning `nodes` as a standalone [`DynNetwork`](@ref), remapping every
cross-reference — target traps, spill paths, merges, culverts, NBS — from `net`'s indices to
local 1-based ones, so the result is solvable without reference to `net`.

`net` is not mutated, though the new paths alias its `cells` vectors by reference.
"""
function _build_subnetwork(net::DynNetwork, nodes::Vector{Int}, np::Int, nbs_coupled)
    nodeset = Set(nodes)
    pids = sort!([n      for n in nodes if n <= np])
    tids = sort!([n - np for n in nodes if n >  np])
    cvids  = _referenced_culverts(net, pids, tids)
    nbsids = [n for n in eachindex(net.nbs) if any(o in nodeset for o in nbs_coupled[n])]

    pmap   = _local_index(pids)
    tmap   = _local_index(tids)
    cvmap  = _local_index(cvids)
    nbsmap = _local_index(nbsids)

    paths = [_localize_path(net.flow_paths[g], tmap, cvmap, nbsmap, pmap) for g in pids]
    traps = [_localize_trap(net.traps[g], pmap, cvmap) for g in tids]
    return DynNetwork(paths, traps, net.culverts[cvids], net.nbs[nbsids])
end

_local_index(global_ids) = Dict(g => l for (l, g) in enumerate(global_ids))

"""
    _referenced_culverts(net, pids, tids) -> Vector{Int}

The sorted, deduplicated ids of every culvert with an inlet or outlet on the given paths
(`pids`) or traps (`tids`) — the culverts a component must carry.  Nothing is mutated.
"""
function _referenced_culverts(net::DynNetwork, pids, tids)
    cv = Int[]
    for g in pids, (c, _) in vcat(net.flow_paths[g].culvert_inlets, net.flow_paths[g].culvert_outlets)
        push!(cv, c)
    end
    for g in tids, c in vcat(net.traps[g].culvert_inlets, net.traps[g].culvert_outlets)
        push!(cv, c)
    end
    return sort!(unique!(cv))
end

# Sentinels pass through unchanged: 0 = none, -1 = out-of-domain (trap spill only).
_relabel(ix, m) = ix <= 0 ? ix : m[ix]

"""
    _localize_path(fp, tmap, cvmap, nbsmap, pmap) -> DynFlowPath

Copy of `fp` with every reference relabelled from global to component-local indices, via the
supplied trap / culvert / NBS / path maps.  `cells` and `departure_point` are terrain
coordinates and carry over unchanged.  Nothing is mutated.
"""
function _localize_path(fp::DynFlowPath, tmap, cvmap, nbsmap, pmap)
    DynFlowPath(fp.cells, fp.departure_point,
                _relabel(fp.target_trap, tmap),
                [(cvmap[c],  pos) for (c, pos) in fp.culvert_inlets],
                [(cvmap[c],  pos) for (c, pos) in fp.culvert_outlets],
                [(nbsmap[n], pos) for (n, pos) in fp.nbs_inlets],
                [(nbsmap[n], pos) for (n, pos) in fp.nbs_outlets],
                [(pmap[m],   j)   for (m, j)   in fp.merges])
end

"""
    _localize_trap(tr, pmap, cvmap) -> DynTrap

Copy of `tr` with its spill path and culvert ids relabelled to component-local indices.
`trap_ix` is a `TrapStructure` index and stays global; `in_count` carries over because a trap's
feeders all live in its own component.  Nothing is mutated.
"""
function _localize_trap(tr::DynTrap, pmap, cvmap)
    lt = DynTrap(tr.trap_ix,
                 _relabel(tr.spill_path, pmap),
                 [cvmap[c] for c in tr.culvert_inlets],
                 [cvmap[c] for c in tr.culvert_outlets])
    lt.in_count = tr.in_count   # component-invariant (a trap's feeders are all in its component)
    return lt
end
