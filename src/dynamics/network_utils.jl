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
    nbs_owners = _nbs_owner_nodes(net, tstruct, np)
    g  = _coupling_graph(net, np, nbs_owners)
    return [_build_subnetwork(net, nodes, np, nbs_owners)
            for nodes in Graphs.connected_components(g)]
end

# Undirected graph over path nodes (1:np) and trap nodes (np+1:np+nt); an edge
# means the two elements must live in the same component.
function _coupling_graph(net::DynNetwork, np::Int, nbs_owners)
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
    for nodes in nbs_owners;                    _link_all!(g, nodes); end
    return g
end

# Connect a set of nodes into one component (a chain suffices for connectivity).
_link_all!(g, nodes) = for i in 2:length(nodes); Graphs.add_edge!(g, nodes[1], nodes[i]); end

# For each culvert, the path/trap nodes hosting its inlet or outlet (both sides
# must share a component, so the culvert unites its owners).
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

# For each NBS, the nodes coupled to it (sharing its mass balance): paths carrying
# one of its outlets (emission), paths emanating from an outflow-boundary cell
# (deduction, detected as the path's first cell), and the trap holding its
# internal-accumulation cells (deduction).
function _nbs_owner_nodes(net::DynNetwork, tstruct, np::Int)
    owners  = [Int[] for _ in net.nbs]
    outflow = [Set(n.footprint_outflow_cells) for n in net.nbs]
    for (p, fp) in enumerate(net.flow_paths)
        for (n, _) in fp.nbs_outlets; push!(owners[n], p); end
        isempty(fp.cells) && continue
        for n in eachindex(net.nbs)
            first(fp.cells) in outflow[n] && push!(owners[n], p)
        end
    end
    _add_accumulation_traps!(owners, net, tstruct, np)
    return [unique!(o) for o in owners]
end

# Unite each NBS with the network trap holding its internal-accumulation cells.
# Cells map to lowest-level regions (direct `regions` lookup; dedup absorbs flat,
# many-celled bottoms), each region to the one network trap in its supertrap
# hierarchy.
function _add_accumulation_traps!(owners, net::DynNetwork, tstruct, np::Int)
    tnode = Dict(t.trap_ix => np + l for (l, t) in enumerate(net.traps))
    for (n, nb) in enumerate(net.nbs)
        for r in unique(tstruct.regions[c] for c in nb.internal_accumulation_cells)
            covering = [tnode[st] for st in tstruct.supertraps_of[r] if haskey(tnode, st)]
            @assert length(covering) == 1 "NBS accumulation region $r maps to $(length(covering)) network traps (expected 1)"
            push!(owners[n], covering[1])
        end
    end
end

# Rebuild one component as a standalone DynNetwork, remapping every cross-reference
# (target traps, spill paths, merges, culverts, NBS) to local 1-based indices.
function _build_subnetwork(net::DynNetwork, nodes::Vector{Int}, np::Int, nbs_owners)
    nodeset = Set(nodes)
    pids = sort!([n      for n in nodes if n <= np])
    tids = sort!([n - np for n in nodes if n >  np])
    cvids  = _referenced_culverts(net, pids, tids) # all culverts referenced by the component's paths and traps
    # an NBS belongs here iff any element it is coupled to (emission / outflow /
    # accumulation, from _nbs_owner_nodes) is in this component
    nbsids = [n for n in eachindex(net.nbs) if any(o in nodeset for o in nbs_owners[n])]

    pmap   = _local_index(pids)
    tmap   = _local_index(tids)
    cvmap  = _local_index(cvids)
    nbsmap = _local_index(nbsids)

    paths = [_localize_path(net.flow_paths[g], tmap, cvmap, nbsmap, pmap) for g in pids]
    traps = [_localize_trap(net.traps[g], pmap, cvmap) for g in tids]
    return DynNetwork(paths, traps, net.culverts[cvids], net.nbs[nbsids])
end

_local_index(global_ids) = Dict(g => l for (l, g) in enumerate(global_ids))

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

function _localize_path(fp::DynFlowPath, tmap, cvmap, nbsmap, pmap)
    DynFlowPath(fp.cells,
                _relabel(fp.target_trap, tmap),
                [(cvmap[c],  pos) for (c, pos) in fp.culvert_inlets],
                [(cvmap[c],  pos) for (c, pos) in fp.culvert_outlets],
                [(nbsmap[n], pos) for (n, pos) in fp.nbs_outlets],
                [(pmap[m],   j)   for (m, j)   in fp.merges])
end

function _localize_trap(tr::DynTrap, pmap, cvmap)
    DynTrap(tr.trap_ix,
            _relabel(tr.spill_path, pmap),
            [cvmap[c] for c in tr.culvert_inlets],
            [cvmap[c] for c in tr.culvert_outlets])
end
