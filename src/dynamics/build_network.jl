# NOTES TO SELF
# - [x] DynNBSPlacement moved to elements.jl (renamed from NBSPlacement).  DynNBS still present;
#       remove it once the distributor re-impl no longer needs the placement_ix overlay handle.
# - [ ] Contract that all NBS footprints should have zero infiltration (to avoid special handling in fill_sequence)
# - [ ] This file is not yet included in the module (src/SurfaceWaterIntegratedModeling.jl); wire it in
#       once the path-construction (culvert/nbs tuples + merges arg) and distributor are ready.


# ----------------------------------------------------------------------------
function setup_network(tstruct, full_traps;
                       dyn_coords::Vector{CartesianIndex{2}}=CartesianIndex{2}[],
                       culverts::Vector{DynCulvert}=DynCulvert[],
                       nbs::Vector{DynNBSPlacement}=DynNBSPlacement[])

    # Input validation (no culverts or dyn_coords inside NBS footprints)
    _validate_network_inputs(tstruct, dyn_coords, culverts, nbs)

    # Compute inflow and outflow cells for each NBS footprint
    _compute_nbs_inflow_outflow_cells!(nbs, tstruct)
    
    # seeds for tracking dynamic flows are the dyn_coords, the culvert outlet points,
    # and all NBS outflow points
    seeds = union(Set(dyn_coords),
                  Set([culvert.outlet for culvert in culverts]),
                  Set(vcat([n.outlets for n in nbs]...)),
                  Set(vcat([n.internal_accumulation_cells for n in nbs]...)),
                  Set(vcat([n.footprint_outflow_cells for n in nbs]...)))

    # Generate map to keep track of all existing network elements on terrain 
    pathmap = Dict{Int, Int}() # terrain cell -> flow path index
    
    # Track all seeds and build up corresponding network of paths and traps
    monolithic_network = DynNetwork(culverts, nbs)
    foreach(s -> _grow_network_from_seed!(monolithic_network,
                                          pathmap, s, tstruct, full_traps), seeds)

    # Split network in connected components
    networks = _split_network_into_connected_components(monolithic_network, tstruct)

    return networks
end

# ----------------------------------------------------------------------------
function _map_inverse_flow(flowgraph)
    inv_flow = Dict{Int, Vector{Int}}()
    for cell in 1:Graphs.nv(flowgraph)
        ds = Graphs.outneighbors(flowgraph, cell)
        for d ∈ ds
            inv_flow[d] = get(inv_flow, d, Int[]) # initialize if not present
            push!(inv_flow[d], cell)
        end
    end
    return inv_flow
end

# ----------------------------------------------------------------------------
function _compute_nbs_inflow_outflow_cells!(nbs, tstruct)
    inv_flow = _map_inverse_flow(tstruct.flowgraph)
    CI = CartesianIndices(tstruct.topography)
    LI = LinearIndices(tstruct.topography)
    tbottoms = Set(LI[tstruct.trap_bottoms])
    
    for n in nbs
        n.footprint_outflow_cells = CI[_footprint_outflow_cells(tstruct, n.footprint)]
        n.footprint_inflow_cells = CI[_footprint_inflow_cells(inv_flow, n.footprint)]
        n.internal_accumulation_cells = CI[collect(intersect(Set(n.footprint), tbottoms))]
    end
end

# ----------------------------------------------------------------------------
function _trace_to_next_trap(tstruct, start_cell, full_traps)

    path = _trace_path(tstruct, start_cell)
    end_reg = tstruct.regions[path[end]]

    if end_reg <= 0
        return path, 0
    end

    # Determine the end trap.  This is the same as the uppermost full trap,
    # unless all its siblings are also full, in which case we go up to the next
    # supertrap.
    supertraps = tstruct.supertraps_of[end_reg]
    highest_full_supertrap = findlast(t -> t in full_traps, supertraps)
    if isnothing(highest_full_supertrap)
        # no full traps in the supertrap hierarchy, so we return the lowest-level trap
        return path, end_reg
    end
    
    # if this is not the highest level trap, and all its siblings are full, return the
    # trap above it.  Otherwise, return the highest full trap.
    target_trap =
        (highest_full_supertrap < length(supertraps) &&
        all(s -> s in full_traps, subtrapsof(tstruct, supertraps[highest_full_supertrap + 1]))) ?
        supertraps[highest_full_supertrap + 1] :
        supertraps[highest_full_supertrap]

    # remove trap footprint cells from path, since they are not part of the flow path. 
    path = filter(c -> !(c in tstruct.footprints[target_trap]), path)
        
    return path, target_trap
end

# ----------------------------------------------------------------------------
function _grow_network_from_seed!(network, pathmap, seed::CartesianIndex{2},
                                 tstruct, full_traps)
    LI = LinearIndices(tstruct.topography)
    CI = CartesianIndices(tstruct.topography)
    terminus = CartesianIndex(0, 0)
    departing_trap_ix = 0 # if the path exited a previous trap, this is the trap index
    while seed != terminus
        path, trap_ix = _trace_to_next_trap(tstruct, LI[seed], full_traps)

        trap_local = 0
        # registering trap if it exists
        if trap_ix > 0
            # switch 'trap' from index into tstruct to index into network.traps
            trap_local = findfirst(dtrap -> dtrap.trap_ix == trap_ix, network.traps)
            if isnothing(trap_local)
                # trap_ix is not yet in network.traps, so add it
                cv_in, cv_out, nbs_out =
                    _intersecting_culverts_and_nbs_outlets(CI[tstruct.footprints[trap_ix]], network.culverts, network.nbs)
                push!(network.traps, DynTrap(trap_ix, 0, cv_in, cv_out)) # @@ handle spill_path index later
                trap_local = length(network.traps)
            end
        end
        # registering path if it exists
        if !isempty(path)
            isect_path, isect_cell =
                _update_pathmap!(pathmap, path, length(network.flow_paths)+1) # may truncate path

            # add newly constructed path to network, with intersection culverts, if any
            cv_in, cv_out, nbs_out = _intersecting_culverts_and_nbs_outlets(CI[path], network.culverts, network.nbs)
            target_trap = isect_path > 0 ? 0 : trap_local
            push!(network.flow_paths, DynFlowPath(CI[path], target_trap, cv_in, cv_out, nbs_out))

            # registering where this path exited, if there is a departing trap
            (departing_trap_ix > 0) && (network.traps[departing_trap_ix].spill_path = length(network.flow_paths))
            
            if isect_path > 0
                # path merged with another path, so we need to register the
                # merge in the other path's merges list
                isect_pt = findfirst(network.flow_paths[isect_path].cells .== CI[isect_cell])
                @assert isect_pt > 0 "Intersection point should be in the other path's cells."
                push!(network.flow_paths[isect_path].merges,
                      (length(network.flow_paths), isect_pt))
            end
        end

        # if path ended in a full trap, it is spilling and we should keep on
        # tracing.  Otherwise, we are done.
        departing_trap_ix = trap_local

        if trap_ix == 0
            seed = terminus
        else
            spoint = tstruct.spillpoints[trap_ix]
            seed = (trap_ix ∈ full_traps) && spoint.downstream_region_cell != spoint.current_region_cell ?
                CI[spoint.downstream_region_cell] : terminus
        end
        
    end
end

# ----------------------------------------------------------------------------
function _update_pathmap!(pathmap, path::Vector{Int}, path_ix)
    # check if path intersects with any existing NBS or culvert
    isect_path, isect_cell = 0, 0
    for (ix, cell) in enumerate(path)
        if haskey(pathmap, cell)
            # truncate path at intersection point
            @assert (ix > 1) "The first point in the path should never be an intersection point."
            #path = path[1:ix-1]
            splice!(path, ix:lastindex(path))  # Mutate `path` by removing elements from `ix` to the end
            isect_path, isect_cell = pathmap[cell], cell
        end
    end
    # register path in pathmap
    for cell in path
        pathmap[cell] = path_ix
    end
    
    return isect_path, isect_cell
end

# ----------------------------------------------------------------------------
function _intersecting_culverts_and_nbs_outlets(footprint, culverts, nbs)
    cv_in = Int[]
    cv_out = Int[]
    nbs_out = Int[]
    for (ix, culvert) in enumerate(culverts)
        if culvert.inlet in footprint
            push!(cv_in, ix)
        end
        if culvert.outlet in footprint
            push!(cv_out, ix)
        end
    end
    for (ix, n) in enumerate(nbs)
        if any(outlet in footprint for outlet in n.outlets)
            push!(nbs_out, ix)
        end
    end
    
    return cv_in, cv_out, nbs_out
end

# ----------------------------------------------------------------------------
# Split the monolithic network into independently-solvable pieces.  Two elements
# share a component iff a flow link, a merge, a culvert, or an NBS couples them,
# so each component can be solved without reference to the others.
#
function _split_network_into_connected_components(net::DynNetwork, tstruct)
    np = length(net.flow_paths)
    g  = _coupling_graph(net, np, tstruct)
    return [_build_subnetwork(net, nodes, np) for nodes in Graphs.connected_components(g)]
end

# Undirected graph over path nodes (1:np) and trap nodes (np+1:np+nt); an edge
# means the two elements must live in the same component.
function _coupling_graph(net::DynNetwork, np::Int, tstruct)
    g = Graphs.SimpleGraph(np + length(net.traps))
    tnode(t) = np + t
    for (p, fp) in enumerate(net.flow_paths)
        fp.target_trap > 0 && Graphs.add_edge!(g, p, tnode(fp.target_trap))
        for (m, _) in fp.merges; Graphs.add_edge!(g, p, m); end
    end
    for (t, tr) in enumerate(net.traps)
        tr.spill_path > 0 && Graphs.add_edge!(g, tnode(t), tr.spill_path)
    end
    for nodes in _culvert_owner_nodes(net, np);      _link_all!(g, nodes); end
    for nodes in _nbs_owner_nodes(net, tstruct, np); _link_all!(g, nodes); end
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
function _build_subnetwork(net::DynNetwork, nodes::Vector{Int}, np::Int)
    pids = sort!([n      for n in nodes if n <= np])
    tids = sort!([n - np for n in nodes if n >  np])
    cvids  = _referenced_culverts(net, pids, tids)
    nbsids = _referenced_nbs(net, pids)

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

function _referenced_nbs(net::DynNetwork, pids)
    ns = Int[]
    for g in pids, (n, _) in net.flow_paths[g].nbs_outlets; push!(ns, n); end
    return sort!(unique!(ns))
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

# ----------------------------------------------------------------------------
function _validate_network_inputs(tstruct, dyn_coords, culverts, nbs)
    # Check that no culvert or dyn_coord is inside an NBS footprint
    nbs_footprints = vcat([n.footprint for n in nbs]...)
    LI = LinearIndices(tstruct.topography)
    for culvert in culverts
        if LI[culvert.inlet] in nbs_footprints || LI[culvert.outlet] in nbs_footprints
            error("Culvert inlet or outlet is inside an NBS footprint.")
        end
    end
    for coord in dyn_coords
        if LI[coord] in nbs_footprints
            error("Dynamic flow coordinate is inside an NBS footprint.")
        end
    end
end

# ----------------------------------------------------------------------------
function _downstream_cell(flowgraph, ix)
    ds = Graphs.outneighbors(flowgraph, ix)
    @assert length(ds) <= 1
    return isempty(ds) ? -1 : ds[1]
end

# ----------------------------------------------------------------------------
function _footprint_outflow_cells(tstruct, footprint::Vector{Int})
    flowgraph = tstruct.flowgraph
    outflow_cells = Set{Int}()
    for cell in footprint
        ds = _downstream_cell(flowgraph, cell)
        if ds != -1 && !(ds in footprint)
            push!(outflow_cells, ds)
        end
    end
    # if there are no outflow cells, that is OK, but if upper NBS layers do not
    # have other, explicitly given outlets, this NBS system will become
    # submerged with no runoff.
    if isempty(outflow_cells)
        @warn "NBS placement footprint has no outflow cells."
    end
    return Vector(outflow_cells)
end

# ----------------------------------------------------------------------------
function _footprint_inflow_cells(inv_flow::Dict{Int, Vector{Int}},
                                 footprint::Vector{Int})
    all_inflow_cells = vcat([inv_flow[cell] for cell in footprint]...)
    return Vector(setdiff(Set(all_inflow_cells), Set(footprint)))
end

