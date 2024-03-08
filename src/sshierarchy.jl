import Graphs
export sshierarchy!

"""
    sshierarchy!(grid, spillregions, spillpoints, boundaries)

Identify higher-level traps and spillpoints, resulting from merger of lower-level
traps that flow into each other.

The vector of _spillpoints_, and the associated _spillregion boundaries_ will be
updated with the higher-level spillpoints and traps.

The arguments `grid` and `spillregions` will remain unmodified.

The function also returns a _graph_ representing the hierarchy of subtraps and
supertraps, and a _vector_ listing the indices of the lowest-level
spillregions contained within each trap.  For lowest-level traps, this list only
references the trap itself, whereas for higher-level traps, it references all
the low-level traps whose merger leads to the higher-level trap.

# Arguments
- `grid::Matrix{<:Real}`: terrain raster grid with height values.  This
                          input argument remains unmodified.
- `spillregions::Matrix{Int}`: the matrix containing the lowest-level 
                               regions, already computed from the 
                               [`spillregions`](@ref) function.  This
                               input argument remains unmodified.
- `spillpoints::Vector{Spillpoint}`: vector of lowest-level spillpoints
                                     (corresponding to the lowest-level 
                                     regions), as computed by the 
                                     [`spillpoints`](@ref) function.  
                                     This vector will be modified, as it 
                                     will be supplemented with the newly
                                     identified, higher-level spillpoints.
- `boundaries::Vector{Vector{Tuple{Int, Int}}}`: 
            vector containing the boundaries for each of the lowest-level spill
            regions.  This vector will be modified, as it will be extended
            with the boundaries of the newly found higher-level spill regions.
            Each boundary is represented by a vector of integer pairs, where 
            each integer pair are the indices of two neighbor cells in the grid;
            one on each side of the region boundary.
"""
function sshierarchy!(
    grid::Matrix{<:Real},              # remains constant
    spillregions::Matrix{Int},         # remains constant
    spillpoints::Vector,               # will be modified (added to)
    boundaries::Vector{Vector{Tuple{Int, Int}}}) # will be modified (added to)

    cur_num_traps = maximum(spillregions)
    
    # Directed graph keeping track of subtraps and supertraps.  In the
    # beginning, only the bottom-level traps are present
    subtrapgraph = Graphs.SimpleDiGraph(cur_num_traps)

    # identify cycles to merge into higher-level traps, keep doing this until
    # there are no more cycles (we have found the highest-level traps)
    spillgraph = Graphs.SimpleDiGraph(Graphs.Edge.(
        [(i, spillpoints[i].downstream_region) for i in 1:length(spillpoints)]))

    lowest_regions = [[r] for r in 1:length(spillpoints)]

    mergers = Graphs.simplecycles(spillgraph)
    
    while !isempty(mergers)
        Graphs.add_vertices!(subtrapgraph, length(mergers))
        Graphs.add_vertices!(spillgraph, length(mergers))

        # identify spillpoints and update update the subtrap/supertrap graph
        for m in mergers
            cur_num_traps = cur_num_traps + 1

            # identify new spillpoints and lowest-level downstream traps
            spoint = _identify_new_spillpoint!(grid, m, spillregions, boundaries)
            push!(spillpoints, spoint)

            # update the graph respresenting the subtrap hierarchy (subtraphraph)
            for subtrap in m
                Graphs.add_edge!(subtrapgraph, Graphs.Edge(subtrap, cur_num_traps))
            end

            # register the lowest-level regions covered by this trap's spill field
            push!(lowest_regions, vcat(lowest_regions[m]...))
        end

        # insert new traps into spillgraph structure
        new_mergers = Set{Vector{Int}}()
        m_ix = cur_num_traps - length(mergers)
        sibling_connections = Set{Tuple{Int, Int}}()
        for m in mergers
            m_ix = m_ix + 1
            spoint = spillpoints[m_ix]
            
            # determine highest-level downstream traps
            top_dstrap = _identify_highest_level_trap_of(spoint.downstream_region,
                                                         subtrapgraph)
            # replace cycles with new traps in graph
            _merge_traps!(spillgraph, m, m_ix, top_dstrap)

            # check if a new cycle was generated from the merger
            if top_dstrap != 0 && (spillpoints[top_dstrap].elevation == spoint.elevation)
                visited = [m_ix]
                next = top_dstrap
                while next != 0 && spillpoints[next].elevation == spoint.elevation
                    push!(sibling_connections, (visited[end], next))
                    (next in visited) ? break : push!(visited, next)
                    next = _identify_highest_level_trap_of(
                        spillpoints[next].downstream_region, subtrapgraph)
                end
            end
        end
        # determine new mergers by identifying connected siblings
        mergers = _identify_sibling_loops(sibling_connections)
    end
    return subtrapgraph, lowest_regions
end

# ----------------------------------------------------------------------------
function _identify_sibling_loops(connections)
    # remap numbering
    unique_numbers = unique(vcat([x[1] for x in connections],
                                 [x[2] for x in connections]))
    num_nodes = length(unique_numbers)
    num2loc = Dict([(unique_numbers[i], i) for i in 1:length(unique_numbers)])
    
    g = Graphs.SimpleDiGraph([Graphs.SimpleEdge{Int64}((num2loc[c[1]], num2loc[c[2]]))
                              for c in connections])
    # The strongly connected components will be the loops in the graph.
    # It is important that the components are strongly connected for the fill_sequence
    # algorithm to function correctly, in particular its _update_spillgraph! function.
    comps = Graphs.strongly_connected_components(g) 

    # keep only components with more than one trap, and remap back to original numbering
    siblings = [unique_numbers[c] for c in filter(x->length(x) > 1, comps)]
    return siblings
end

# ----------------------------------------------------------------------------
function _merge_traps!(spillgraph, traps, new_trap_ix, top_dstrap)

    # delete all incoming and outgoing edges.  Keep track of the incoming edges,
    # as these will be connected to the new trap.
    all_incoming = Set{Int}()
    for t in traps
        cur_incoming = Set(Graphs.inneighbors(spillgraph, t))
        union!(all_incoming, cur_incoming)
        outgoing = Graphs.outneighbors(spillgraph, t)
        for i in cur_incoming
            # note, we do not remove the vertices themselves, as we do not want
            # renumbering of remaining vertices
            Graphs.rem_edge!(spillgraph, i, t)
        end
        for o in outgoing
            Graphs.rem_edge!(spillgraph, t, o)
        end
    end
    
    # removing incoming nodes that are themselves part of the nodes to be
    # discarded from the list of incoming nodes
    setdiff!(all_incoming, traps)

    # Connect all incoming edges to the new trap, as well as the single outgoing
    for n in all_incoming
        Graphs.add_edge!(spillgraph, (n, new_trap_ix)); # incoming edges
    end
    
    Graphs.add_edge!(spillgraph, (new_trap_ix, top_dstrap)); # outgoing edge
end   

# ----------------------------------------------------------------------------
function _identify_highest_level_trap_of(lowest_trap, subtrapgraph)

    if lowest_trap <= 0 # the out-of-domain region doesn't have a higher level trap
        return 0
    end
    
    cur_trap = lowest_trap

    neighs = Graphs.outneighbors(subtrapgraph, cur_trap)

    while !isempty(neighs)
        @assert(length(neighs)==1, "A subtrap should only have _one_ immediate supertrap.")
        cur_trap = neighs[1]
        neighs = Graphs.outneighbors(subtrapgraph, cur_trap)
    end

    # we climbed as high as we could along the subtrap/supertrap hierarchy
    return cur_trap
end

# ----------------------------------------------------------------------------
function _identify_new_spillpoint!(grid, subtraps, spillregions, boundaries)
    # - `boundaries` will be updated, the other arguments remain constant.
    # NOTE: `subtraps` and `boundaries` are here made with reference to the
    # highest-level working spillgraph, while `spillregions` is made with reference
    # to the lowest-level spill regions.

    # compute the outer boundary of the new region.  We do that by merging all
    # the boundaries of the subregions, and removing those parts that are shared
    # ("internal" to the new region)

    # @@ by not copying, the boundaries for these subregions will no longer be
    # correct, but they are presumably not needed anymore
    outer_bnd = [boundaries[s] for s in subtraps]; 
                                                   
    # remove shared parts
    for i = 1:length(subtraps)-1
        for j = i+1:length(subtraps)
            bnd1 = boundaries[subtraps[i]]
            bnd2 = boundaries[subtraps[j]]

            outer_bnd[i] = setdiff(outer_bnd[i], [(x[2], x[1]) for x in bnd2])
            outer_bnd[j] = setdiff(outer_bnd[j], [(x[2], x[1]) for x in bnd1])
        end
    end
    
    boundary = vcat(outer_bnd[:]...)

    zvals1 = grid[[x[1] for x in boundary]]
    zvals2 = grid[[x[2] for x in boundary]]

    zvals = max.(zvals1, zvals2)

    push!(boundaries, boundary)    
    if isempty(zvals)
        # No outer boundary of trap.  This may happen if the merged traps are
        # enclosed within building footprints.
        return Spillpoint();
    end

    ix = argmin(zvals)
    
    # spillpoint consists or neighbor region index, spillpoint cell in current region,
    # spillpoint cell in other region, and zvalue.

    # check if we are spilling out of the boundary, which is detected by both
    # entries of the boundary element pointing to the same (inside) cell,
    # instead of pointing to one cell in the region and one in the neighbor
    # region
    reg_ix = boundary[ix][1] == boundary[ix][2] ? 0 : spillregions[boundary[ix][2]]
    
    return Spillpoint(reg_ix, boundary[ix][1], boundary[ix][2], zvals[ix])
end
