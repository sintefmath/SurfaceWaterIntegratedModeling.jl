import Graphs

export watercourses, saturated_spillgraph

"""
    watercourses(tstruct, full_traps, precipitation, infiltration)

Compute the total intensity of flow of water across all non-lake grid cells in
the terrain.

Determine the areas containing high amounts of flowing water, and return a
grid with the corresponding values.

# Arguments
- `tstruct::TrapStructure`: a TrapStructure object containing all the information from 
                            the topographical analysis, as returned by the
                            `spillanalysis` function.
- `full_traps::Vector{Bool}`: A vector with one entry per (sub)trap, indicating 
                              whether this trap has been filled yet or not.  Only full
                              traps will overflow and forward the water from upstream.
- `precipitation::Matrix{<:Real}`: A grid expressing precipitation rate per grid cell.
                                    If left empty, it will be substituted by a grid 
                                    filled with ones.
- `infiltration::Matrix{<:Real}`: A grid expressing the maximum infiltration rate 
                                  per grid cell.  If left empty, it will be substituted
                                  by a grid filled with zeros.

# Returns
- `runoff::Matrix{Float64}`: Grid expressing infiltration excess runoff rate 
                             (positive values) 
                             or remaining infiltration capacity (negative values)
- `region_accum::Vector{Float64}`: Vector with one entry per spill region, indicating
                          net water accumulation rate for that region (also including
                          water spilling into the region from upstream filled traps). 
                          Note that if the region is associated with a filled trap, 
                          the accumulation value is set to zero, unless it spills into
                          a filled 'sibling' trap.  The logic is in any case such that 
                          the accumulation of any still unfilled trap should equal
                          the sum of accumulation of its regions.
- `offregion_runoff::Float64`: Total rate of water flowing off the domain
- `used_infiltration::Float64`: Total infiltration rate across the terrain
"""
function watercourses(tstruct::TrapStructure{<:Real},
                      full_traps::Vector{Bool};
                      precipitation::Union{Matrix{<:Real}, Real} = 1.0, 
                      infiltration::Union{Matrix{<:Real}, Real} = 0.0)

    # expand `precipitation` and `infiltration` to matrices if necessary
    gridres = size(tstruct.topography)
    (typeof(precipitation) <: Real) && (precipitation = fill(precipitation, gridres...))
    (typeof(infiltration) <: Real) && (infiltration = fill(infiltration, gridres...))
    
    # The following three variables constitute the return values
    runoff = precipitation - infiltration;
    region_accum = zeros(Float64, maximum(tstruct.regions))
    offregion_runoff = Ref(0.0)
    
    # Compute basic flow field intensity, as if all traps were empty
    g = compute_spillfield_graph(tstruct.spillfield)
    sortedg = Graphs.topological_sort_by_dfs(g) # sort nodes to avoid invalidating
                                                # earlier work
    for cur_node in sortedg

        if runoff[cur_node] >= 0
            # there is infiltration excess flow across this node.  Propagate downstream
            ds_node = Graphs.outneighbors(g, cur_node)
            
            if !isempty(ds_node)
                @assert length(ds_node) == 1
                ds_node = ds_node[1]
                runoff[ds_node] += runoff[cur_node]
            else
                # we are at a bottom point.  Register accumulated runoff to
                # the corresponding spill region
                reg = tstruct.regions[cur_node]
                if reg > 0
                    region_accum[reg] += runoff[cur_node]
                else
                    offregion_runoff[] += runoff[cur_node]
                end
            end
        end
    end

    # Adjust flow according to traps that are full (if any)
    topmost_full_traps = any(full_traps) ? 
        filter(x -> !any(full_traps[Graphs.outneighbors(tstruct.agglomerations, x)]),
               findall(full_traps)) : []

    if !isempty(topmost_full_traps)
        # Adjust flow according to the traps that are full

        # we sort these traps according to height to avoid invalidating earlier
        # work as we process each trap in turn
        order = _traps_sortorder(topmost_full_traps, full_traps, tstruct)

        for tr in order
            output, spoint = _compute_trap_output(tstruct, tr, runoff, region_accum)
            _update_runoff!(runoff, region_accum, offregion_runoff,
                            spoint, output, tstruct.regions, g)
            region_accum[tstruct.lowest_subtraps_for[tr]] .= 0.0
        end
    end

    
    # compute total infiltration
    remaining_infil_capacity = max.(-runoff, 0.0);
    for tt in topmost_full_traps
        remaining_infil_capacity[tstruct.footprints[tt]] .= 0.0;
    end
    used_infiltration = sum(infiltration - remaining_infil_capacity)
    
    return runoff, region_accum, offregion_runoff[], used_infiltration

end

# ----------------------------------------------------------------------------
"""
    saturated_spillgraph(tstruct)

Create the spillgraph corresponding to the provided trapping structure, considering
that all traps are already filled.

The result will be returned as a `SimpleDiGraph` showing which traps are
spilling into which traps or regions.  

The trapping structure `tstruct` may have regions with negative indices
(referring to regions that spill out of the terrain).  Since `SimpleDiGraph`
does not handle nodes with negative indices, a remapping is necessary between
the indices used in `tstruct` and the set 1:N used in the returned graph.
These mappings are returned as the second and third return value:
- The second return value is a `Vector{Int}` mapping from 1:N to the indices 
  in `tstruct`
- The third return value is a `Dict` mapping from the indices in `tstruct` to 
  1:N.
"""
function saturated_spillgraph(tstruct::TrapStructure{<:Real})

    function target(trap)
        parent = Graphs.outneighbors(tstruct.agglomerations, trap)
        dsreg = tstruct.spillpoints[trap].downstream_region
        return !isempty(parent) ? parent[1] : dsreg
    end
    
    edges = [(i, target(i)) for i in 1:length(tstruct.spillpoints)]

    return _make_simple_graph(edges)
end

# ----------------------------------------------------------------------------
function _traps_sortorder(top_full_traps, is_full, tstruct)

    toptrapmap = collect(1:length(is_full)) # default is to map to itself
    for tt in top_full_traps
        # all subtraps of
        toptrapmap[tstruct.lowest_subtraps_for[tt]] .= tt
    end

    find_targets = sources ->
        [x > 0 ? toptrapmap[x] : -1 for x in
             [tstruct.spillpoints[y].downstream_region for y in sources]]
    
    # first, we identify all "sibling" subtraps spilling into each other.
    # These do not need to be processed, so we remove them and make another,
    # acyclic graph
    edges = collect(zip(top_full_traps, find_targets(top_full_traps)))
    g, vmap = _make_simple_graph(edges) # this graph may contain cycles
    cyclenodes = unique(vcat(Graphs.simplecycles(g)...))
    top_traps_nocycles = setdiff(top_full_traps, vmap[cyclenodes])
    
    # we can now generate a proper, acyclic graph that can be topologically sorted
    edges = collect(zip(top_traps_nocycles, find_targets(top_traps_nocycles)))
    if isempty(edges)
        # no traps to process
        return Vector{Int64}([])
    end
    g, vmap = _make_simple_graph(edges)

    # sort the full top traps so that processing of later traps do not interfere
    # with the processing of earlier.  We do not need to process leaf nodes
    # (which may not be part of top_full_traps anyway)
    sortedg = Graphs.topological_sort_by_dfs(g)
    is_leaf = [isempty(Graphs.outneighbors(g, x)) for x in sortedg]
    return vmap[sortedg[.!is_leaf]]
end

# ----------------------------------------------------------------------------    
function _compute_trap_output(tstruct::TrapStructure{<:Real},
                             trap_ix::Int64,
                             runoff::Matrix{<:Real}, 
                             region_accum::Vector{<:Real})
    # the runoff from a trap is the sum of the runoff entering all its bottommost
    # subtraps, less the remaining infiltration capacity within the footprint
    # of the trap (since this function should only be called on filled traps, and
    # the footprint of the trap is thus submerged).

    output =
        sum(region_accum[tstruct.lowest_subtraps_for[trap_ix]]) -
        sum(-min.(0, runoff[tstruct.footprints[trap_ix]]))

    # there should not be negative flow out of the trap
    return max(output, 0.0), tstruct.spillpoints[trap_ix]
    
end

# ----------------------------------------------------------------------------
function _update_runoff!(runoff::Matrix{<:Real},         # modified in-place
                         region_accum::Vector{<:Real},   # modified in-place
                         offregion_runoff::Ref{<:Real},  # modified in-place
                         spoint::Spillpoint,
                         output::Real,
                         regions::Matrix{Int64},
                         dstream::Graphs.SimpleDiGraph)
    # Propagate the additional runoff from the trap downstream until it reaches the
    # bottom of the next downstream trap, or exits the domain.

    if spoint.downstream_region == 0
        # trap is spilling directly out of domain
        offregion_runoff[] += output
        return
    end

    target = spoint.downstream_region_cell
    
    # if we got here, trap is spilling to a cell inside the domain
    while target > 0 && output > 0

        infiltrated = -min(runoff[target], 0.0)
        runoff[target] += output
        output = max(output - infiltrated, 0.0)

        reg = regions[target]
        target = Graphs.outneighbors(dstream, target)
        if isempty(target)
            if (reg > 0)
                region_accum[reg] += output
            else
                offregion_runoff[] += output
            end
            target = 0
        else
            @assert length(target) == 1
            target = target[1]
        end
    end
end

# ----------------------------------------------------------------------------
function _make_simple_graph(edges::Vector{Tuple{Int, Int}})
    
    # make mapping between integers encountered in edges and 1...N
    vmap = sort(unique(vcat([x[1] for x in edges], [x[2] for x in edges])))
    vmap_inv = Dict(val => i for (i, val) in enumerate(vmap))
    
    # creating graph, using the consecutive node values 1...length(vmap)

    g = Graphs.SimpleDiGraph(
        [Graphs.SimpleEdge{Int64}(vmap_inv[e[1]], vmap_inv[e[2]]) for e in edges]
    )

    return g, vmap, vmap_inv
end

