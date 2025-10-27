export spillanalysis

# all outside regions represented as -1 if true
# @@ NB: Domain currently only implemented for spillfield!
"""
    spillanalysis(grid, usediags=true, building_mask=nothing, sinks=nothing,
                  lengths=nothing, domain=nothing, merge_outregions=false, 
                  verbose=false)

Analyse a terrain and compute all key information regarding its trap structure.

This information includes the spillfield, all the spill regions, traps with their 
volumes, spillpoints and footprints, the upstream/downstream trap hierarchy, 
and the supertrap/subtrap hierarchy.  

All computed information is returned as a [`TrapStructure`](@ref).  Refer to its
documentation for details.

# Arguments
- `grid::Matrix{<:Real}`: topograpical grid to analyse
- `usediags::Bool`: if true, also consider slopes along diagonals
- `building_mask::Union{Matrix{<:Bool}, BitMatrix, Nothing}`: 
      if present, provides a mask that specifies the footprint of buildings.  
      These parts of the domain will be clipped away.
- `sinks::Union{Vector{Tuple{Int, Int}}, Matrix{Bool}, Nothing}`:
      vector containing (i, j) grid coordinates of any point sinks in the grid, if any.
      Can also be a Matrix{Bool} of same size as `grid`, indicating the sink locations.
- `lengths::Union{Tuple{<:Real, <:Real}, Nothing}`: 
      tuple expressing the length and width of the grid (used to compute aspect ratios)
- `domain::Union{Domain2D, Nothing}`: 
      restrict computation to the specified domain of the grid.  @@ Note that this is not
      fully supported yet for this function.
- `merge_outregions::Bool`: if `true`, all "outside" regions will be merged and 
      represented as region -1.   Otherwise, each "outside" region will be represented
      by its own negative integer.
- `verbose::Bool`: if `true`, print information showing progress in the computation 
                   along the way.

See also [`TrapStructure`](@ref), [`fill_sequence`](@ref).
"""
function spillanalysis(grid::Matrix{<:Real};
                  usediags::Bool=true,
                  building_mask::Union{Matrix{<:Bool}, BitMatrix, Nothing}=nothing,
                  sinks::Union{Vector{Tuple{Int, Int}}, Matrix{Bool}, Nothing}=nothing,
                  lengths::Union{Tuple{<:Real, <:Real}, Nothing}=nothing,
                  domain::Union{Domain2D, Nothing}=nothing,
                  merge_outregions::Bool=false, 
                  verbose::Bool=false)
    


    verbose && println("Entering spillfield")

    if typeof(sinks) <: Matrix
        sinks = [(i[1], i[2]) for i in findall(sinks)]
    end

    field, slope = spillfield(grid, usediags=usediags,
                              lengths=lengths, domain=domain,
                              sinks=sinks,
                              building_mask=building_mask)

    verbose && println("Entering spillregions")
    regions = spillregions(field, usediags=usediags)

    verbose && println("Entering spillpoints")
    spoints, regbnd = spillpoints(grid, regions, usediags=usediags)
    
    verbose && println("Entering sshierarchy")
    subtrapgraph, lowest_regions = sshierarchy!(grid, regions, spoints, regbnd)
    
    toptraps = []
    for i in 1:Graphs.nv(subtrapgraph)
        if Graphs.outdegree(subtrapgraph, i) == 0
            push!(toptraps, i)
        end
    end

    supertraps_of = _compute_supertraps_of(lowest_regions)

    trapvols = trapvolumes(grid, regions, spoints, lowest_regions)

    footprints = _compute_trap_footprints(grid, lowest_regions, regions,
                                         spoints, maximum(regions))

    subvolumes = _compute_subvolumes(trapvols, subtrapgraph) 

    if merge_outregions
        regions[regions .< 0] .= -1
        for i = 1:length(spoints)
            if spoints[i].downstream_region <= 0
                spoints[i] = Spillpoint(-1, spoints[i].current_region_cell,
                                        spoints[i].downstream_region_cell,
                                        spoints[i].elevation)
            end
        end
    end
    
    return TrapStructure{eltype(grid)}(copy(grid),
                                       field,
                                       regions,
                                       spoints,
                                       trapvols,
                                       subvolumes, 
                                       footprints,
                                       lowest_regions,
                                       supertraps_of,
                                       subtrapgraph,
                                       building_mask,
                                       sinks)
end
    
# ----------------------------------------------------------------------------
function _compute_supertraps_of(lowest_regions)
    # produce a vector with one entry per lowest-level trap, giving the indices
    # of itself and all the supertraps that contain it

    if isempty(lowest_regions)
        return Vector{Vector{Int64}}()
    end
    num_lowlevel_regions = maximum(vcat(lowest_regions...))
    result = [Vector{Int64}() for _ in 1:num_lowlevel_regions]

    for i in 1:length(lowest_regions)
        for k in lowest_regions[i]
            push!(result[k], i)
        end
    end
    
    return result
end

# ----------------------------------------------------------------------------
function _compute_subvolumes(trapvols, subtrapgraph)
    # compute the volume of a trap that is wholly contained within its subtraps
    
    svols = zeros(length(trapvols))
    
    # from each trap, subtract the volumes of its subtraps
    for i in 1:length(svols)
        subtrap_ixs = Graphs.inneighbors(subtrapgraph, i)
        for j in subtrap_ixs
            svols[i] += trapvols[j]
        end
    end
    return svols
end

# ----------------------------------------------------------------------------
function _compute_trap_footprints(heights, lowest_subtraps_for, regions,
                                 spillpoints, num_regions)

    # the footprints vector has one entry per trap, 
    num_traps = length(spillpoints)
    footprints = [Vector{Int64}() for i in 1:num_traps]

    # Determine all traps that are influenced by each region
    supertraps_for = [Vector{Int64}() for i in 1:num_regions]
    for i in 1:num_traps
        for j in lowest_subtraps_for[i]
            push!(supertraps_for[j], i)
        end
    end

    # place highest-level traps first
    reverse!.(supertraps_for)

    for i in LinearIndices(regions)
        z = heights[i]

        (regions[i] <= 0) && continue
        
        for tr in supertraps_for[regions[i]]
            if z <= spillpoints[tr].elevation
                push!(footprints[tr], i)
            else
                # if it is not within the supertrap footprint, it will not be
                # within any remaining subtrap's footprints either
                break
            end
        end
    end
    return footprints
end
