import Graphs
export spillanalysis
import Images

# all outside regions represented as -1 if true
# @@ NB: Domain currently only implemented for spillfield!
"""
    spillanalysis(grid; usediags=true, building_mask=nothing, sinks=Vector{CartesianIndex{2}}(),
                  lengths=nothing, domain=nothing, merge_outregions=false, verbose=false,
                  culverts=Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}}(),
                  barriers=Vector{Vector{CartesianIndex{2}}}())

Analyse a terrain and compute all key information regarding its trap structure.

This information includes the spillfield, all the spill regions, traps with their 
volumes, spillpoints and footprints, the upstream/downstream trap hierarchy, 
and the supertrap/subtrap hierarchy.  

All computed information is returned as a [`TrapStructure`](@ref).  Refer to its
documentation for details.

# Arguments
- `grid::Matrix{<:Real}`: topograpical grid to analyse
- `usediags::Bool=true`: if true, also consider slopes along diagonals
- `building_mask::Union{Matrix{<:Bool}, BitMatrix, Nothing}=nothing`: 
      if present, provides a mask that specifies the footprint of buildings.  
      These parts of the domain will be clipped away.
- `sinks::Union{Vector{CartesianIndex{2}}, Matrix{Bool}}=Vector{CartesianIndex{2}}()`:
      vector containing (i, j) grid coordinates of any point sinks in the grid, if any.
      Can also be a Matrix{Bool} of same size as `grid`, indicating the sink locations.
- `lengths::Union{Tuple{<:Real}, Nothing}=nothing`: 
      tuple expressing the length and width of the grid (used to compute aspect ratios)
- `domain::Union{Domain2D, Nothing}=nothing`: 
      restrict computation to the specified domain of the grid.  @@ Note that this is not
      fully supported yet for this function.
- `merge_outregions::Bool=false`: if `true`, all "outside" regions will be merged and 
      represented as region -1.   Otherwise, each "outside" region will be represented
      by its own negative integer.
- `verbose::Bool=false`: if `true`, print information showing progress in the computation 
                   along the way.
- `culverts::Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}}=Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}}()`:
      vector of culverts, each defined by a pair of grid coordinates. Culverts allow flow
      between two cells that would otherwise be blocked by terrain.
- `barriers::Vector{Vector{CartesianIndex{2}}}=Vector{Vector{CartesianIndex{2}}}()`:
      vector of barriers, where each barrier is a polyline defined by a sequence of grid 
      coordinates. Barriers block flow between cells along the polyline.

See also [`TrapStructure`](@ref), [`fill_sequence`](@ref).
"""
function spillanalysis(grid::Matrix{<:Real};
                  usediags::Bool=true,
                  building_mask::Union{Matrix{<:Bool}, BitMatrix, Nothing}=nothing,
                  sinks::Union{Vector{CartesianIndex{2}}, Matrix{Bool}}=Vector{CartesianIndex{2}}(),
                  waterbodies::Union{Matrix{<:Bool}, Nothing}=nothing,
                  lengths::Union{Tuple{<:Real}, Nothing}=nothing,
                  domain::Union{Domain2D, Nothing}=nothing,
                  merge_outregions::Bool=false, 
                  verbose::Bool=false,
                  culverts::Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}}=
                           Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}}(),
                  barriers::Vector{Vector{CartesianIndex{2}}}=
                           Vector{Vector{CartesianIndex{2}}}())

    verbose && println("Entering spillfield")

    # a copy of the grid will be returned by the TrapStructure, and there may be
    # modifications to it if waterbodies are provided
    gridcpy = copy(grid) 
    
    # ensure `sinks` is a vector of CartesianIndex
    if typeof(sinks) <: Matrix
        sinks = [CartesianIndex(i[1], i[2]) for i in findall(sinks)]
    elseif sinks == nothing
        sinks = Vector{CartesianIndex{2}}()
    end

    # if waterbodies are provided, ensure that the regions they cover are of uniform height
    if waterbodies !== nothing
        _flatten_waterbody_regions!(gridcpy, waterbodies)
    end
    
    # ensure culverts are directed from higher to lower elevation
    directed_culverts = [ gridcpy[x[1]] >= gridcpy[x[2]] ? x : (x[2], x[1]) for x in culverts ]
    
    # determine any connections cut by barriers
    cut_edges = Set{Tuple{CartesianIndex{2}, CartesianIndex{2}}}()
    for b in barriers
        @assert length(b) >= 2 "Each barrier must have at least two points"
        e = polyline_grid_intersections([p.I for p in b], usediags=usediags)
        union!(cut_edges, e)
    end
    cut_edge_dict = edgeset2dict(cut_edges)
    
    field, slope = spillfield(gridcpy, usediags=usediags,
                              lengths=lengths, domain=domain,
                              sinks=sinks,
                              blocked_edges=cut_edge_dict,
                              building_mask=building_mask)

    verbose && println("Entering spillregions")
    regions, flowgraph, trap_bottoms =
        spillregions(field, usediags=usediags,
                     grid=gridcpy, # for eliminating 'flat' traps
                     cut_edges=cut_edge_dict, culverts=directed_culverts)
    
    verbose && println("Entering spillpoints")
    spoints, regbnd = spillpoints(gridcpy, regions, usediags=usediags,
                                  cut_edges=cut_edge_dict)

    # verbose && println("Eliminating flat traps")
    # _eliminate_flat_traps!(regions, spoints, flowgraph, trap_bottoms, gridcpy)
    
    verbose && println("entering sshierarchy")
    subtrapgraph, lowest_regions = sshierarchy!(gridcpy, regions, spoints, regbnd)
    
    toptraps = []
    for i in 1:Graphs.nv(subtrapgraph)
        if Graphs.outdegree(subtrapgraph, i) == 0
            push!(toptraps, i)
        end
    end

    supertraps_of = _compute_supertraps_of(lowest_regions)

    trapvols = trapvolumes(gridcpy, regions, spoints, lowest_regions)

    footprints = _compute_trap_footprints(gridcpy, lowest_regions, regions,
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

    wbody_cells = isnothing(waterbodies) ? Vector{CartesianIndex{2}}() : findall(waterbodies)
    
    return TrapStructure{eltype(grid)}(gridcpy,
                                       flowgraph,
                                       trap_bottoms,
                                       regions,
                                       spoints,
                                       trapvols,
                                       subvolumes, 
                                       footprints,
                                       lowest_regions,
                                       supertraps_of,
                                       subtrapgraph,
                                       building_mask,
                                       sinks,
                                       wbody_cells,
                                       cut_edge_dict)
end

# ----------------------------------------------------------------------------
function _flatten_waterbody_regions!(grid, waterbodies)
    # for each contiguous waterbody region, set the elevation of all cells in that
    # region to the same value (the minimum elevation within the region)
    wbody_labels = Images.label_components(waterbodies)
    for label in 1:maximum(wbody_labels)
        region_cells = findall(wbody_labels .== label)
        min_elev = minimum(grid[region_cells])
        grid[region_cells] .= min_elev
    end
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

# # ----------------------------------------------------------------------------
# # Helper function to eliminate flat traps (traps with spill points at the same
# # elevation as their bottoms). This is not necessary for the spillpoint
# # algorithms to work, but reduces unnecessary complexity in the trap structure
# # and can speed up subsequent computations.
# # In this function, the following arguments are mutable:
# # - regions
# # - spillpoints
# # - flowgraph
# # - trap_bottoms
# # Note that the `grid` argument is not mutable, but is needed to determine the
# # elevation of the spill points and trap bottoms.
# function _eliminate_flat_traps!(regions, spillpoints, flowgraph, trap_bottoms, grid)

#     lowest_elevations = fill(-Inf, length(spillpoints))
#     for tb in trap_bottoms
#         reg = regions[tb]
#         if reg > 0
#             lowest_elevations[reg] = grid[tb]
#         end
#     end
#     @assert all(lowest_elevations .> -Inf) "All traps should be covered."
#     @assert length(spillpoints) == length(lowest_elevations) "One spill point per trap."
    
#     traps_to_eliminate = []
#     # loop through regions and identify those that are 'flat' and should be eliminated
#     for reg in 1:length(spillpoints)
#         if isapprox(lowest_elevations[reg], spillpoints[reg].elevation; atol=1e-6)
#              push!(traps_to_eliminate, reg)
#         end
#     end
#     println("Number of traps to eliminate: ", length(traps_to_eliminate)) # @@
    
#     # loop through connections in flowgraph, and reconnect all bottom cells in
#     # flat traps directly to the corresponding spillpoint cell (ensure that existing
#     # connections are eliminated)
#     linear = LinearIndices(size(grid))
#     tb_to_remove = []
#     for tb_ix in 1:length(trap_bottoms)
#         tb = trap_bottoms[tb_ix]
#         tb_linear = linear[tb]
#         reg = regions[tb]
#         if reg > 0 && reg ∈ traps_to_eliminate
#             sp_cell = spillpoints[reg].current_region_cell
#             target = (tb_linear == sp_cell) ? spillpoints[reg].downstream_region_cell : sp_cell

#             # delete all outgoing edges from tb
#             outneighs = Graphs.outneighbors(flowgraph, tb_linear)
#             for neighbor in outneighs
#                 Graphs.rem_edge!(flowgraph, tb_linear, neighbor)
#             end
#             # add the new, unique edge
#             # @@ If we want contiguous spill paths, we can replace the below line with
#             # an algorithm determining the correct path
#             Graphs.add_edge!(flowgraph, tb_linear, target)

#             push!(tb_to_remove, tb_ix)
#         end
#     end
#     deleteat!(trap_bottoms, tb_to_remove)
        
#     # Loop through traps to eliminate, set the associated spill region equal to
#     # the downstream region, update flowgraph, remove the spillpoint, and
#     # renumerate all traps/regions.
#     regcells = regioncells(regions) # make lookups quicker
    
#     for reg in traps_to_eliminate
#         downstream_reg = spillpoints[reg].downstream_region

#         # assign cells from region 'reg' to downstream region
#         append!(regcells[downstream_reg], regcells[reg])
#         empty!(regcells[reg])

#     end

#     # delete entries correpsonding to eliminated cells, and regenerate (positive
#     # entries of) regions matrix
#     deleteat!(regcells, traps_to_eliminate)
#     for r_ix = 1:length(regcells)
#         regions[regcells[r_ix]] .= r_ix
#     end

#     # update spillpoint list (including recomputing the downstream regions of
#     # spillpoints, as they might have changed
#     deleteat!(spillpoints, traps_to_eliminate)
#     for sp_ix = 1:length(spillpoints)
#         target_region = regions[spillpoints[sp_ix].downstream_region_cell]
#         if spillpoints[sp_ix].downstream_region != target_region
#             spillpoints[sp_ix] = Spillpoint(target_region,
#                                             spillpoints[sp_ix].current_region_cell,
#                                             spillpoints[sp_ix].downstream_region_cell,
#                                             spillpoints[sp_ix].elevation)
#         end
#     end
# end
