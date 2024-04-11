import Graphs
export Spillpoint, spillpoints

# --------------------------- Spillpoint structure ---------------------------
"""
    struct Spillpoint

A struct representing the spillpoint of a trap.  It has the following fields:
- `downstream_region::Int`: index of the downstream region.
- `current_region_cell::Int`: index of the cell in the current region bordering on 
                              the spillpoint
- `downstream_region_cell::Int`: index of the cell in the downstream region 
                                 bordering on the spillpoint
- `elevation::Real`: the elevation (vertical height) of the spillpoint

"""
struct Spillpoint
    downstream_region::Int
    current_region_cell::Int
    downstream_region_cell::Int
    elevation::Real

end

""" 
Default constructor of a Spillpoint instance sets all indices to zero
and elevation to infinity.
"""
Spillpoint() = Spillpoint(0, 0, 0, Inf)
    


"""
    spillpoints(grid, spillregions, usediags)

Compute the spillpoint positions for each (low-level) spill region in the grid.

Also computes the topology of the lowest-level spill network (the graph
describing the upstream/downstream relationship between spill regions).

Returns a vector of [`Spillpoint`](@ref) with one entry per spill region
(numbered from 1 upwards.  Each [`Spillpoint`](@ref) contains information on
downstream spill region index, cell indices for the cell in the current region
and downstream region that border on the spill point, and the heigh volume of
the spill point.

In addition, returns a vector of 'boundaries', one per spill point.  The boundary for 
spill region 'i' is defined by all the cell pairs (i, j) such that cell 'i', and 'j' 
are neighbors, cell 'i' is in spill region 'i', and cell 'j' is in a different spill
region 'j'.

# Arguments
- `grid::Matrix{<:Real}` - the terrain grid (matrix) with height values
- `spillregions::Matrix{Int}`: the matrix containing the regions already computed 
                                from the [`spillregions`](@ref) function
- `usediags::Bool` : if `true` (default), diagonal connections between cells will also be
                     considered
- `tiling::Union{Tuple{Int, Int}, Nothing}`: 
        tuple specifying number of 'tiles' to subdivide surface in for parallel
        processing.  Default is (1,1), which means the whole surface is treated
        as a single tile (no parallel processing).

See also [`Spillpoint`](@ref).
"""
function spillpoints(grid::Matrix{<:Real}, spillregions::Matrix{Int};
                     usediags::Bool=true, tiling=nothing)

    domain = Domain2D(1:size(grid,1), 1:size(grid,2))
    
    if tiling == nothing
        spillpoints, boundaries =
            _process_domain(grid, spillregions, usediags, domain)
    else
        tiles, = tiledomain(domain, tiling...)
        spoints_vec = Vector{Vector{Spillpoint}}(undef, length(tiles))
        boundaries_vec = Vector{Vector{Vector{Tuple{Int, Int}}}}(undef, length(tiles))

        Threads.@threads for i = 1:length(tiles)
            spoints_vec[i], boundaries_vec[i] =
                _process_domain(grid, spillregions, usediags, tiles[i])
        end

        spillpoints = spoints_vec[1]
        boundaries = boundaries_vec[1]

        # patch together boundaries from the different processed parts of the
        # domain, and determine which spillpoints to keep.
        for i = 2:length(tiles)
            if return_boundaries
                append!(boundaries, boundaries_vec[i])
            end
            spillpoints = map((s1, s2) -> s1.elevation < s2.elevation ? s1 : s2,
                              spillpoints, spoints_vec[i])
        end
    end

    # add domain boundaries to boundaries if user wants them returned
    _add_outer_bounderies!(boundaries, spillpoints, grid, spillregions)
    
    return spillpoints, boundaries
end

# ----------------------------------------------------------------------------
function _process_domain(grid::Matrix{<:Real}, spillregions::Matrix{Int},
                         usediags::Bool, domain::Domain2D)
    # each tuple contains (second region, cell in this region, cell in second
    # region, z-value).

    # Note: we do not care about the number of distinct spillregions, only the
    # number of distinct _positive_ spillregions, since negative spill regions
    # are spilling out of the domain anyway.
    result = fill(Spillpoint(), max(maximum(spillregions), 0))
    
    LI = LinearIndices(spillregions)
    IFirst = CartesianIndex(domain.xrange[1], domain.yrange[1])
    ILast = CartesianIndex(min(domain.xrange[end]+1, size(grid, 1)),
                           min(domain.yrange[end]+1, size(grid, 2)))
    boundaries = [Vector{Tuple{Int, Int}}() for i in 1:maximum(spillregions)]

    function _update_spillpointlist!(pos, shift)
        
        reg1 = spillregions[pos]
        reg2 = spillregions[pos+shift]

        if reg1 == 0 || reg2 == 0
            # encountered a masked value (inactive part of grid)
            return
        end
        
        if reg1 != reg2
            lpos1 = LI[pos]
            lpos2 = LI[pos+shift]

            # register this point as part of the respective region boundaries
            reg1 > 0 ? push!(boundaries[reg1], (lpos1, lpos2)) : nothing
            reg2 > 0 ? push!(boundaries[reg2], (lpos2, lpos1)) : nothing

            # check if the point is a spillpoint
            zval = max(grid[pos], grid[pos+shift])
            if reg1 > 0 && zval < result[reg1].elevation
                result[reg1] = Spillpoint(reg2, lpos1, lpos2, zval)
            end
            if reg2 > 0 && zval < result[reg2].elevation
                result[reg2] = Spillpoint(reg1, lpos2, lpos1, zval)
            end
        end
    end
    
    # look for spillpoints along main axes and (optionally) diagonals
    shifts = [CartesianIndex(1, 0), CartesianIndex(0, 1)]
    if usediags
        append!(shifts, [CartesianIndex(1,1), CartesianIndex(1, -1)])
    end
    for shift in shifts

        startix = IFirst - min(shift, CartesianIndex(0, 0))
        endix   = ILast - max(shift, CartesianIndex(0, 0))
                               
        for ix in startix:endix
            _update_spillpointlist!(ix, shift)
        end
    end

    return result, boundaries
end

# ----------------------------------------------------------------------------
function _add_outer_bounderies!(boundaries, spillpoints, grid, spillregions)
    # 'boundaries' and 'spillpoints' will be updated

    LI = LinearIndices(spillregions)
    NI, NJ = size(spillregions)

    function _update(region, cell)
        if region <= 0
            # We are not in a region, so there is no boundary to update
            return
        end
        push!(boundaries[region], (cell, cell))
        zval = grid[cell]
        if zval < spillpoints[region].elevation
            # Repeating the index means that the second index (which should be in
            # the neighbour region) is out of the grid boundary, i.e. out of the domain
            spillpoints[region] = Spillpoint(0, cell, cell, zval)
        end
    end
        
    for i = 1:NI
        _update(spillregions[i, 1],  LI[CartesianIndex(i,  1)])
        _update(spillregions[i, NJ], LI[CartesianIndex(i, NJ)])
    end
    
    for j = 1:NJ
        _update(spillregions[1, j], LI[CartesianIndex(1, j)])
        _update(spillregions[NI, j], LI[CartesianIndex(NI, j)])
    end
end
