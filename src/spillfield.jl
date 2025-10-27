import SparseArrays
import OffsetArrays: OffsetArray, center

export spillfield, update_spillfield!, vconcat_spillfields, hconcat_spillfields

# ----------------------------------------------------------------------------
# Note: explicit loops rather than matrix operations have been used several
# places in the code below.  This was to avoid matrix slicing and copying
# operations, which profiler revealed to be quite expensive when applied on very
# large grids.
# ----------------------------------------------------------------------------
"""
    spillfield(grid, usediags=true, lengths=nothing, domain=nothing, 
               tiling=nothing, building_mask=nothing)

Compute the spillfield of a raster terrain, represented by `grid`.  

The spillfield is returned an integer array of same shape as `grid`, to be
interpreted as follows: 
* -3 : sink (any passing streamline is terminated here)
* -2 : covered by building / clipped away
* -1 : no downward slope (gridcell is a trap)
*  0 : steepest slope towards (i-1, j)
*  1 : steepest slope towards (i+1, j)
*  2 : steepest slope towards (i, j-1)
*  3 : steepest slope towards (i, j+1)
*  4 : steepest slope towards (i-1, j-1)
*  5 : steepest slope towards (i+1, j+1)
*  6 : steepest slope towards (i+1, j-1)
*  7 : steepest slope towards (i-1, j+1)

Trap bottoms, i.e. cells for which there exists no downward slope, are given the
value -1.

In addition, a matrix with the value of the steepest slope in each point is
returned as a second argument.   

# Arguments
- `grid::Matrix{<:Real}`: terrain raster grid with height values
- `usediags::Bool`: if true, also consider slopes along diagonals
- `building_mask::Union{Matrix{Bool}, BitMatrix, Nothing}`: 
      a grid of logicals, specifying which cells are masked by buildings (true), 
      and thus inactive. These cells will be assigned a spill field value of -2
      (see list above).
- `sinks::Union{Vector{Tuple{Int, Int}}, Nothing}`: 
      vector containing (i, j) grid coordinates of any point sinks in the grid, if any
- `lengths::Union{Tuple{<:Real, <:Real}, Nothing}`: 
      tuple expressing the length and width of the grid (used to compute aspect ratios)
- `domain::Union{Domain2D, Nothing}`: restrict computation to the specified domain
                                      of the grid
- `tiling::Union{Tuple{Int, Int}, Nothing}`: 
      tuple specifying number of 'tiles' to subdivide surface in for parallel
      processing.  Default is (1,1), which means the whole surface is treated
      as a single tile (no parallel processing).

See also [`update_spillfield!`](@ref).
"""
function spillfield(grid::Matrix{<:Real};
                    usediags::Bool=true,
                    building_mask::Union{Matrix{<:Bool}, BitMatrix, Nothing}=nothing,
                    sinks::Union{Vector{Tuple{Int, Int}}, Nothing}=nothing,
                    lengths::Union{Tuple{<:Real, <:Real}, Nothing}=nothing,
                    domain::Union{Domain2D, Nothing}=nothing,
                    tiling::Union{Tuple{Int, Int}, Nothing}=nothing)

    if domain == nothing
        domain = Domain2D((1:x for x in size(grid))...)
    end
    
    xlen = length(domain.xrange)
    ylen = length(domain.yrange)

    dir   = Matrix{Int8}(undef, xlen, ylen)
    slope = Matrix{Float64}(undef, xlen, ylen)

    if building_mask != nothing
        # we are going to locally modify the grid, so we must copy it
        building_mask = building_mask .|> Bool  # ensure type is Bool
        grid = copy(grid)
        grid[building_mask] .= Inf
    end
    
    if tiling == nothing
        _spillfield!(dir, slope, grid,
                     usediags=usediags,
                     lengths=lengths,
                     domain=domain);

    else
        # divide domain in tiles, and use parallel processing
        tiles, = tiledomain(domain, tiling[1], tiling[2]);

        Threads.@threads for i = 1:prod(size(tiles))
            _spillfield!(dir, slope, grid,
                        usediags=usediags,
                        lengths=lengths,
                        domain=tiles[i]);
        end
    end

    # fill in any buildings
    _fill_in_buildings_and_sinks!(dir, slope, building_mask, sinks)
    
    return dir, slope
end

# ----------------------------------------------------------------------------
"""
    update_spillfield!(dir, slope, grid, domain, usediags=true, lengths=nothing)

Update an existing spill field in-place within a specific rectangular domain (where 
the topography grid has presumably changed).

# Arguments
- `dir::Matrix{Int}` : the spillfield, as described in the documentation of the
                       `spillfield` function.  Will be updated within the
                       specified domain.
- `slope::Array{<:Real}`: the steepest slope in each grid point, as returned by
                          the `spillfield` function.  Will be updated within the
                          specified domain.
- `grid::Matrix{<:Real}`: terrain raster grid with height values.  This grid has
                          presumably already been changed within the specified 
                          domain.
- `domain::Domain2D` : the domain in which to update the information in `dir` 
                       and `slope`
- `usediags::Bool=true`: if true, also consider slopes along diagonals

- `lengths::Union{Tuple{<:Real, <:Real}, Nothing}`: 
      tuple expressing the length and width of the grid (used to compute aspect ratios)
- `building_mask::Union{Matrix{Bool}, Nothing}`: 
      a grid of logicals, specifying which cells are masked by buildings (true), 
      and thus inactive. These cells will be assigned a spill field value of -2 
      (see list of possible fieldvalues in documentation of [`spillfield`](@ref).)

- `sinks::Vector{Union{Tuple{Int, Int}, Nothing}}`: 
      vector containing (i, j) grid coordinates of any point sinks in the grid, if any.

See also [`spillfield`](@ref).
"""
function update_spillfield!(dir::Matrix{Int8}, slope::Array{<:Real}, # output
                            grid::Matrix{<:Real}, domain::Domain2D;
                            building_mask=nothing, sinks=nothing,
                            usediags::Bool=true, lengths=nothing)
    # expand domain by one gridcell in all direction, since the gridcell closest
    # to the modified area may also have changed spill direction (which depends
    # on the immediate neighbors)
    domain2 = Domain2D(
        max(domain.xrange[1]-1, 1):min(domain.xrange[end]+1, size(grid, 1)),
        max(domain.yrange[1]-1, 1):min(domain.yrange[end]+1, size(grid, 2)));

    if building_mask != nothing
        # we are going to locally modify the grid, so we must copy it
        building_mask = building_mask .|> Bool  # ensure type is Bool
        grid = copy(grid)
        grid[building_mask] .= Inf
    end
    
    _spillfield!(dir, slope, grid, domain=domain2, usediags=usediags, lengths=lengths)

    # fill in any buildings
    _fill_in_buildings_and_sinks!(dir, slope, building_mask, sinks)
end

# ----------------------------------------------------------------------------
# Flag the cells in the spillfield raster grid that are covered by buildings, or 
# constitute sinks.
function _fill_in_buildings_and_sinks!(dir, slope, building_mask, sinks)

    # fill in any buildings
    if building_mask != nothing
        dir[building_mask] .= -2 # these gridcells are clipped away by buildings
        slope[building_mask] .= NaN
    end

    # fill in any sinks
    if sinks != nothing
        for pos in sinks
            dir[pos...] = -3; # flag this cell as a sink
            slope[pos...] = NaN
        end
    end
end

# ----------------------------------------------------------------------------
"""
    _spillfield!(dir, slope, grid, usediags=true, lengths=nothing)

Mutating version of spillfield.  The results 'dir' and 'slope' are not returned,
but passed as arguments.  See documentation of 'spillfield' for a general
description.  A key difference in the mutating version is that the results (dir
and slope) retain the full size of the grid, even though only a subdomain is
addressed.  (In the non-mutating version, the returned grids are limited to the
size of the addressed subdomain.)
"""
function _spillfield!(dir::Matrix{Int8}, slope::Array{<:Real}, # output args
                     grid::Matrix{<:Real}; # input arg
                     usediags::Bool=true,
                     lengths=nothing,
                     domain=nothing)
    if domain == nothing
       domain = Domain2D(1:size(grid,1), 1:size(grid, 2)); 
    end
    if lengths==nothing
        lengths = size(grid);  # default aspect ratio is 1
    end

    # physical dimensions (lengths) of each pixel (aspect ratio may matter)
    dx, dy = lengths ./ size(grid);

    # ensure the result grids have the same size as the input grid
    _setsamesize!(dir, grid);  
    _setsamesize!(slope, grid);
    
    # compute the steepest downslopes along cardinal directions
    tmpdir, tmpslope = _find_downslopes(grid, dx, dy, :axes, domain);

    if usediags
        tmpdirD, tmpslopeD = _find_downslopes(grid, dx, dy, :diags, domain);
        dix = tmpslopeD .< tmpslope; 
        tmpdir[dix] = tmpdirD[dix] .+ 4;
        tmpslope[dix] = tmpslopeD[dix]; # possible directions now in [0, 7]
    end

    tmpdir[tmpslope .>= 0] .= -1; # flag trap bottoms

    view(dir, domain.xrange, domain.yrange) .= tmpdir;
    view(slope, domain.xrange, domain.yrange) .= tmpslope; 
    
end

# ----------------------------------------------------------------------------
"""
    vconcat_spillfields(dir1, slope1, grid1, dir2, slope2, grid2, usediags=true, 
                        lengths=nothing)

Concatenate two spill fields along the 'vertical' direction (adding rows).

In addition to the spill fields `dir1` and `dir2` to concatenate, the
corresponding slopes and original terrain grids are also given (`slope1`, `slope2`, 
`grid1`, `grid2`).  These are used to re-compute the spill directions for 
gridcells on the 'seam' between the two spill fields.

The function returns the combined spill field (corresponding to `[dir1; dir2]`), as 
well as the associated slopes (corresponding to `[slope1; slope2]`)

# Arguments
- `dir1::Matrix{Int8}`: first spill field
- `slope1::Matrix{<:Real}`: matrix with local slopes for first spill field
- `grid1::Matrix{<:Real}`: topography grid from which `dir1` was computed
- `dir2::Matrix{Int8}`:  second spill field
- `slope2::Matrix{<:Real}`: matrix with local slopes for second spill field
- `grid2::Matrix{<:Real}`: topography grid from which `dir2` was computed
- `usediags::Bool`: if true, also consider slopes along diagonals
- `lengths::Union{Tuple{<:Real}, Nothing}`: 
      tuple expressing the length and width of the combined grid (used to compute
      aspect ratios)

See also [`spillfield`](@ref), [`hconcat_spillfields`](@ref).
"""
function vconcat_spillfields(dir1::Matrix{Int8}, # upper grid spillfield info
                             slope1::Matrix{<:Real},
                             grid1::Matrix{<:Real}, 
                             dir2::Matrix{Int8},
                             slope2::Matrix{<:Real},
                             grid2::Matrix{<:Real}; # lower grid spillfield info
                             usediags::Bool=true,
                             lengths::Union{Tuple{<:Real, <:Real}, Nothing}=nothing)
                                 
    # Check that grids are of compatible sizes
    @assert(size(dir1) == size(slope1) == size(grid1));
    @assert(size(dir2) == size(slope2) == size(grid2));
    @assert(size(dir1, 2) == size(dir2, 2));
    @assert(size(grid1, 1) > 1) # to enable recomputing the seam easily
    @assert(size(grid2, 1) > 1) # to enable recomputing the seam easily

    # concatenation
    dir = [dir1; dir2];
    slope = [slope1; slope2];

    # recompute the "seam"
    tmpgrid = [grid1[end-1:end, :]; grid2[1:2, :]];
    seamdomain = Domain2D(2:3, 1:size(tmpgrid,2));

    if lengths != nothing
        # readjust lengths to correspond to the seamgrid
        dx, dy = lengths ./ [size(dir1, 1) + size(dir2, 1), size(dir1, 2)]
        lengths = (dx, xy) .* size(tmpgrid)
    end

    dir_seam, slope_seam = spillfield(tmpgrid;
                                      usediags=usediags,
                                      lengths=lengths,
                                      domain=seamdomain);
        
    # overwrite the seam with the recomputed values
    istart = size(dir1, 1);
    dir[istart:istart+1, :] = dir_seam;
    slope[istart:istart+1, :] = slope_seam;

    return dir, slope
end

# ----------------------------------------------------------------------------
"""
    hconcat_spillfields(dir1, slope1, grid1, dir2, slope2, grid2, usediags=true, 
                        lengths=nothing)

Concatenate two spill fields along the 'horizontal' direction (adding columns).

In addition to the spill fields `dir1` and `dir2` to concatenate, the
corresponding slopes and original terrain grids are also given (`slope1`, `slope2`, 
`grid1`, `grid2`).  These are used to re-compute the spill directions for 
gridcells on the 'seam' between the two spill fields.

The function returns the combined spill field (corresponding to `[dir1, dir2]`), as 
well as the associated slopes (corresponding to `[slope1, slope2]`)

# Arguments
- `dir1::Matrix{Int8}`: first spill field
- `slope1::Matrix{<:Real}`: matrix with local slopes for first spill field
- `grid1::Matrix{<:Real}`: topography grid from which `dir1` was computed
- `dir2::Matrix{Int8}`:  second spill field
- `slope2::Matrix{<:Real}`: matrix with local slopes for second spill field
- `grid2::Matrix{<:Real}`: topography grid from which `dir2` was computed
- `usediags::Bool`: if true, also consider slopes along diagonals
- `lengths::Union{Tuple{<:Real}, Nothing}`: 
      tuple expressing the length and width of the combined grid (used to compute
      aspect ratios)

See also [`spillfield`](@ref), [`vconcat_spillfields`](@ref).
"""
function hconcat_spillfields(dir1::Matrix{Int8}, # upper grid spillfield info
                             slope1::Matrix{<:Real},
                             grid1::Matrix{<:Real}, 
                             dir2::Matrix{Int8},
                             slope2::Matrix{<:Real},
                             grid2::Matrix{<:Real}; # lower grid spillfield info
                             usediags::Bool=true,
                             lengths::Union{Tuple{<:Real, <:Real}, Nothing}=nothing)
                                 
    # Check that grids are of compatible sizes
    @assert(size(dir1) == size(slope1) == size(grid1));
    @assert(size(dir2) == size(slope2) == size(grid2));
    @assert(size(dir1, 1) == size(dir2, 1));
    @assert(size(grid1, 2) > 1) # to enable recomputing the seam easily
    @assert(size(grid2, 2) > 1) # to enable recomputing the seam easily

    # horizontal concatenation
    dir = [dir1 dir2];
    slope = [slope1 slope2];

    # recompute the "seam"
    tmpgrid = [grid1[:, end-1:end]  grid2[:, 1:2]];
    
    seamdomain = Domain2D(1:size(tmpgrid,1), 2:3);

    if lengths != nothing
        # readjust lengths to correspond to the seamgrid
        dx, dy = lengths ./ [size(dir1, 1), size(dir1, 2) + size(dir2, 2)]
        lengths = (dx, xy) .* size(tmpgrid)
    end
    
    dir_seam, slope_seam = spillfield(tmpgrid;
                                      usediags=usediags,
                                      lengths=lengths,
                                      domain=seamdomain);
        
    # overwrite the seam with the recomputed values
    jstart = size(dir1, 2);
    dir[:, jstart:jstart+1] = dir_seam;
    slope[:, jstart:jstart+1] = slope_seam;

    return dir, slope
end

# ----------------------------------------------------------------------------
# Identify and compute the steepest downwards slope along the four cardinal, or
# the four diagonal directions.
function _find_downslopes(grid, dx, dy, direction, domain)


    @assert(direction == :axes || direction == :diags)
    diag = direction == :diags;
    dxy = sqrt(dx^2 + dy^2);

    delta1, delta2 = diag ? (dxy, dxy)        : (dx, dy);
    step1, step2   = diag ? ((1, 1), (-1, 1)) : ((1, 0), (0, 1));
 
    d1grid = _diffgrid(grid, delta1, step1, domain);
    d2grid = _diffgrid(grid, delta2, step2, domain);
    
    orient1, slopes1 = _compare_slopes(d1grid, d1grid, -1.0, step1);
    orient2, slopes2 = _compare_slopes(d2grid, d2grid, -1.0, step2);

    orient, slope = _compare_slopes(slopes1, slopes2, 1.0, (0, 0));

    dir = similar(orient, Int8);
    dir[.!orient] = orient1[.!orient];
    dir[orient] = orient2[orient] .+ 2; # possible directions now in [0, 3]

    return dir, slope;
    
end

# ----------------------------------------------------------------------------
# For two grids g1 and g2 of equal size, containing local slopes (as computed by
# differentiating a raster grid in one given direction), identify the smallest
# local slope of the two for each (i, j) coordinate.  Return which one was smallest,
# as well as its value.
function  _compare_slopes(g1, g2, fac1, offset2)

    minslope = Array{eltype(g1)}(undef, size(g1) .- abs.(offset2));

    shape = size(minslope);  
    choice = Array{Bool}(undef, shape...);

    # identify lower-left corner in original grids
    corner1 = map(x->firstindex(g1, x), (1, 2)) .+ max.(offset2, 0) .- 1;
    corner2 = map(x->firstindex(g2, x), (1, 2)) .+ max.(offset2, 0) .- 1;
    
    for col = 1:shape[2]
        for row = 1:shape[1]
            lval = g1[row + corner1[1] - offset2[1],
                      col + corner1[2] - offset2[2]] * fac1;
            rval = g2[row + corner2[1], col + corner2[2]];

            pick = rval < lval;
            choice[row, col] = pick;
            minslope[row, col] = pick ? rval : lval;
        end
    end
    
    return choice, minslope
end

# ----------------------------------------------------------------------------
# Return an OffsetArray d such that the forward/backward finite difference
# derivative at grid[i,j] are given by d[i-shift[1],j-shift[2]] and d[i,j].  If
# a given dimension 'k' in 'grid' is indexed (a:b), then the corresponding
# dimension in the resulting derivative grid will be (c:d), where
# c = min(a, a-shift[k]) and d = max(b, b - shift[k]).
function _diffgrid(grid, delta, shift, dom)

    domsize = size(dom);
    gsize = size(grid);
    deltainv = 1.0/delta;

    # range of the resulting derivative grid
    resultrangeI = _rangeunion(dom.xrange, dom.xrange .- shift[1]);
    resultrangeJ = _rangeunion(dom.yrange, dom.yrange .- shift[2]);

    # required grid points (ignoring boundary issues)
    requiredI = (dom.xrange[1]-abs(shift[1])) : (dom.xrange[end] + abs(shift[1]));
    requiredJ = (dom.yrange[1]-abs(shift[2])) : (dom.yrange[end] + abs(shift[2]));

    # the part of the required grid points that are actually available
    # when taking boundary into account
    availableI = max(requiredI[1], 1) : min(requiredI[end], gsize[1]);
    availableJ = max(requiredJ[1], 1) : min(requiredJ[end], gsize[2]);
    
    # establishing target storage grid
    target = Array{eltype(grid)}(undef, domsize .+ abs.(shift));
    target[:] .= NaN

    # The range of result values that are computable given the available grid points
    internalI =
        (max(availableI[1], availableI[1] - shift[1]) : 
        min(availableI[end], availableI[end] - shift[1]))
    internalJ =
        (max(availableJ[1], availableJ[1] - shift[2]) : 
        min(availableJ[end], availableJ[end] - shift[2]))

    # Establish target grid with proper indexing
    result = OffsetArray(target, resultrangeI[1]-1, resultrangeJ[1]-1);

    # compute finite differences where possible (e.g. in the interior)
    for j in internalJ
        for i in internalI
            result[i, j] = (grid[i+shift[1], j+shift[2]] - grid[i, j]) * deltainv;
        end
    end
    
    # extrapolate where necessary
    sgn = prod(shift);    
    if resultrangeI[1] < internalI[1]
        Jtarget = (resultrangeJ[1] - min(sgn, 0), resultrangeJ[end] - max(sgn, 0));
        Jsource = Jtarget .+ sgn;
        result[resultrangeI[1], range(Jtarget...)] =
            result[resultrangeI[1] + 1, range(Jsource...)];
    end
    if resultrangeI[end] > internalI[end]
        Jtarget = (resultrangeJ[1] + max(sgn, 0), resultrangeJ[end] + min(sgn, 0))
        Jsource = Jtarget .- sgn;
        result[resultrangeI[end], range(Jtarget...)] =
            result[resultrangeI[end] - 1, range(Jsource...)];
    end
    if resultrangeJ[1] < internalJ[1]
        Itarget = (resultrangeI[1] - min(sgn, 0), resultrangeI[end] - max(sgn, 0));
        Isource = Itarget .+ sgn;
        result[range(Itarget...), resultrangeJ[1]] =
            result[range(Isource...), resultrangeJ[1] + 1];
    end
    if resultrangeJ[end] > internalJ[end]

        Itarget = (resultrangeI[1] + max(sgn, 0), resultrangeI[end] + min(sgn, 0))
        Isource = Itarget .- sgn;
        result[range(Itarget...), resultrangeJ[end]] =
            result[range(Isource...), resultrangeJ[end] - 1];
    end
    
    return result;
end

# ----------------------------------------------------------------------------
function _rangeunion(r1, r2)
    return min(r1[1], r2[1]) : max(r1[end], r2[end]);
end

# ----------------------------------------------------------------------------
function _rangeintersect(r1, r2)
    return max(r1[1], r2[1]) : min(r1[end], r2[end]);
end

# ----------------------------------------------------------------------------
# Ensure target grid has same size as modelgrid.  If it is already the case,
# leave it alone.
function _setsamesize!(targetgrid, modelgrid)
    if size(targetgrid) != size(modelgrid)
       targetgrid = Matrix{eltype(targetgrid)}(undef, size(modelgrid)...)
    end
end        

