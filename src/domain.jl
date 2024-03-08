export Domain2D, tiledomain, tiles

# ----------------------------------------------------------------------------
"""
    struct Domain2D

Type representing a Cartesian (sub)domain of a grid, defined by all
gridcells indexed (ix, iy) such that ix ∈ xrange and iy ∈ yrange.

Fields
- `xrange::UnitRange{Int}`: range of domain in the first coordinate direction
- `yrange::UnitRange{Int}`: range of domain in the second coordinate direction
"""
struct Domain2D
    xrange::UnitRange{Int}
    yrange::UnitRange{Int}
end

"""
   size(d::Domain2D)

Return size of 2D domain d
"""
function Base.size(d::Domain2D)
    return (Base.size(d.xrange)[1], Base.size(d.yrange)[1])
end

"""
    min(d::Domain2D, dir::Int)

Return minimum value of 2D domain `d`, in first or second direction
"""
function Base.min(d::Domain2D, dir::Int)
    return dir == 1 ? d.xrange[1] : d.yrange[1];
end

"""
    max(d::Domain2D, dir::Int)

Return maximum value of 2D domain `d`, in first or second direction
"""
function Base.max(d::Domain2D, dir::Int)
    return dir == 1 ? d.xrange[end] : d.yrange[end];
end

# ----------------------------------------------------------------------------
"""
    tiledomain(dom, tilenum_x, tilenum_y)

Divide a given domain up into a cartesian grid of subdomains ('tiles')

The domain 'dom' is split up into a grid of domains, with 'tilenum\\_x' number
of domain in the first coordinate direction, and 'tilenum\\_y number of domains
in the second.  

Returns a Vector{Domain2D} with all the resulting subdomains.  A second return
argument gives a tuple of integer vectors specifying where the splits happened
in the first and second coordinate directions.

See also [`Domain2D`](@ref)

"""
function tiledomain(dom::Domain2D, tilenum_x, tilenum_y)

    @assert(tilenum_x <= length(dom.xrange))
    @assert(tilenum_y <= length(dom.yrange))

    xsplit = range(dom.xrange[1], stop=dom.xrange[end]+1, length=1+tilenum_x) |>
        collect .|> floor .|> Int
    ysplit = range(dom.yrange[1], stop=dom.yrange[end]+1, length=1+tilenum_y) |>
        collect .|> floor .|> Int

    return tiles(xsplit, ysplit), (xsplit, ysplit);
end

# ----------------------------------------------------------------------------
"""
    tiles(xsplit::Vector{Int}, ysplit::Vector{Int})

Divide a domain up into a cartesian grid of subdomains ('tiles')

The vectors `xsplit` and `ysplit` should be monotonously increasing integer
vectors that specify where to split an underlying 2D region into subdomains.

The function will return a vector of [`Domain2D`](@ref) that tile the 
domain.  If `xsplit` is of length Nx and `ysplit` of length Ny, the number of tiles
generated will be (Nx-1) x (Ny-1).

"""
function tiles(xsplit::Vector{Int}, ysplit::Vector{Int})::Vector{Domain2D}

    tilenum_x = length(xsplit) - 1;
    tilenum_y = length(ysplit) - 1;

    return [Domain2D(xsplit[ix]:xsplit[ix+1]-1,
                     ysplit[iy]:ysplit[iy+1]-1)
            for ix=1:tilenum_x, iy=1:tilenum_y];
end

