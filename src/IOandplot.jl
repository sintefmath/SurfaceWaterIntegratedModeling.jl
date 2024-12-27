import Makie
import DelimitedFiles
import Base: size, min, max
using GeometryBasics: Point3f0, Vec2f0, decompose, QuadFace, Tesselation, Rect, Mesh
using Colors, ColorSchemes

export loadgrid, savegrid, plotgrid, drape_surface, set_camerapos

# ----------------------------------------------------------------------------
"""
    loadgrid(filename, delimiter)

Load a 2D grid saved in csv format, with a specified delimiter separating 
individual entries.  Returns a Matrix{Real}.
"""
function loadgrid(filename::String; delim=',')::Matrix{Float64}
    try
        result = DelimitedFiles.readdlm(filename, delim);
    catch e
        println(e.msg)
    end
end

# ----------------------------------------------------------------------------
"""
    savegrid(filename, matrix, delimiter)

Save the 2D grid represented by `matrix` to the csv file with name `filename`, 
using a specified character as delimiter.
"""
function savegrid(filename::String, matrix; delim=',')
    try
        DelimitedFiles.writedlm(filename, matrix, delim);
    catch e
        println(e.msg)
    end
end

# ----------------------------------------------------------------------------
"""
    plotgrid(grid, texture=nothing, colormap=:lightrainbow, wireframe=false, 
             downsamplefac=1, heightfac=1, colorrange=nothing)

Plot a texture-mapped 2½D surface of a given grid, using Makie.

The function returns three values in the following order: 
- the surface object
- the figure
- the scene

# Arguments
- `grid::AbstractArray{<:Real, 2}`: 2D array representing the grid
- `texture::Union{Matrix{<:Union{Real, Colorant}}, Nothing}`: 
     optional texture for the grid.  This can be presented as a matrix of
     Colorant (directly representing colors), or of numbers (in which case the
     prevailing `colormap` and `colorrange` are used to specify the actual
     colors).  If `texture` is left to `nothing`, the height values of `grid` 
     will be directly used as a substitute.
- `colormap::Union{Symbol, ColorScheme}=:lightrainbow`: 
     specifies the colormap to use if the texture is provided as a matrix 
     of numbers.  It can either be a symbol referring to a specific colormap
     used by Makie (call `Makie.available_gradients()` for a list of available
     colormap symbols), or directly provided as a `ColorScheme`.  Note that 
     if `texture` is given as a `Matrix{Colorant}`, then `colormap` has no 
     effect.
- `wireframe::Bool`: option to display wireframe grid
- `downsamplefac::Real`: scaling factor for downsampling the grid. Useful for
                         rapid visualization of very large grids.
- `heightfac::Real`: Scaling factor for height values of the grid
- `colorrange::Tuple{<:Real, <:Real}`: specify the values representing the start
                                       and end points of `colormap`

See also [`drape_surface`](@ref)
"""
function plotgrid(grid::AbstractArray{<:Real, 2};
                  texture::Union{Matrix{<:Union{<:Real, <:Colorant}}, Nothing}=nothing,
                  colormap::Union{Symbol, ColorScheme}=:lightrainbow,
                  wireframe::Bool=false,
                  downsamplefac::Real=1.0,
                  heightfac::Real=1.0,
                  colorrange::Union{Tuple{<:Real, <:Real}, Nothing}=nothing)
    z = grid;

    if downsamplefac > 1
        downsample(range) =
            Int.(floor.(unique(ceil.(range / downsamplefac)) * downsamplefac));

        ix1 = min.(downsample(1:size(grid, 1)), size(grid, 1));
        ix2 = min.(downsample(1:size(grid, 2)), size(grid, 2));

        z = Base.view(grid, ix1, ix2);
    end

    if texture==nothing
        @show typeof(z)
        texture = z;
    end

    # define points
    res = size(z);
    x = LinRange(0, size(grid, 1), res[1]);
    y = LinRange(0, size(grid, 2), res[2]);
    points = [Point3f0(x[i], y[j], heightfac * z[i,j]) for i in 1:res[1], j in 1:res[2]];

    # define faces
    faces = decompose(QuadFace{Makie.GLIndex}, Tesselation(Rect(0, 0, 1, 1), res));

    # define parameterization
    uv = map(points) do p
        tup = ((p[1], p[2])) ./ size(grid);
        return Vec2f0(tup);
    end;

    # define normals
    normals = vec(_computenormals(x[2] - x[1], y[2] - y[1], z));

    # create mesh
    #glmesh = Mesh(Makie.meta(vec(points); uv=Makie.Buffer(vec(uv)), normals), faces);
    glmesh = Mesh(Makie.meta(vec(points); uv=vec(uv), normals), faces);

    # GLMakie handles texture coordinates differently
    textransform = _get_tex_transform(texture)
    
    if typeof(colormap) == ColorScheme
        colormap = Makie.cgrad(colormap, length(colormap), categorical=true)
    end

    # plot mesh
    if isnothing(colorrange)
        colorrange = eltype(texture) <: Real         ? extrema(texture) :
                     typeof(colormap) == ColorScheme ? (1,length(colormap)) :
                                                       (1, 256)
    end

    fig, ax, plt =
        Makie.mesh(glmesh,
                   color=textransform(texture), colormap=colormap,
                   interpolate=false, colorrange=colorrange)

    if _using_cairo()
        # get rid of transparent background
        fig.attributes[:backgroundcolor]=(:white, 1)
        ax.scene.attributes[:backgroundcolor]=(:white, 1)
    end
    
    if wireframe
        Makie.wireframe!(ax, glmesh, color=(:black, 0.1),
                           linewidth=1, transparency=true);      
    end
    return plt, fig, ax.scene
end

# ----------------------------------------------------------------------------
function _using_cairo()
    string(Makie.current_backend()) == "CairoMakie"
    #string(typeof(Makie.current_backend[])) == "CairoMakie.CairoBackend"
end

# ----------------------------------------------------------------------------
function _using_glmakie()
    string(Makie.current_backend()) == "GLMakie"
    #string(typeof(Makie.current_backend[])) == "GLMakie.GLBackend"
end

# ----------------------------------------------------------------------------
function _get_tex_transform(tex)
    if _using_cairo() && eltype(tex) <: Real
        x -> x
    else
        x -> reverse(transpose(x), dims=1)
    end
    #_using_glmakie() ? x -> reverse(transpose(x), dims=1) : x -> x
end

# ----------------------------------------------------------------------------
function _computenormals(dx, dy, z)
    diz = z[2:end,:] - z[1:end-1,:];
    djz = z[:, 2:end] - z[:, 1:end-1];
    
    # add padding
    diz = diz[[1, (1:end)..., end], :];
    djz = djz[:, [1, (1:end)..., end]];
    
    nx = (diz[1:end-1, :] + diz[2:end, :]) / (2dx);
    ny = (djz[:, 1:end-1] + djz[:, 2:end]) / (2dy);
    
    # compute normals and normalize them
    normals = Makie.normalize.([Point3f0(nx[i, j], ny[i,j], 1) for i in 1:size(z, 1), j in 1:size(z,2)]);
    return normals;
end

# ----------------------------------------------------------------------------
"""
    drape_surface(surf, tex)

Update the texture of a surface drawn with [`plotgrid`](@ref).
Here, `surf` is a mesh object (first return value of [`plotgrid`](@ref), 
whereas `tex` can be a matrix of numbers or of `Colorant` (see `plotgrid` 
documentation for specifics).

See also [`plotgrid`](@ref).
"""
function drape_surface(surf, tex)
    surf.attributes.color = _get_tex_transform(tex)(tex)
end

# ----------------------------------------------------------------------------
"""
   set_camerapos(figure, scene, cpos, ctarget, czoom)

Set the camera position, target and zoom level for a given scene.

This function is provided as a workaround to smooth over different
idiosyncracies in camera handling for the Makie backends GLMakie and 
CairoMakie.
"""
function set_camerapos(fig, scene, cpos, ctarget, czoom)
    cam = Makie.cameracontrols(scene)
    upvec = Makie.Vec3f0(0.0, 0.0, 1.0)
    # if _using_glmakie()
    #     cam.zoom_mult[] = czoom
    # elseif _using_cairo()
    cam.fov[] = 45.0 * czoom
    # else
    #     error("Unimplemented for backend:" * string(typeof(Makie.current_backend[])))
    # end        
    Makie.update_cam!(scene, cpos, ctarget, upvec)

    #scene.center[] = false
    if _using_cairo()
        return fig
    end
end
