# # Handling flat areas

# This is a simple demonstration of how large flat surfaces (like large water
# bodies) can be recognized, and straightened or excluded from analysis if
# desired.

# ## Imporing packages and loading the surface.
using SurfaceWaterIntegratedModeling
import GLMakie, Images # for visualization and loading of textures
import ColorSchemes
using Pkg.Artifacts
using DisplayAs # for image generation in documentation

# The package with SWIM testdata is provided as a Julia artifact, which can be
# accessed using the function `datapath_testdata`.  We subsequently load and
# display the synthetic grid.
datapath = joinpath(datapath_testdata(), "data", "small")
grid = loadgrid(joinpath(datapath, "bay.txt"));

## for ease of use, we create our own label of key colors found in the
## colorsheme `:Paired_12` used below
cmap = Dict(:blue => 2, :green => 4, :red => 6, :orange => 8,
            :lilac => 10, :bright => 11);

# The terrain used for demonstration is a varied surface that includes hills,
# a mountainside, a road and a significant part consisting of ocean.
tex = fill(cmap[:green], size(grid))
sf, fig, ax = plotgrid(grid, texture= tex,
                       colormap=ColorSchemes.:Paired_12,
                       colorrange=(1,12))
cam = GLMakie.cameracontrols(ax)
cam_lookat = GLMakie.Vec(285, 317, 23.9) # set observer target point
cam_eyepos = GLMakie.Vec(767, 613, 414) # set observer position
GLMakie.update_cam!(ax.scene, cam_eyepos, cam_lookat, GLMakie.Vec3f0(0, 0, 1))
DisplayAs.PNG(ax.scene)

# ## Identifying flat areas and running spill analysis
# 
# While the ocean part is clearly visible on the surface, it has not been
# clearly identified.  For this purpose, the [`identify_flat_areas`](@ref)
# function can be used:

rel_tol = 1e-3 # relative tolerance for an angle between two cells to be considered zero
max_cluster_size = 1 # threshold for excluding too small areas
isflat = identify_flat_areas(grid, rel_tol, max_cluster_size)

tex[isflat] .= cmap[:blue]
drape_surface(sf, tex)
DisplayAs.PNG(ax.scene)

# We can see that the ocean part was correctly identified as flat, but so were
# many parts that ought not be included.  We can increase the cluster size limit
# to filter these parts out.
max_cluster_size=100
isflat = identify_flat_areas(grid, rel_tol, max_cluster_size)

tex[isflat] .= cmap[:blue]
tex[.!isflat] .= cmap[:green]
drape_surface(sf, tex)
DisplayAs.PNG(ax.scene)
# With this value for the threshold, only the ocean remains identified as flat.
#
# For the surface used in this example, the surface is actually flat.  However,
# it is sometimes the case that even regions that should be conceptually flat
# contains small irregularities that may cause noisy output from the spill
# analysis.  One way to handle this, is to use the [`flatten_grid!`](@ref)
# function, which ensures that flat parts are exactly that.  Even if not 
# strictly necessary for our sample terrain, we demonstrate the use below:
flatten_grid!(grid, isflat, :min)

# Another way to exclude flat areas from the analysis is to declare them as
# "sinks".  For oceans, this makes conceptual sense, since they can be
# considered traps with "infinite" capacity.  For other large flat areas,
# e.g. parking lots, this approach does not work so well, and `flatten_grid!`
# should be used.  For our terrain, though, declaring the flat areas as sinks
# should be fine.  We pass along `isflat` as the sink map in our call to
# [`spillanalysis`](@ref) below:
tstruct = spillanalysis(grid, sinks=isflat)

# The spill anaysis was able to identify 895 spill regions and 1118 traps:
print("Number of spill regions identified: ", numregions(tstruct), '\n')
print("Number of traps identified: ", numtraps(tstruct), '\n')

# We can visualize the traps, spill regions, rivers and oceans on the
# terrain surface by creating a texture with [`show_region_selection`](@ref):
tex = show_region_selection(tstruct,
                            region_color=cmap[:green],
                            river_color=cmap[:blue]-1,
                            trap_color=cmap[:blue]);

# terrain that spills directly to the ocean, or out of the domain is
# assigned a light green color, to distinguish it from terrain that spills to
# identified traps, while the ocean is set to blue:
tex[tex.==0] .= cmap[:green]-1;
tex[isflat] .= cmap[:blue];

drape_surface(sf, tex)
DisplayAs.PNG(ax.scene)

# ## Adding more sinks
# 
# The ridge in front of the large lake in the upper left part of the plot is
# actually a railroad.  The regions immediately in front of it and behind it
# should be considered drained, so water does not massively accumulated like the
# above plot indicates.  By adding a drain near its lowest point, the lake goes
# away:
sinks = [(119,193), (180, 193)]; # two drain locations
for s in sinks
    ## add the sinks to the map that we used to represent sinks (and flat areas)
    isflat[s...] = true;
end

## run the spill analysis again
tstruct = spillanalysis(grid, sinks=isflat)

# The result of this analysis no longer indicate a large lake in this region:
tex = show_region_selection(tstruct,
                            region_color=4,
                            river_color=2,
                            trap_color=2);
tex[tex.==0] .= cmap[:green]-1
tex[isflat] .= cmap[:blue]

drape_surface(sf, tex)
DisplayAs.PNG(ax.scene)
#

# This concludes this simple demo, and we close all graphical windows
GLMakie.closeall()
