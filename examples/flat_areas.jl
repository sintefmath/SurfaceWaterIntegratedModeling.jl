# # Handling flat areas

# This is a simple demonstration of how large flat surfaces (like large water
# bodies) can be recognized, and straightened or excluded from analysis if
# desired [^1].

# ## Imporing packages and loading the surface.
using SurfaceWaterIntegratedModeling
import GLMakie, Images # for visualization and loading of textures
import ColorSchemes
using Pkg.Artifacts

# The package with SWIM testdata is provided as a Julia artifact, which can be
# accessed using the function `datapath_testdata`.  We subsequently load and
# display the synthetic grid.
datapath = joinpath(datapath_testdata(), "data", "small")
grid = loadgrid(joinpath(datapath, "bay.txt"));

## for ease of use, we create our own label of key colors found in the
## colorsheme `:Paired_12` used below
cmap = Dict(:blue => 2, :green => 4, :red => 6, :orange => 8,
            :lilac => 10, :bright => 11);

## Also define a preset viewpoint
view1 = (GLMakie.Vec(1059, 783, 622), GLMakie.Vec(280, 311, 58), 0.68);

# The terrain used for demonstration is a varied surface that includes hills,
# a mountainside, a road, and ocean areas.
tex = fill(cmap[:green], size(grid))
sf, fig, sc = plotgrid(grid, texture = tex, 
                       colormap=ColorSchemes.:Paired_12,
                       colorrange=(1,12))
set_camerapos(sc, view1...)

# ## Identifying flat areas and running spill analysis
# 
# While the ocean part is clearly visible on the surface, its cells in the
# terrain grid have not been identified.  For this purpose, the
# [`identify_flat_areas`](@ref) function can be used:

rel_tol = 1e-3 # relative tolerance for an angle between two cells to be considered zero
max_cluster_size = 1 # threshold for excluding too small areas
isflat = identify_flat_areas(grid, rel_tol, max_cluster_size)

tex[isflat] .= cmap[:blue]
drape_surface(sf, tex)
set_camerapos(sc, view1...)


# We can see that the ocean part was correctly identified as flat, but so were
# many parts that ought not be included.  We can increase the cluster size limit
# to filter these parts out.
max_cluster_size=100
isflat = identify_flat_areas(grid, rel_tol, max_cluster_size)

tex[isflat] .= cmap[:blue]
tex[.!isflat] .= cmap[:green]
drape_surface(sf, tex)
set_camerapos(sc, view1...)

# With this value for the threshold, only the ocean remains identified as flat.
#
# For the surface used in this example, the ocean surface is already exactly
# flat.  However, it is sometimes the case with terrain input data that 
# conceptually flat regions still contains small irregularities that
# may cause noisy output from the spill analysis.  One way to handle this, is to
# use the [`flatten_grid!`](@ref) function, which ensures that flat parts are
# exactly that.  Even if not strictly necessary for our sample terrain, we
# demonstrate the use below:
flatten_grid!(grid, isflat, :min)

# Another way to exclude flat areas from the analysis is to declare them as
# "sinks".  For oceans, this makes conceptual sense, since they can be
# considered traps with "infinite" capacity.  For other large flat areas,
# e.g. parking lots, this approach does not work so well, and `flatten_grid!`
# should be used.  For our terrain, declaring the flat areas as sinks should be
# fine.  We pass along `isflat` as the sink map in our call to
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

# Terrain that spills directly to the ocean, or out of the domain is
# assigned a light green color, to distinguish it from terrain that spills to
# identified traps, while the ocean is set to blue:
tex[tex.==0] .= cmap[:green]-1;
tex[isflat] .= cmap[:blue];

drape_surface(sf, tex)
set_camerapos(sc, view1...)

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

## Run the spill analysis again
tstruct = spillanalysis(grid, sinks=isflat)

# The result of this analysis no longer indicate a large lake in this region:
tex = show_region_selection(tstruct,
                            region_color=4,
                            river_color=2,
                            trap_color=2);
tex[tex.==0] .= cmap[:green]-1
tex[isflat] .= cmap[:blue]

drape_surface(sf, tex)
set_camerapos(sc, view1...)

# ## Inserting a culvert
#
# Instead of using a sink in the depression behind the railroad, we can connect it
# to the other side by a culvert.  We specify the culvert as a connection between
# two grid cells on either side of the railroad:
for s in sinks
    isflat[s...] = false; # remove the sink previously added
end

# Specify a culvert connecting the two sides of the railroad
culvert = (CartesianIndex(119, 193), CartesianIndex(175, 193))

# Run the spill analysis again, now with the culvert specified
tstruct = spillanalysis(grid, sinks=isflat, culverts=[culvert])

# In the result, the lake behind the railroad is now drained through the culvert:
tex = show_region_selection(tstruct,
                            region_color=4,
                            river_color=2,
                            trap_color=2);
tex[tex.==0] .= cmap[:green]-1
tex[isflat] .= cmap[:blue]

tex[119, 193] = cmap[:red] # mark the culvert location in red
tex[175, 193] = cmap[:red]
tex[180, 193] = cmap[:red]

drape_surface(sf, tex)
set_camerapos(sc, view1...)
# Since there is no sink in the the depression in front of the railroad, there is water
# accumulation both in front of and behind the railroad now, with the same water level
# on both sides due to the connection provided by the culvert. 
#
# If we examine the upstream area of the trap bottom of the front depression, we see
# that this now includes the spill region behind the railroad as well:
pt_ix = LinearIndices(size(grid))[CartesianIndex(180, 193)] # trap bottom in front of railroad
upstream_cells = upstream_area(tstruct, pt_ix, local_only=false)
tex[upstream_cells] .= cmap[:orange] # mark upstream area in orange

drape_surface(sf, tex)
set_camerapos(sc, view1...)
# [^1]:
#     The data used in this example was originally obtained from
#     [Kartverket](https://kartverket.no/) (the Norwegian Mapping Authority) under the
#     [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license.
