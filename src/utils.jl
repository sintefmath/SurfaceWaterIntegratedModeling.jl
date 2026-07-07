import Graphs
import Roots
#using Infiltrator @@@ for debugging

export flatten_grid!, raise_buildings!, identify_flat_areas, toplevel_traps,
    show_region_selection, all_subtraps_of, interpolate_timeseries,
    trap_states_at_timepoints, compute_spillfield_graph, all_upstream_regions,
    upstream_area, current_upstream_area,
    polyline_grid_intersections, edgeset2dict, reconstruct_spillfield,
    flatten_small_traps, waterheight

# ----------------------------------------------------------------------------
"""
    interpolate_timeseries(tstruct, seq, timepoints; 
                           filled_color=1, trap_color=2, river_color=3)

Compute the exact terrain fill states for a sequence of timepoints, given a
trapping structure and sequence of `SpillEvent`s.

Each generated fill state is represented as an integer matrix
(`Matrix{Int}`), where submerged cells are set to `filled_color`, the
parts of traps that are not (yet) submerged set to `trap_color` and intermittent
streams set to `river_color`.  Other cells are attributed the value zero.

The result is returned as a `Vector{Matrix{Int}}`, of the same length as
`timepoints`.  A second return value provides a corresponding `Vector{Int}` that
for each timepoint provides the index for the latest preceding event in `seq`.

# Arguments
- `tstruct::TrapStructure{<:Real}`: trap structure object describing the terrain traps
- `seq::Vector{SpillEvent}`: the sequence of events, as computed by
                              [`fill_sequence`](@ref) for a given weather scenario
- `timepoints::Vector{<:Real}`: the timepoints for which we want to compute the exact
                                terrain fill states.  Should be given in ascending
                                order.
- `filled_color::Int`: The 'color' value to represent filled areas (default: 1).
- `trap_color::Int`: The 'color' to represent unfilled parts of traps (default: 2).
- `river_color::Int`: The 'color' to represent intermittent rivers (default: 3).
- `verbose::Bool`: Whether to print progess information during computations (default: true).

See also [`spillanalysis`](@ref), [`fill_sequence`](@ref), 
[`trap_states_at_timepoints`](@ref).
"""
function interpolate_timeseries(tstruct::TrapStructure{<:Real},
                                seq::Vector{SpillEvent},
                                timepoints::Vector{<:Real};
                                filled_color::Int=1,
                                trap_color::Int=2,
                                river_color::Int=3,
                                verbose::Bool=true)

    (issorted(timepoints) && (seq[1].timestamp <= timepoints[1])) || 
        error("Timepoint sequence should be strictly increasing and within bounds.")

    tstates = trap_states_at_timepoints(tstruct, seq, timepoints, verbose=verbose)

    result = Vector{Matrix{Int}}()
    tix = Vector{Int}()
    for i in 1:length(tstates)
        verbose && println("Generating timepoint: ", i)
        filled_traps = tstates[i][1]
        push!(result, _fill_state_to_terrainmap(tstruct, filled_traps, tstates[i][2],
                                                filled_color, trap_color, river_color))
        push!(tix, tstates[i][3])
    end
    return result, tix
end

# ----------------------------------------------------------------------------
"""
    trap_states_at_timepoints(tstruct, seq, timepoints)

Compute the exact amount of water in all traps for a specific set of timepoints,
given a trapping structure and sequence of `SpillEvent`s.

The result is returned as a vector of three-element tuples, with one entry per
timepoint.  
- The first element of the tuple is a `Vector{Bool}` with one entry
  per trap, indicating if the trap is filled or not at the specified timepoint.
- The second element is a `Vector{Real}` with one entry per trap, giving the 
  amount of water in that trap at the specified timepoint.
- The third element is an index into `seq`, pointing to the last `SpillEvent` to
  occur before the timepoint.

# Arguments
- `tstruct::TrapStructure{<:Real}`: trap structure object describing the terrain traps
- `seq::Vector{SpillEvent}`: the sequence of events, as computed by
                              [`fill_sequence`](@ref) for a given weather scenario
- `timepoints::Vector{<:Real}`: the timepoints for which we want to compute the exact
                                water amounts of the traps.  Should be given in 
                                ascending order.
- `verbose::Bool`: Whether to print progess information during computations (default: true).

See also [`spillanalysis`](@ref), [`fill_sequence`](@ref), 
[`interpolate_timeseries`](@ref).
"""
function trap_states_at_timepoints(tstruct::TrapStructure{<:Real},
                                   seq::Vector{SpillEvent},
                                   timepoints::Vector{<:Real};
                                   verbose::Bool=true)

    num_traps, num_regions = numtraps(tstruct), numregions(tstruct)
    
    # various checks on input arguments
    @assert length(seq[1].filled) == num_traps
    @assert length(seq[1].inflow) == num_traps
    @assert length(seq[1].amount) == num_traps

    (issorted(timepoints) && (seq[1].timestamp <= timepoints[1])) || 
        error("Illegal timepoint sequence entered by user." *
              "It should be strictly increasing and within bounds.")
    
    Smin, Smax = zeros(num_traps), zeros(num_traps)
    z_vol_tables = SurfaceWaterIntegratedModeling._compute_z_vol_tables(tstruct)
    result = []
    ix = 1
    
    for tp in timepoints
        verbose && println("Handling timepoint: ", tp)
        # find the sequence point of departure
        while ix < length(seq) && seq[ix+1].timestamp <= tp
            ix += 1
        end
        @assert seq[ix].timestamp <= tp
        filled = filled_at(seq, ix)
        transitory_traps = findall(.!filled) # if they are filled at start of
                                             # interval [i, i+1], they will
                                             # remain so at least until the end
                                             # of it compute exact fill level
                                             # for all traps that are not full
        cur_runoff = runoff_at(seq, ix)
        cur_inflow = inflow_at(seq, ix)
        rateinfo = RateInfo(cur_runoff, Smax, Smin, cur_inflow)
        cur_amount = amount_at(seq, ix)
        
        for trap ∈ transitory_traps
            # we only need Smin and Smax for the target trap and its children.
            # Leave the rest uncomputed
            _compute_Smin_Smax_for_specific_trap!(rateinfo, trap, tstruct)

            cur_amount[trap] = FilledAmount(
                SurfaceWaterIntegratedModeling._compute_exact_fill(
                    rateinfo, cur_amount, trap, filled,
                    tstruct, tp, z_vol_tables, false), tp)
        end
        push!(result, (filled, [ca.amount for ca ∈ cur_amount], ix))
    end
    return result
end

# ----------------------------------------------------------------------------
# helper function for trap_states_at_timepoints

function _compute_Smin_Smax_for_specific_trap!(rateinfo, trap, tstruct)
    children = subtrapsof(tstruct, trap)
    for c in vcat(children, trap)
        # we need these in order to compute Smin for the requested trap.  `_ponding_infiltration`
        # excludes the cells at or above each trap's spillpoint (they never pond as part of
        # that trap), matching `_update_Smin_Smax!` and the dynamic network solver.
        setsmax!(rateinfo, c, _ponding_infiltration(rateinfo, tstruct, c))
    end
    setsmin!(rateinfo, trap, sum(getsmax(rateinfo, children)))
end

# ----------------------------------------------------------------------------
function _ensure_subtraps_filled!(step_filled, tstruct, last_step)

    function fill_subtraps_recur(ix)
        subs = Graphs.inneighbors(tstruct.agglomerations, ix)
        for s in subs
            if step_filled[s] == 0 # the sub has not yet been filled
                step_filled[s] = last_step
                fill_subtraps_recur(s)
            end
        end
    end
    
    just_filled = findall(step_filled .== last_step)
    
    for ix in just_filled
        fill_subtraps_recur(ix)
    end
end
    
# ----------------------------------------------------------------------------
function _find_trapcells(tstruct, trap)
    regs = tstruct.lowest_subtraps_for[trap]
    zspill = tstruct.spillpoints[trap].elevation
    trapcells = Vector{Vector{Int64}}()
    for r in regs
        regcells = findall(tstruct.regions[:] .== r)
        regtrapcells = tstruct.topography[regcells] .< zspill
        push!(trapcells, regcells[regtrapcells])
    end
    
    return vcat(trapcells...)
end

# ----------------------------------------------------------------------------
function _regionmap(tstruct)
    # we only map positive regions (those that do not spill out of domain)
    result = [Vector{Int64}() for i = 1:maximum(tstruct.regions)]
    
    for i in LinearIndices(tstruct.regions)
        val = tstruct.regions[i]
        if val > 0
            push!(result[val], i)
        end
    end
    return result
end

# ----------------------------------------------------------------------------
function _fill_state_to_terrainmap(tstruct::TrapStructure{<:Real},
                                   filled::Vector{Bool}, 
                                   trapwater::Vector{<:Real},
                                   filled_color::Int, trap_color::Int, river_color::Int)

    # if a trap is filled, show it along with the river running out of it.
    # Ignore any "trap" with zero volume
    selection = findall(filled .&& tstruct.trapvolumes .> 0)
    result = show_region_selection(tstruct, selection=selection,
                                   region_color=0, trap_color=filled_color,
                                   river_color=river_color)

    # indicate toplevel traps if requested (whether or not filled)
    if trap_color != 0
        selection = findall(tstruct.trapvolumes .> 0)
        toplevel_traps =
            show_region_selection(tstruct, selection=selection, region_color=0,
                                  trap_color=trap_color, river_color=trap_color)
        ixs = (result .== 0)
        result[ixs] = toplevel_traps[ixs]
    end
    # # if a trap is not filled, but all subtraps are filled, we need to plot its
    # # exact fill trap_is_filled
    for trap in findall(.!filled)
        if trapwater[trap] <= sqrt(eps())
            # there is no water in the trap yet - nothing to plot
            continue
        end
        subtraps = Graphs.inneighbors(tstruct.agglomerations, trap)
        if all(filled[subtraps])
            trapcells = tstruct.footprints[trap]
            z_cur = waterheight(tstruct, trap, trapwater[trap])
            
            covered_trapcells = trapcells[tstruct.topography[trapcells] .< z_cur]
            result[covered_trapcells] .= filled_color
        end
    end
    
    return result
end

# # ----------------------------------------------------------------------------
# This earlier implemenation of `waterheight` runs into trouble if z_spill=Inf.
# Keep in here for now, as it might be more robust in some cases?
# function waterheight(tstruct, trap, current_watervolume; tol=sqrt(eps()))
#     ownvolume = tstruct.trapvolumes[trap] - tstruct.subvolumes[trap]
#     remaining = ownvolume - current_watervolume
    
#     trapcells = tstruct.footprints[trap]
#     z_spill = tstruct.spillpoints[trap].elevation
#     bottom = minimum(tstruct.topography[trapcells])
#     bracket = (bottom, z_spill)
#     fac = 1 - tol # avoid problem with bracketing interval due to roundoff error
#     volfun = z -> sum(z_spill .- max.(z, tstruct.topography[trapcells]))

#     z_cur = Roots.find_zero(z -> volfun(z) - remaining*fac, bracket)
#     return z_cur
# end

# ----------------------------------------------------------------------------
"""
        waterheight(tstruct, trap, current_watervolume)
Compute the water height in a trap, given the current water volume.

The water volume provided should not exceed the total volume of the trap, minus
the volume of its subtraps. The function will compute the water height by
finding the level at which the volume of water in the trap (net of any subtraps)
matches the provided `current_watervolume`.
"""
function waterheight(tstruct, trap, current_watervolume; tol=sqrt(eps()))
    total_watervolume = current_watervolume + tstruct.subvolumes[trap]
    trapcells = tstruct.footprints[trap]
    z_spill = tstruct.spillpoints[trap].elevation
    bottom = minimum(tstruct.topography[trapcells])
    bracket = (bottom, z_spill)
    fac = 1 - tol # avoid problem with bracketing interval due to roundoff error
    volfun = z -> sum(max.(z .- tstruct.topography[trapcells], 0.0))

    z_cur = Roots.find_zero(z -> volfun(z) - total_watervolume*fac, bracket)
    return z_cur
end


# ----------------------------------------------------------------------------
"""
    flatten_grid!(grid::Matrix{<:Real}, mask::Matrix{<:Bool}, height_choice::Symbol)

Flatten indicated areas of a terrain grid.

The grid `grid` represents height values of a terrain, and will be modified by
the function.  The `mask` is a similarly sized, boolean grid that identifies
which areas should be flattened (all cells with value `true`).  Each such
connected region is then assigned a fixed height value and `grid` is modified
accordingly.  

There are three ways of computing the height, indicated by the `height_choice`
argument.  Valid options are:

- `:min` - the height of each flat area becomes the minimum value of its grid cells
- `:max` - the height of each flat area becomes the  maximum value of its grid cells
- `:mean` - the height of each flat area becomes the mean value of its grid cells

See also: [`identify_flat_areas`](@ref)

"""
function flatten_grid!(grid::Matrix{<:Real},  # will be modified
                       mask::Matrix{<:Bool},  # constant
                       height_choice::Symbol) # constant
    cl = _identify_clusters(mask)

    hchoice = Dict([(:min, minimum),
                    (:max, maximum),
                    (:mean, x -> sum(x)/length(x))])

    hfun = hchoice[height_choice]
    
    for c in cl
        grid[c] .= hfun(grid[c])
    end
    
end

# ----------------------------------------------------------------------------
"""
    raise_buildings!(grid, building_mask, elevation_above_roof=10.0)

Raise terrain cells covered by buildings to represent flat-roofed buildings in
the elevation model.

For each contiguous building footprint, all cells within the footprint are set
to the same elevation: the maximum terrain elevation found within that footprint
plus `elevation_above_roof`.  This models each building as having a flat roof
at a uniform height across its entire base.

Modifies `grid` in place.  This is a preprocessing step for including buildings
in the terrain analysis rather than clipping them away.  After calling this
function the terrain can be passed directly to [`spillanalysis`](@ref) without
a `clip_mask`, and water will be routed around and over the raised building
cells.

# Arguments
- `grid::Matrix{<:Real}`: terrain raster grid with height values.  Modified in place.
- `building_mask::Matrix{<:Bool}`: boolean grid of the same shape as `grid`, where
      `true` marks cells covered by a building footprint.
- `elevation_above_roof::Real=5.0`: height (in the same units as `grid`) added
      above the highest terrain point within each building footprint, giving the
      uniform roof elevation for that building.

See also [`spillanalysis`](@ref), [`flatten_grid!`](@ref).
"""
function raise_buildings!(grid::Matrix{<:Real},
                          building_mask::Matrix{<:Bool},
                          elevation_above_roof::Real=5.0)
    clusters = _identify_clusters(building_mask)
    for c in clusters
        grid[c] .= maximum(grid[c]) + elevation_above_roof
    end
end

# ----------------------------------------------------------------------------
function _identify_clusters(mask::Matrix{<:Bool})

    edges = Vector{Tuple{Int, Int}}()
    Dx = CartesianIndex(1, 0)
    Dy = CartesianIndex(0, 1)
    lind = LinearIndices(mask)
    count = 0;
    for I in CartesianIndex(1,1):CartesianIndex(size(mask) .- 1)
        if mask[I]
            count = count + 1;
            mask[I + Dx] && push!(edges, (lind[I], lind[I + Dx]))
            mask[I + Dy] && push!(edges, (lind[I], lind[I + Dy]))
        end
    end
    g = Graphs.SimpleGraph([Graphs.SimpleEdge{Int64}(e) for e in edges])
    clusters = Graphs.connected_components(g)
    # filter out clusters with only one cell, since otherwise, all cells
    # belong to a cluster, regardless of the mask.
    clusters = filter(c -> length(c) > 1, clusters)
    
    return clusters
end



# ----------------------------------------------------------------------------
"""
    identify_flat_areas(grid, rel_tol, min_cluster_size, lengths=nothing)

Identify areas in the grid that are flat within a specified tolerance.  Can be
used to detect lakes in a terrain.

Returns a boolean grid of same shape as the input grid, where cells that are 
considered to belong to flat areas are flagged as 'true' (rest are false).

# Arguments
- `grid::Matrix{<:Real}`: terrain raster grid with height values
- `rel_tol::Real`: tolerance to use in determine when a slope between two grid
                   cells may be considered 'flat'.  The tolerance is specified
                   relatively to the maximum slope present in the grid.
- `min_cluster_size::Int`: Specify how large an agglomeration of "flat" cells
                           needs to be (in terms of number of cells) in order 
                           to be included in the result.  This is used to 
                           filter out small fragments that are not considered 
                           important.
- `lengths::Tuple{Int, Int}`: Specify the length of the grid in x and y 
                             directions.  This is needed to correctly handle
                             the aspect ratio when computing slopes.  If no 
                             argument provided, the resolution of the grid 
                             is used as a substitute.

See also: [`flatten_grid!`](@ref)
"""
function identify_flat_areas(grid::Matrix{<:Real}, rel_tol::Real,
                             min_cluster_size::Int,
                             lengths::Union{Tuple{Int, Int}, Nothing}=nothing)

    result = fill(false, size(grid))
        
    if lengths==nothing
        lengths=size(grid)
    end
    slopes_x = abs.(grid[2:end, :] - grid[1:end-1, :])./lengths[1]
    slopes_y = abs.(grid[:, 2:end] - grid[:, 1:end-1])./lengths[2]

    # identify maximum and minimum slopes, and threshold for "flat" regions
    minsl, maxsl = extrema((extrema(slopes_x)..., extrema(slopes_y)...))
    threshold = maxsl * rel_tol
    if minsl > threshold
        # nothing to do.  No slope is below thresholdd
        return result
    end

    # Building graph connecting neighboring gridcells where slope is below
    # threshold
    (Nx, Ny) = size(grid)
    linind = LinearIndices(grid)
    edges = Vector{Tuple{Int, Int}}()
    for iy = 1:Ny
        for ix = 1:Nx
            if (ix != Nx && slopes_x[ix, iy] < threshold)
                push!(edges, (linind[ix, iy], linind[ix+1, iy]))
            end
            if (iy != Ny && slopes_y[ix, iy] < threshold)
                push!(edges, (linind[ix, iy], linind[ix, iy+1]))
            end
        end
    end
    
    g = Graphs.SimpleGraph([Graphs.SimpleEdge{Int64}(e) for e in edges])

    components = Graphs.connected_components(g)

    # sufficiently large components are identified as flat areas
    for c in components
        if length(c) > min_cluster_size
            for ix in c
                result[ix] = true
            end
        end
    end
        
    return result
    
end

# ----------------------------------------------------------------------------
"""
    toplevel_traps(subtrap_graph)

Identify all top-level traps in the given subtrap graph.

Top-level traps are those that those that are not subtraps of yet larger lakes.
The result is returned as a vector of indices to those traps.

# Arguments
- `subtrap_graph::Graphs.SimpleDiGraph`: oriented graph expressing the trap/subtrap 
                                         graph structure.  If A is a subtrap of B, there 
                                         will be an edge pointing from A to B.  As such, 
                                         the toplevel traps are those with not outwards-
                                         pointing edges.
See also: [`sshierarchy!`](@ref)
"""
function toplevel_traps(subtrapgraph::Graphs.SimpleDiGraph)
    toptraps = []
    for i in 1:Graphs.nv(subtrapgraph)
        if Graphs.outdegree(subtrapgraph, i) == 0
            push!(toptraps, i)
        end
    end
    return toptraps
end

# ----------------------------------------------------------------------------
"""
    all_subtraps_of(subtrap_graph, trap_ixs)

Return all subtraps (at any level) of a specified set of traps.

The graph representing the supertrap/subtrap tree is given by `subtrap_graph`. 
The indices to the set of traps is specified by `trap_ixs`.  This argument
may be an integer refering to a single trap, or a `Vector{Int}` referring to 
one or more traps.

The result is returned as a `Vector{Int}` giving the indices of all traps that 
are subtraps of the trap(s) referred to by `trap_ixs`.
"""
function all_subtraps_of(subtrap_graph::Graphs.SimpleDiGraph,
                         trap_ixs::Union{Int, Vector{Int}})

    result = Vector{Int64}()
    active_set = Set(trap_ixs)
    while !isempty(active_set)
        found = Vector{Int64}()
        for i in active_set
            union!(found, Graphs.inneighbors(subtrap_graph, i))
        end
        union!(result, found)
        active_set = found
    end
    return result
end


# ----------------------------------------------------------------------------
"""
    show_region_selection(tstruct; selection=nothing, 
                          region_color=1, trap_color=2, river_color=3)

Create a texture identifying the requested regions, traps and related rivers.

The information about the terrain, including topography, spillpoints, etc., is
provided by the structure `tstruct`, which can be obtained from the
[`spillanalysis`](@ref) function.  The selected regions are identified by
`selection`, a vector with the indices of the selected regions.  If left empty,
all toplevel traps/regions will be selected.  

The result is returned as a matrix of integers, where the integers can be
thought of as 'colors'. The integer values used to indicate region, trap or
river are given by the arguments `region_color`, `trap_color` and `river_color`,
respectively.  If `region_color` is a negative integer, each region will be
assigned an unique color, starting from `abs(region_color)` and incrementing by
one per region.

# Arguments
- `selection`: integer vector containing indices of the selected regions
- `region_color`: integer to assign to the cells of all selected regions
- `trap_color`: integer to assign to all selected traps 
- `river_color`: integer to assign to all river cells related to the selected
                 regions

See also: [`spillanalysis`](@ref)
"""
function show_region_selection(tstruct::TrapStructure{<:Real};
                               selection::Union{Nothing, Vector{Int}}=nothing, 
                               region_color::Int=1, trap_color::Int=2,
                               river_color::Int=3)::Matrix{Int}

    subtrapgraph = tstruct.agglomerations
    regions      = tstruct.regions
    grid         = tstruct.topography
    spillpoints  = tstruct.spillpoints
    #spillfield   = tstruct.spillfield
    flowgraph    = tstruct.flowgraph
    regmap       = tstruct.lowest_subtraps_for

    if selection == nothing
        # by default, show top-level traps
        selection = toplevel_traps(subtrapgraph)
    end

    result = zeros(Int64, size(regions)...)

    # return early if no internal regions
    if maximum(regions) < 1
        return result
    end

    # trace rivers
    wbodies = Set(tstruct.waterbodies) # convert to set for faster lookup
    if river_color != 0
        CI = CartesianIndices(size(grid))
        for reg in selection
            if spillpoints[reg].downstream_region_cell > 0
                #rix = CI[spillpoints[reg].downstream_region_cell]; # river startpoint
                rix = spillpoints[reg].downstream_region_cell; # river startpoint, linear ix
                finished = false
                while !finished
                    result[rix] = river_color
                    cell_old = CI[rix]
                    rix, finished = _downstream_cell(flowgraph, rix)
                    # in case we are crossing a 'flat' area, the downstream
                    # cell might not be directly connected to the previous
                    # cell.  In this case, trace the line connecting the two cells.
                    cell_new = CI[rix]
                    if cell_new ∈ wbodies
                        finished = true # always finished when reaching a water body
                    elseif !_are_connected(cell_old, cell_new)
                        line_ixs = _connect_cells(cell_old, cell_new)
                        result[line_ixs] .= river_color
                    end
                end
            end
        end
    end

    supertrap_lookup = fill(Int(0), maximum(regions))
    for s in selection
        supertrap_lookup[regmap[s]] .= s
    end

    selected_region_lookup = fill(Int(0), maximum(regions))

    if region_color >= 0
        # one common color for all selected regions
        all_selected_subregions = vcat(regmap[selection]...)
        selected_region_lookup[all_selected_subregions] .= region_color
    else
        # one unique color per region
        cur_color = abs(region_color)
        for reg in selection
            subregs = regmap[reg]
            selected_region_lookup[subregs] .= cur_color
            cur_color += 1
        end
    end

    for r in CartesianIndices(regions)

        if regions[r] > 0
            strap_ix = supertrap_lookup[regions[r]]
            if strap_ix > 0 && grid[r] <= spillpoints[strap_ix].elevation
                result[r] = trap_color
            else
                if result[r] == 0 # no river drawn on this cell
                    result[r] = selected_region_lookup[regions[r]]
                end
            end
        end
        
    end
    
    return result
end

# ----------------------------------------------------------------------------
@inline function _connect_cells(cell1::CartesianIndex, cell2::CartesianIndex)
    # return the linear indices of the cells along the line connecting cell1 and cell2
    # using Bresenham's line algorithm
    x1, y1 = Tuple(cell1)
    x2, y2 = Tuple(cell2)
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = x1 < x2 ? 1 : -1
    sy = y1 < y2 ? 1 : -1
    err = dx - dy

    line_ixs = Vector{CartesianIndex{2}}()
    
    while true
        push!(line_ixs, cell1)
        if cell1 == cell2
            break
        end
        err2 = err * 2
        if err2 > -dy
            err -= dy
            cell1 += CartesianIndex(sx, 0)
        end
        if err2 < dx
            err += dx
            cell1 += CartesianIndex(0, sy)
        end
    end
    
    return line_ixs
end

# ----------------------------------------------------------------------------
@inline function _are_connected(cell1::CartesianIndex, cell2::CartesianIndex)
    # check if the cells are neighbors, along any of the eight cardinal directions
    return all(abs.(Tuple(cell1) .- Tuple(cell2)) .<= 1)
end

# ----------------------------------------------------------------------------
@inline function _downstream_cell(flowgraph::Graphs.SimpleDiGraph, lix::Int)
    outneighs = Graphs.outneighbors(flowgraph, lix)
    if isempty(outneighs)
        return lix, true
    else
        return outneighs[1], false
    end
end

# ----------------------------------------------------------------------------
"""
        show_downstream_path(tstruct, startpoint; trap_color=2, river_color=3)

Show the downstream path from a given point in the terrain grid, as traced on a raster.

The path is traced by following the flow directions from the given starting
point, as given by the flow graph in `tstruct`.  The path is shown on a raster
grid, where the cells along the path are colored with `river_color`.  If the
path passes through any traps, those are colored with `trap_color`.
The result is returned as a matrix of integers, where the integers can be
thought of as 'colors'.  Cells along the downstream path are colored with
`river_color`, while cells belonging to traps along the path are colored with
`trap_color`.  Other cells are attributed the value zero.

# Arguments
- `tstruct::TrapStructure{<:Real}`: trap structure object describing the terrain
- `startpoint::Int`: linear index to the terrain grid cell from which the downstream path should be traced
- `trap_color::Int`: integer to represent trap cells along the path (default: 2)
- `river_color::Int`: integer to represent river cells along the path (default: 3)
"""
function show_downstream_path(tstruct, startpoint; trap_color=2, river_color=3)

    paths, traps = flow_path_from(tstruct, startpoint)

    result = zeros(Int64, size(tstruct.topography)...)

    for p in paths
        for cix in p
            result[cix] = river_color
        end
    end
    for t in traps
        footprint = tstruct.footprints[t]
        result[trapcells] .= trap_color
    end
    return result
end


# # ----------------------------------------------------------------------------
# function current_upstream_area(tstruct::TrapStructure{<:Real}, point, tstates; recur=1)
#     if recur >= 10000 @@@ stack overflow guard
#         return []
#     end
#     submerged, puddlecells = _is_submerged(tstruct, point, tstates)
#     if !submerged
#         # if we are not submerged, the only upstream area is the one consisting of
#         # cells spilling directly into the current cell
#         result = findall(Graphs.dfs_parents(tstruct.flowgraph, point, dir=:in) .> 0)
#     else
#         # all cells that immediately spills into the puddle, i.e. the whole local region
#         puddle_regions = Set(tstruct.regions[puddlecells])
#         result = collect(findall(in(puddle_regions), tstruct.regions[:]))
#     end
#     # add contribution of upstream spill regions that actively spill into the current
#     utraps = _in_dynamic_spillpath(tstruct, result, tstates)
#     for ut in utraps
#         @show ut
#         spillpoint = tstruct.spillpoints[ut].current_region_cell
#         append!(result, current_upstream_area(tstruct, spillpoint, tstates, recur=recur+1))
#     end
#     return unique(result)
# end

function _trap_stack_of(tstruct, pt)

    # check if the cell is in a region at all.  If not, it cannot be in a trap,
    # and we can return nothing.
    region = tstruct.regions[pt]
    if region <= 0
        return []
    end

    region_supertraps = tstruct.supertraps_of[region]

    # return stack of traps whose spill points are above the current point
    pt_z = tstruct.topography[pt]
    keep = pt_z .<= [x.elevation for x in tstruct.spillpoints[region_supertraps]]

    return region_supertraps[keep]
end

function _is_submerged(tstruct, pt, tstates)

    ltrap = _trap_stack_of(tstruct, pt)

    if isempty(ltrap)
        # the cell is not within a trap, so it cannot be submerged
        return false, nothing
    end    
    # we are in a trap, but is the cell submerged?  We can check this by looking
    # at the current fill state of the trap, and comparing it to the current
    # water level in it.
    filled_traps = tstates[1]
    trap_water_content = tstates[2]

    # find the uppermost trap in the stack that contains water
    ix = findlast(trap_water_content[ltrap] .> 0)
    if ix == nothing
        # no trap in the stack contains water, so we are not submerged.
        return false, nothing
    end

    target_trap = ltrap[ix]

    if filled_traps[target_trap]
        # we are submerged in this trap, and no higher trap contains water, so
        # puddle is limited to the footprint of this trap
        return true, tstruct.footprints[target_trap]
    end
    
    # If we got here, we are in a trap that is partially filled with water.  We
    # need to determine whether the waterlevel submerges the point in question.
    z_cur = waterheight(tstruct, target_trap, trap_water_content[target_trap])

    if tstruct.topography[pt] < z_cur
        # the cell is submerged, and the 'puddle' consists of all cells in the trap
        # with elevation below z_cur
        footprint = tstruct.footprints[target_trap]
        puddle_cells = footprint[tstruct.topography[footprint] .< z_cur]
        return true, puddle_cells
    else
        # the cell is not submerged, so we are done
        return false, nothing
    end
end

# function _in_dynamic_spillpath(tstruct, cells, tstates)
    
#     # Check if there are any filled traps that are spilling into one of the
#     # cells, while not part of the cells themselves.

#     # @ Note: the recurring loop through spillpoints may be an efficiency
#     # bottleneck if the recursion runs deeply, and there are many spillpoints.
#     # If this turns out to be a problem, we should optimize by precomputing
#     # trap spillgraph and pass it along.

#     # To avoid having to search through all cells each time (they may be numerous)
#     # we first identify the regions they represent, so we can do a quick screening
#     # first
#     regions = unique(tstruct.regions[cells])

#     # Now we can check for each filled trap if it is spilling into any of the regions
#     filled_traps = tstates[1]
#     spillpoints = tstruct.spillpoints
#     utraps = []
#     for i in 1:length(spillpoints)

#         if (filled_traps[i] &&
#             spillpoints[i].downstream_region > 0 &&
#             spillpoints[i].downstream_region ∈ regions &&
#             spillpoints[i].downstream_region_cell ∈ cells &&
#             spillpoints[i].current_region_cell ∉ cells)
#             push!(utraps, i)
#         end
#     end
#     return utraps
# end
# ----------------------------------------------------------------------------
"""
    upstream_area(tstruct, point, local_only=true)

Determine all grid cells belonging to the upstream area of a given point location
in the terrain grid.

Result is returned as a `Vector{Int}`, giving the linear indices of the terrain 
grid cells that constitute the upstream area of the specified point.

# Arguments
- `tstruct::TrapStructure{<:Real}`: trap structure object describing the terrain traps
- `point::Int`: linear index to the terrain grid cell for which the upstream area
                is requested
- `local_only::Bool`: if true, consider only the immediate spill region.  If 
                      otherwise, consider also upstream spill regions.

See also [`all_upstream_regions`](@ref).
"""
function upstream_area(tstruct::TrapStructure{<:Real},
                       point::Int;
                       local_only::Bool=true)
    # identify region
    region = tstruct.regions[point]
    
    # If the point is within a trap, the whole spill region is considered upstream
    in_trap = (region > 0) &&
        (tstruct.topography[point] <=
        tstruct.spillpoints[tstruct.supertraps_of[region][end]].elevation)
    # convert point to cartesian coordinate
    cpt = CartesianIndices(tstruct.topography)[point]
    is_sink = (cpt ∈ tstruct.sinks)
    if in_trap
        if local_only
            return findall(tstruct.regions[:] .== region)
        else
            # topmost subtrap of the current region
            toptrap = maximum(tstruct.supertraps_of[region])

            # identify regions belonging to the upstream area of the topmost trap
            num_regions = length(tstruct.supertraps_of)
            lookup = fill(false, num_regions) # used to keep trap of upstream regions
            
            # loop over all subtraps of the topmost trap
            sgraph = saturated_spillgraph(tstruct) # cache this for speedup
            for s in tstruct.lowest_subtraps_for[toptrap]
                lookup[all_upstream_regions(tstruct, s; sgraph=sgraph)] .= true
            end
            
            # identify all grid cells belonging to any of the identified upstream regions
            tmp = fill(false, size(tstruct.topography))
            for i=1:length(tmp[:])
                r = tstruct.regions[i]
                tmp[i] = (r > 0) ? lookup[r] : false
            end
            return findall(tmp[:])
        end
    elseif is_sink
        # @@@ ideally this should be restructured together with the outside-trap
        # case, as there is some redundancy.
        if local_only
            return findall(tstruct.regions[:] .== region)
        else
            utraps = filter(x -> tstruct.spillpoints[x].downstream_region == region,
                            1:length(tstruct.spillpoints))
            ucells = findall(tstruct.regions[:] .== region)
            for ut in utraps
                append!(ucells,
                        upstream_area(tstruct,
                                      tstruct.spillpoints[ut].current_region_cell,
                                      local_only=false))
            end
        end
        return unique(ucells)
    else
        
        # outside a trap, the only upstream area is the one consisting of cells
        # spilling directly into the current cell

        # upstream cells in region
        ucells = findall(Graphs.dfs_parents(tstruct.flowgraph, point, dir=:in) .> 0) 
        
        if local_only
            return ucells
        else
            # check if we are directly in the path of a spill path from an upstream trap
            utraps = filter(x -> tstruct.spillpoints[x].downstream_region == region,
                            1:length(tstruct.spillpoints))
            
            ucells_loc = copy(ucells)
            for ut in utraps
                if tstruct.spillpoints[ut].downstream_region_cell ∈ ucells_loc
                    append!(ucells,
                            upstream_area(tstruct,
                                          tstruct.spillpoints[ut].current_region_cell,
                                          local_only=false))
                end
            end
        end
        return unique(ucells)
    end
end

# ----------------------------------------------------------------------------
"""
    all_upstream_regions(tstruct, region)

Identify all spill regions that will eventually become part of the extended
watershed of a specified region.

When traps are empty, water accumulates locally.  As traps fill up, water starts
to spill from upstream to downstream regions.  This function identifies all the
spill regions that will flow into a specified region/trap when all upstream
traps have been completely filled.

The result is returned as a `Vector{Int}`, giving the indices of all regions 
that will be upstream of the specified region when all traps are filled.

# Arguments
- `tstruct::TrapStructure{<:Real}`: trap structure object describing the terrain traps
- `region::Int`: index of the region in question

# Optional argument
- `sgraph`::Tuple{Graphs.SimpleDiGraph, Vector{Int}, Dict{Int, Int}}: the spill graph,
  as computed by `saturated_spillgraph`.  Providing it can speed up repeated calls to this
  function.
See also [`upstream_area`](@ref)
"""
function all_upstream_regions(tstruct::TrapStructure{<:Real},
                              region::Int;
                              sgraph::Union{Nothing,
                                            Tuple{Graphs.SimpleDiGraph{Int},
                                                  Vector{Int},
                                                  Dict{Int, Int}}}=nothing
                              )

    g, vmap, vmap_inv = (sgraph == nothing) ? saturated_spillgraph(tstruct) : sgraph
    
    reg_graphnode_ix = vmap_inv[tstruct.supertraps_of[region][end]]
    upstream_graph_nodes = findall(Graphs.dfs_parents(g, reg_graphnode_ix, dir=:in) .> 0)
    upstream_regions = vmap[upstream_graph_nodes]

    # only consider the lowest-level regions
    num_regions = length(tstruct.supertraps_of)
    filter!(x -> x <= num_regions, upstream_regions)

    return upstream_regions
end

# ----------------------------------------------------------------------------
function _compute_direction(ix1::CartesianIndex{2}, ix2::CartesianIndex{2}, undef_val)
    dix = ix2 - ix1

    dir = dix == CartesianIndex(-1, 0)  ? 0 :
          dix == CartesianIndex(1, 0)   ? 1 :
          dix == CartesianIndex(0, -1)  ? 2 :
          dix == CartesianIndex(0, 1)   ? 3 :
          dix == CartesianIndex(-1, -1) ? 4 :
          dix == CartesianIndex(1, 1)   ? 5 :
          dix == CartesianIndex(1, -1)  ? 6 :
          dix == CartesianIndex(-1, 1)  ? 7 :
          undef_val
end
# ----------------------------------------------------------------------------
function reconstruct_spillfield(tstruct::TrapStructure{<:Real})
    xmax, ymax = size(tstruct.topography)
    undef_val = -4
    spillfield = fill(Int8(undef_val), xmax, ymax)

    # fill in regular flow directions (0..7) and trap cells (-1)
    LI = LinearIndices(size(spillfield))
    CI = CartesianIndices(size(spillfield))

    for cix in CartesianIndices(size(spillfield))
        outneigh = Graphs.outneighbors(tstruct.flowgraph, LI[cix])
        if length(outneigh) == 0
            spillfield[cix] = -1
        else
            @assert length(outneigh) == 1 # for now, only single downstream cell supported
            spillfield[cix] = _compute_direction(cix, CI[outneigh[1]], undef_val)
        end
        
    end
        
    # fill in sinks (-3)
    for sink in tstruct.sinks
        spillfield[sink] = -3
    end

    # fill in clipped-away areas (-2)
    if tstruct.clip_mask != nothing
        for I in CartesianIndices(tstruct.clip_mask)
            if tstruct.clip_mask[I]
                spillfield[I] = -2
            end
        end
    end

    return spillfield
end

# ----------------------------------------------------------------------------
"""
    compute_spillfield_graph(spillfield::Matrix{Int8})

Compute a Graphs.SimpleDiGraph representation from a spillfield matrix.

The `spillfield` variable is a matrix that can be generated by the 
[`spillfield`](@ref) function.
"""
function compute_spillfield_graph(spillfield::Matrix{Int8})

    xmax, ymax = size(spillfield)
    spilldomain = view(spillfield, 1:xmax, 1:ymax)

    # As sizehint, we use a value slightly larger than the number of cells in
    # the domain
    edges = sizehint!(Vector{Tuple{Int, Int}}(),
        spilldomain |> size |> prod |> (x)->x*1.1 |> ceil |> Int)

    # Identify grid edges constituting "streamlines"
    SurfaceWaterIntegratedModeling._spillfield_flow_edges!(edges, spilldomain)

    edge_filt = filter(a -> a[1] != a[2], edges)

    # Create directed graph of connections
    g  = Graphs.SimpleDiGraph([Graphs.SimpleEdge{Int64}((e[1], e[2])) for e in edge_filt])

    # Ensure all cells are present as vertices, including isolated trap bottoms
    # surrounded by masked cells that produce no flow edges.
    Graphs.add_vertices!(g, xmax * ymax - Graphs.nv(g))

    return g
end

# ----------------------------------------------------------------------------
"""
Take a polyline and subdivide lines into segments starting and ending precisely
on grid nodes, thus creating an equivalent polyline with potentially more segments.
No segments on the resulting polyline should run through grid nodes except at start
and end points.
"""
function _subdivide_polyline(pl::Vector{CartesianIndex{2}})
    @assert length(pl) >= 2
    result = [pl[1]]
    for i in 1:length(pl)-1
        Dx = pl[i+1][1] - pl[i][1]
        Dy = pl[i+1][2] - pl[i][2]
        if Dx == 0 && Dy == 0
            continue
        end
        subdiv = gcd(abs(Dx), abs(Dy))
        Dx_step = Dx ÷ subdiv
        Dy_step = Dy ÷ subdiv
        for j in 1:subdiv
            newpt = CartesianIndex(pl[i][1] + j*Dx_step, pl[i][2] + j*Dy_step)
            push!(result, newpt)
        end
    end
    return result
end

# ----------------------------------------------------------------------------
function polyline_grid_intersections(pl::Vector{CartesianIndex{2}}; usediags=true)
    @assert length(pl) >= 2
    result = Set{Tuple{CartesianIndex{2}, CartesianIndex{2}}}()
    subpl = _subdivide_polyline(pl)

    # intesections excluding segment endpoints
    for i in 1:length(subpl)-1
        seg_edges = _grid_segment_interior_intersections(subpl[i], subpl[i+1], usediags)
        union!(result, seg_edges)
    end

    # intersections at segment endpoints, excluding start of first and end of last
    isects = Set{Tuple{CartesianIndex{2}, CartesianIndex{2}}}()
    for i in 2:length(subpl)-1
        isects = _intersections_at_node(subpl[i-1], subpl[i], subpl[i+1], isects, usediags)
        union!(result, isects)
    end
    
    return result
end

# ----------------------------------------------------------------------------
function _intersections_at_node(prevpt, curpt, nextpt, prev_isects, usediags)
    @assert prevpt != curpt
    @assert nextpt != curpt
    v1 = (prevpt[1] - curpt[1], prevpt[2] - curpt[2])
    v2 = (nextpt[1] - curpt[1], nextpt[2] - curpt[2])

    result = Set{Tuple{CartesianIndex{2}, CartesianIndex{2}}}()
    
    flags = fill(false, usediags ? 8 : 4)
    ix_start = _dir2subquadrant(v1, usediags) # from 1 to 8 or 4
    ix_end   = _dir2subquadrant(v2, usediags)
    if ix_start == ix_end
        # neighbor segments enter and exit the same subquadrant - no
        # new intersections to process
        return result
    elseif ix_end < ix_start
        ix_end += usediags ? 8 : 4
    end

    for ix = ix_start+1:ix_end
        # process subquadrant ix
        ix_mod = ((ix-1) % (usediags ? 8 : 4)) + 1
        flags[ix_mod] = true
    end

    # Choose where to add edges.  If there are more 'false' than true, add
    # intersection edges for the false ones, otherwise for the true ones.
    # However, if the previous neighbor is cardinal, it should not be included
    # among the intersection edges to add, and we must choose the other
    # intersection edges consistent with how those for the neighbor were
    # chosen. If next neighbor is cardinal, we should make sure not to introduce
    # an intersection edge to it.

    # check whether any neighbor is one step in a cardinal direction away from
    # curpt
    prev_is_cardinal = (abs(v1[1]) <= 1) && (abs(v1[2]) <= 1)
    next_is_cardinal = (abs(v2[1]) <= 1) && (abs(v2[2]) <= 1)

    if prev_is_cardinal
        # will flip flags if necessary
        _set_flags_consistent_with_prev_point!(flags, curpt, prevpt, prev_isects, usediags)
    elseif sum(flags) > length(flags)/2
        flags = .!flags
    end
    prev_is_cardinal && (flags[ix_start] = false)
    next_is_cardinal && (flags[((ix_end-1) % (usediags ? 8 : 4)) + 1] = false)

    result = Set{Tuple{CartesianIndex{2}, CartesianIndex{2}}}()
    for ix = 1:length(flags)
        if flags[ix]
            offset = _subquad2offset(ix, usediags)
            push!(result, (CartesianIndex(curpt), CartesianIndex(curpt) + offset))
        end
    end
    return result
end

# ----------------------------------------------------------------------------
function _set_flags_consistent_with_prev_point!(flags, curpt, prevpt, prev_isects, usediags)
    @assert prev_isects != nothing
    # find which edge was added for the previous point
    prev_dir = (curpt[1] - prevpt[1], curpt[2] - prevpt[2])
    prev_subquad = _dir2subquadrant(prev_dir, usediags)
    p_offset = _subquad2offset(prev_subquad%(usediags ? 8 : 4) + 1, usediags)

    pflag = (CartesianIndex(prevpt), CartesianIndex(prevpt) + p_offset) in prev_isects
    
    cur_subquad = (prev_subquad - 1 + (usediags ? 4 : 2)) % (usediags ? 8 : 4) + 1

    if flags[(cur_subquad - 2 + (usediags ? 8 : 4)) % (usediags ? 8 : 4) + 1] != pflag
        # flip flags to make them consistent
        flags .= .!flags
    end
    
end

# ----------------------------------------------------------------------------
function _subquad2offset(subquad::Int, usediags::Bool)
    result = nothing
    if usediags
        result = 
            (subquad == 1) ? CartesianIndex(1, 0) :
            (subquad == 2) ? CartesianIndex(1, 1) :
            (subquad == 3) ? CartesianIndex(0, 1) :
            (subquad == 4) ? CartesianIndex(-1, 1) :
            (subquad == 5) ? CartesianIndex(-1, 0) :
            (subquad == 6) ? CartesianIndex(-1, -1) :
            (subquad == 7) ? CartesianIndex(0, -1) :
            (subquad == 8) ? CartesianIndex(1, -1) :
            error("conceptually unreachable")
    else
        result = 
            (subquad == 1) ? CartesianIndex(1, 0) :
            (subquad == 2) ? CartesianIndex(0, 1) :
            (subquad == 3) ? CartesianIndex(-1, 0) :
            (subquad == 4) ? CartesianIndex(0, -1) :
            error("conceptually unreachable")
    end
    return result
end

# ----------------------------------------------------------------------------
function _dir2subquadrant(dir::Tuple{Int, Int}, usediags::Bool)
    # return a number from 1 to 8 (or 4) indicating which subquadrant the direction
    # vector is contained in
    xflag = sign(dir[1]) # -1, 0, or 1
    yflag = sign(dir[2]) # -1, 0, or 1
    dflag = sign(abs(dir[2])-abs(dir[1])) # -1, 0, or 1

    @assert abs(xflag) + abs(yflag) != 0 # zero vector not allowed
    return (xflag > 0 && yflag >= 0) ?
              # quadrant 1
              (usediags ? (dflag == -1 ? 1 : 2) : 1 ) :
           (xflag <= 0 && yflag > 0) ?
              # quadrant 2
              (usediags ? (dflag ==  1 ? 3 : 4) : 2) :
           (xflag < 0 && yflag <= 0) ?        
              # quadrant 3
              (usediags ? (dflag == -1 ? 5 : 6) : 3) :
           (xflag >= 0 && yflag < 0) ?
              # quadrant 4
              (usediags ? (dflag ==  1 ? 7 : 8) : 4) :
              error("conceptually unreachable")
end

# ----------------------------------------------------------------------------
function _grid_segment_interior_intersections(segment_start, segment_end, usediags)
    
    @assert eltype(segment_start) <: Integer
    @assert eltype(segment_end) <: Integer
    startpt = CartesianIndex(segment_start[1], segment_start[2])
    endpt = CartesianIndex(segment_end[1], segment_end[2])
    
    # An edge is described by its tw o neighbor gridcells.
    edges = Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}}()
    
    # vertical edges intersected, not counting start and end point
    n = segment_end[1] - segment_start[1] 
    irange = (1:abs(n)-1).*sign(n)
    
    for ix in irange
        dydx = (segment_end[2] - segment_start[2]) / (segment_end[1] - segment_start[1])
        ypos = dydx * ix
        iy = floor(Int, ypos)
        push!(edges, (startpt + CartesianIndex(ix, iy),
                      startpt + CartesianIndex(ix, iy+1)))
    end

    # horizontal edges intersected, not counting start and endpoint
    n = segment_end[2] - segment_start[2]
    irange = (1:abs(n)-1).*sign(n)

    for iy in irange
        dxdy = (segment_end[1] - segment_start[1]) / (segment_end[2] - segment_start[2])
        xpos = dxdy * iy
        ix = floor(Int, xpos)
        push!(edges, (startpt + CartesianIndex(ix, iy),
                      startpt + CartesianIndex(ix+1, iy)))
    end

    # diagonal intersections (\)
    dx = segment_end[1] - segment_start[1]
    dy = segment_end[2] - segment_start[2]
    n = dx + dy # an integer here
    irange = (1:abs(n)-1).*sign(n)

    for i in irange
        t = i / n
        xpos = t * dx
        ypos = t * dy
        ix, iy = floor(Int, xpos), floor(Int, ypos)
        push!(edges, (startpt + CartesianIndex(ix+1, iy),
                      startpt + CartesianIndex(ix, iy+1)))
    end

    # diagonal intersections (/)
    n = dx - dy # an integer here
    irange = (1:abs(n)-1).*sign(n)
    for i in irange
        t = i / n
        xpos = t * dx
        ypos = t * dy
        ix, iy = floor(Int, xpos), floor(Int, ypos)
        push!(edges, (startpt + CartesianIndex(ix, iy),
                      startpt + CartesianIndex(ix+1, iy+1)))
    end

    return edges
end

# ----------------------------------------------------------------------------
"""
    edgeset2dict(edges)
Convert a set of edges into a dictionary mapping each node to its neighbors.
"""
function edgeset2dict(edges::Set{Tuple{CartesianIndex{2}, CartesianIndex{2}}})
    result = Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}}()
    for e in edges
        # one way
        if haskey(result, e[1])
            push!(result[e[1]], e[2])
        else
            result[e[1]] = [e[2]]
        end
        # other way
        if haskey(result, e[2])
            push!(result[e[2]], e[1])
        else
            result[e[2]] = [e[1]]
        end
    end
    # ensure no duplicity in any vector
    for (k, v) in result
        result[k] = unique(v)
    end

    # return the finished dictionary
    return result
end

# ----------------------------------------------------------------------------
function flatten_small_traps(topography::Matrix{<:Real}, vol_threshold::Real) 
    # Function to flatten/eliminate traps smaller than a given volume threshold, by
    # raising their elevation to the spill point elevation.  

    tstruct = spillanalysis(topography, verbose=true)
    new_topography = copy(topography)

    for i in 1:length(tstruct.trapvolumes)
        if tstruct.trapvolumes[i] < vol_threshold
            new_topography[tstruct.footprints[i]] .= tstruct.spillpoints[i].elevation
        end
    end

    return new_topography
end
    
# ----------------------------------------------------------------------------
""" _active_region_spilltree(tstruct, filled) Given a set of filled traps,
construct a tree graph expressing which regions are actively spilling onto
others. (Regions associated with unfilled traps do not spill anywhere).
"""
function _active_region_spilltree(tstruct, filled)
    @assert length(filled) == length(tstruct.spillpoints) # should be one per trap
    uncovered = copy(filled) # keep track of which traps are already covered
    targets = zeros(Int, numregions(tstruct))

    # loop backwards,  since the spill of higher-level traps take precedence
    for i in length(tstruct.spillpoints):-1:1
        if uncovered[i]
            # this trap is filled and has not already been covered by a higher-level trap
            targets[tstruct.lowest_subtraps_for[i]] .= tstruct.spillpoints[i].downstream_region

            # mark all subtraps of this trap as covered
            subtraps = findall(Graphs.dfs_parents(tstruct.agglomerations, i, dir=:in) .> 0)

            uncovered[subtraps] .= false
        end
    end
            
    # making edges for the tree graph
    edges = [(i, targets[i]) for i in eachindex(targets) if targets[i] > 0]

    # construct the graph 
    g = Graphs.SimpleDiGraph([Graphs.SimpleEdge{Int64}((e[1], e[2])) for e in edges])

    # ensure the graph contains as many nodes as there are regions, even if some are disconnected
    Graphs.add_vertices!(g, numregions(tstruct) - Graphs.nv(g))

    return g
end

function current_upstream_area(tstruct::TrapStructure{<:Real}, point, tstates)

    # first, determine if location is submerged or not
    submerged, puddlecells = _is_submerged(tstruct, point, tstates)

    # determine involved, non-submerged cells
    cells = submerged ?
        [] :
        findall(Graphs.dfs_parents(tstruct.flowgraph, point, dir=:in) .> 0)

    # determine first-order full region contributions
    full_regions = submerged ?
        Set(tstruct.regions[puddlecells]) :
        unique(vcat(tstruct.lowest_subtraps_for[
            findall(x -> x ∈ cells, [x.downstream_region_cell for x in tstruct.spillpoints])]...))

    # determine higher-order full region contributions
    flowtree = _active_region_spilltree(tstruct, tstates[1])
    tull = [findall(Graphs.dfs_parents(flowtree, x, dir=:in) .> 0) for x in full_regions]
    
    upstream_trapset = vcat([ findall(Graphs.dfs_parents(flowtree, x, dir=:in) .> 0)
                              for x in full_regions ]...)

    # list of all the regions whose full area contributes to the current point,
    # either directly or through active spill paths
    all_full_regions = union(upstream_trapset, full_regions)
    
    # determine all cells belonging to any of the contributing full regions

    #regcells = findall(in(all_full_regions).(tstruct.regions[:])) @@@ Slow! 
    mask = in.(tstruct.regions, Ref(Set(all_full_regions)))
    regcells = findall(mask[:])
    
    return union(cells, regcells)
end

