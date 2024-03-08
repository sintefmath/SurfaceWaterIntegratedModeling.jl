import Graphs
import Roots

export flatten_grid!, identify_flat_areas, toplevel_traps,
    show_region_selection, all_subtraps_of, interpolate_timeseries,
    trap_states_at_timepoints, compute_spillfield_graph, all_upstream_regions,
    upstream_area

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
        ownvolumes = tstruct.trapvolumes - tstruct.subvolumes
        remaining =  ownvolumes - tstates[i][2]
        push!(result, _fill_state_to_terrainmap(tstruct, filled_traps, remaining,
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
        # we need these in order to compute Smin for the requested trap
        setsmax!(rateinfo, c, -sum(min.(getrunoff(rateinfo, tstruct.footprints[c]), 0.0)))
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
function _map_all_trapcells(tstruct)
    num_traps = length(tstruct.spillpoints)
    result = [Vector{Int64}() for i = 1:num_traps]
    
    rmap = _regionmap(tstruct)
    num_regions = length(rmap)
    
    for rix = 1:num_regions
        rcells = rmap[rix]
        trap_ix = rix
        while trap_ix > 0
            spill_height = tstruct.spillpoints[trap_ix].elevation
            push!(result[trap_ix], rcells[tstruct.topography[rcells] .<= spill_height]...)

            immediate_supertrap = Graphs.outneighbors(tstruct.agglomerations, trap_ix)
            @assert length(immediate_supertrap) < 2 # should be max. 1 immediate supertrap
            trap_ix = length(immediate_supertrap) == 0 ? -1 : immediate_supertrap[1]
        end
    end
    return result
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
                                   remaining::Vector{<:Real},
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
    trapmap = _map_all_trapcells(tstruct)
    for trap in findall(.!filled)
        if remaining[trap] >= tstruct.trapvolumes[trap] -
                              tstruct.subvolumes[trap] - sqrt(eps())
            # there is no water in the trap yet - nothing to plot
            continue
        end
        subtraps = Graphs.inneighbors(tstruct.agglomerations, trap)
        if all(filled[subtraps])
            trapcells = trapmap[trap]
            FAC = 1-sqrt(eps()) # Avoid problem with bracketing interval due to
                                # roundoff error
            if !isempty(trapcells)
                z_spill = tstruct.spillpoints[trap].elevation
                bottom = minimum(tstruct.topography[trapcells])
                bracket = (bottom, z_spill)
                volfun = z -> sum(z_spill .- max.(z, tstruct.topography[trapcells]))

                z_cur = Roots.find_zero(z -> volfun(z) - remaining[trap]*FAC, bracket)
                covered_trapcells = trapcells[tstruct.topography[trapcells] .< z_cur]
                result[covered_trapcells] .= filled_color
            end
        end
    end
    
    return result
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
function _identify_clusters(mask::Matrix{<:Bool})

    edges = Vector{Tuple{Int, Int}}()
    Dx = CartesianIndex(1, 0)
    Dy = CartesianIndex(0, 1)
    lind = LinearIndices(mask)
    for I in CartesianIndex(1,1):CartesianIndex(size(mask) .- 1)
        if mask[I]
            mask[I + Dx] && push!(edges, (lind[I], lind[I + Dx]))
            mask[I + Dy] && push!(edges, (lind[I], lind[I + Dy]))
        end
    end
    g = Graphs.SimpleGraph([Graphs.SimpleEdge{Int64}(e) for e in edges])
    clusters = Graphs.connected_components(g)
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
    spillfield   = tstruct.spillfield
    regmap       = tstruct.lowest_subtraps_for
   
    if selection == nothing
        # by default, show top-level traps
        selection = toplevel_traps(subtrapgraph)
    end

    result = zeros(Int64, size(regions)...)
    
    # trace rivers
    if river_color != 0
        CI = CartesianIndices(size(grid))
        for reg in selection
            if spillpoints[reg].downstream_region_cell > 0
                rix = CI[spillpoints[reg].downstream_region_cell]; # river startpoint
                finished = false
                while !finished
                    result[rix] = river_color
                    rix, finished = _downstream_cell(spillfield, rix)
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
@inline function _downstream_cell(spillfield::Matrix{Int8}, cix::CartesianIndex)

    dir = spillfield[cix]
    if dir < 0
        return cix, true; # finished
    end

    cix =
        (dir == 0) ? cix + CartesianIndex(-1, 0) :
        (dir == 1) ? cix + CartesianIndex(1, 0) :
        (dir == 2) ? cix + CartesianIndex(0, -1) :
        (dir == 3) ? cix + CartesianIndex(0, 1) :
        (dir == 4) ? cix + CartesianIndex(-1, -1) :
        (dir == 5) ? cix + CartesianIndex(1, 1) :
        (dir == 6) ? cix + CartesianIndex(1, -1) :
        (dir == 7) ? cix + CartesianIndex(-1, 1) :
        cix

    finished = !( (1 <= cix[1] <= size(spillfield, 1)) && (1 <= cix[2] <= size(spillfield, 2)) )

    return cix, finished
    
end

# ----------------------------------------------------------------------------
@inline function _downstream_cell(spillfield::Matrix{Int8}, lix::Int)
    cix, finished = _downstream_cell(spillfield, CartesianIndices(spillfield)[lix])
    return (finished ? -1 : LinearIndices(spillfield)[cix], finished)
end

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
    if in_trap
        if local_only
            return findall(tstruct.regions[:] .== region)
        else
            uregs = all_upstream_regions(tstruct, region)
            num_regions = length(tstruct.supertraps_of)
            lookup = fill(false, num_regions)
            lookup[uregs] .= true
            tmp = fill(false, size(tstruct.topography))
            for i=1:length(tmp[:])
                r = tstruct.regions[i]
                tmp[i] = (r > 0) ? lookup[r] : false
            end
            return findall(tmp[:])            
        end
    else
        # outside a trap, the only upstream area is the one consisting of cells
        # spilling directly into the current cell
        
        # @@ This computation is wasteful and should be restructured later, as
        # we are only interested in a single region.
        sgraph = compute_spillfield_graph(tstruct.spillfield)

        # upstream cells in region
        ucells = findall(Graphs.dfs_parents(sgraph, point, dir=:in) .> 0) 
        
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

See also [`upstream_area`](@ref)
"""
function all_upstream_regions(tstruct::TrapStructure{<:Real},
                              region::Int)
    
    g, vmap, vmap_inv = saturated_spillgraph(tstruct)
    reg_graphnode_ix = vmap_inv[tstruct.supertraps_of[region][end]]
    upstream_graph_nodes = findall(Graphs.dfs_parents(g, reg_graphnode_ix, dir=:in) .> 0)
    upstream_regions = vmap[upstream_graph_nodes]

    # only consider the lowest-level regions
    num_regions = length(tstruct.supertraps_of)
    filter!(x -> x <= num_regions, upstream_regions)

    return upstream_regions
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

    return g
end

