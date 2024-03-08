using SurfaceWaterIntegratedModeling

function test_spillfield(grid, num_trapcells, num_buildingcells, num_sinkcells;
                         usediags=true, lengths=nothing, domain=nothing, tiling=nothing)

    field, slope = spillfield(grid, usediags=usediags, lengths=lengths, domain=domain, tiling=tiling)

    return length(findall(field.==-1)) == num_trapcells &&
           length(findall(field.==-2)) == num_buildingcells &&
           length(findall(field.==-3)) == num_sinkcells
end


function test_spillregions(grid, minreg, maxreg, outregsize;
                           usediags=true, tiling=nothing)

    field, slope = spillfield(grid, usediags=usediags)

    regions = spillregions(field, usediags=usediags, tiling=tiling)

    exreg = extrema(regions)
    oreg = findall(regions .< 0)

    return minreg == exreg[1] &&
           maxreg == exreg[2] &&
           length(oreg) == outregsize
end

    
function test_update_spillfield_changes(grid, domain)

    # compute solution on unmodified grid
    field, slope = spillfield(grid)
    field_hash = hash(field)

    # modify grid
    modif = view(grid, domain.xrange, domain.yrange)
    modif .= modif'

    update_spillfield!(field, slope, grid, domain)
    field2_hash = hash(field)

    return field2_hash != field_hash
    
end

function test_update_spillfield_sameanswer(grid, domain)

    # compute solution on unmodified grid
    field, slope = spillfield(grid)
    field_hash = hash(field)

    # modify grid
    modif = view(grid, domain.xrange, domain.yrange)
    modif .= modif'

    xrange = domain.xrange
    yrange = domain.yrange
    domain = Domain2D(xrange[1]-1:xrange[end]+1,
                      yrange[1]-1:yrange[end]+1)
    
    update_spillfield!(field, slope, grid, domain)
    field2_hash = hash(field)

    # compute result again, globally
    field, slope = spillfield(grid)

    field3_hash = hash(field)

    return field3_hash == field2_hash
    
end

function test_update_regions_changes(grid, domain)

    # compute solution on unmodified grid
    field, slope = spillfield(grid)
    regions = spillregions(field)

    reg_hash = hash(regions)
    
    # modify grid and spillfield
    modif = view(grid, domain.xrange, domain.yrange)
    modif .= modif'
    update_spillfield!(field, slope, grid, domain)

    # update regions with modification
    update_spillregions!(regions, field, domain)
    reg2_hash = hash(regions)

    return reg_hash != reg2_hash
end

function test_update_regions_sameanswer(grid_orig, domain)

    # compute solution on unmodified grid

    grid_modif = copy(grid_orig)
    modif = view(grid_modif, domain.xrange, domain.yrange)
    modif .= modif'; # transposing
    
    # computing original spill field and region
    field_orig, slope_orig = spillfield(grid_orig)
    regions_orig = spillregions(field_orig)
    
    reg_hash = hash(regions_orig)
    
    # computing modified spillfield
    field_modif, slope_modif = spillfield(grid_modif)

    # compute update locally
    regions_modif = copy(regions_orig)
    update_spillregions!(regions_modif, field_modif, domain)
    reg2_hash = hash(regions_modif)

    # compute update globally
    regions_modif_global = spillregions(field_modif)
    reg3_hash = hash(regions_modif_global)
    
    return reg3_hash == reg2_hash
end

function test_return_region_reindex(grid_orig, domain)

    # compute solution on unmodified grid
    grid_modif = copy(grid_orig)
    modif = view(grid_modif, domain.xrange, domain.yrange)
    modif .= modif'; # transposing
    
    # computing original spill field and region
    field_orig, slope_orig = spillfield(grid_orig)
    regions_orig = spillregions(field_orig)
    
    # computing modified spillfield
    field_modif, slope_modif = spillfield(grid_modif)

    # compute update locally
    regions_modif = copy(regions_orig)
    reindex = update_spillregions!(regions_modif, field_modif, domain,
                                   usediags=true, return_region_reindex=true)
    
    # compute reindex by brute force
    reindex_bf = unique(hcat(regions_orig[:], regions_modif[:]), dims=1)

    return sort(reindex, dims=1) == sort(reindex, dims=1)
    
end

function test_spillpoint_tiling(grid; usediags=true, tiling=nothing)

    field, slope = spillfield(grid, usediags=usediags)
    regions = spillregions(field, usediags=usediags)

    spoints, = spillpoints(grid, regions, usediags=usediags)
    spoints_tiling, = spillpoints(grid, regions, usediags=usediags, tiling=tiling)

    elevation = [x.elevation for x in spoints]
    elevation_tiling = [x.elevation for x in spoints_tiling]

    cell = [x.current_region_cell for x in spoints]
    cell_tiling = [x.current_region_cell for x in spoints_tiling]


    elev_diff = unique(elevation .- elevation_tiling)
    cell_diff = unique(cell .- cell_tiling)

    # verify that the outcome is the same in the two cases (i.e. that there's
    # only a single value for the difference, and that value is zero).
    return length(elev_diff) == 1 && elev_diff[1] == 0.0 &&
           length(cell_diff) == 1 && cell_diff[1] == 0
end

function test_buildings_and_flattening(grid, mask, target_nums; usediags=true)

    water_surfaces = identify_flat_areas(grid, 1e-3, 1000)

    grid_flat = copy(grid)
    flatten_grid!(grid_flat, water_surfaces, :min)

    tstruct  = spillanalysis(grid, usediags=usediags);
    tstruct2 = spillanalysis(grid_flat, building_mask=mask, usediags=usediags);

    return length(unique(tstruct.regions)) == target_nums[1] &&
           length(unique(tstruct2.regions)) == target_nums[2]
end

function test_sinks(grid, sinks; usediags=true)

    tstruct = spillanalysis(grid, usediags=usediags);
    tstruct2 = spillanalysis(grid, sinks=sinks, usediags=usediags);

    tvols = trapvolumes(grid, tstruct.regions, tstruct.spillpoints, tstruct.lowest_subtraps_for);
    tvols2 = trapvolumes(grid, tstruct2.regions, tstruct2.spillpoints, tstruct2.lowest_subtraps_for);

    tot_vol = sum(tvols[toplevel_traps(tstruct.agglomerations)]);
    tot_vol2 = sum(tvols[toplevel_traps(tstruct2.agglomerations)]);

    return tot_vol > tot_vol2
end


function test_spillanalysis(grid, num_spoints, tot_trapvol, tot_subvols, tot_footprint_area,
                            num_agglom_edges;
                            usediags=true, building_mask=nothing, merge_outregs=false)

    tstruct = spillanalysis(grid, usediags=usediags, building_mask=building_mask,
                            merge_outregions=merge_outregs)

    length(tstruct.spillpoints) == num_spoints                         || return false
    sum(tstruct.trapvolumes) ≈ tot_trapvol                             || return false
    sum(tstruct.subvolumes) ≈ tot_subvols                              || return false
    sum([length(x) for x in tstruct.footprints]) == tot_footprint_area || return false
    tstruct.agglomerations.ne ==  num_agglom_edges                     || return false
    
    return true
end

function test_upstream_area(grid, point, area; local_only=true)

    tstruct = spillanalysis(grid)
    return area == length(upstream_area(tstruct, point, local_only=local_only))
end


function test_sequencing(grid, use_infiltration, seqlength, maxtime; mask=nothing)

    tstruct = spillanalysis(grid, building_mask=mask)

    infiltration = nothing
    weather_events = [WeatherEvent(0.0, 1.0)]

    if use_infiltration
        infiltration = 0.5 * ones(size(grid)...)
        weather_events = [WeatherEvent(0.0, 1.0), WeatherEvent(1.0, 0.0)]
    end
    
    seq = fill_sequence(tstruct, weather_events, infiltration=infiltration, verbose=false)

    length(seq) == seqlength || return false
    seq[end].timestamp ≈ maxtime || return false

    return true
end

