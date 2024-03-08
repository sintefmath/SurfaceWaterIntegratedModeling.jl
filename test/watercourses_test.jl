using SurfaceWaterIntegratedModeling

function test_watercourses(grid, precipitation, filled_traps, infiltration)

    tstruct = spillanalysis(grid, usediags=true);

    precipitation = typeof(precipitation) <: Matrix ?
        precipitation :
        ones(size(grid)) * precipitation

    infiltration = typeof(infiltration) <: Matrix ?
        infiltration :
        ones(size(grid)) * infiltration
    

    
    runoff_total, regs, outside, infil =
        watercourses(tstruct, Vector{Bool}(filled_traps), precipitation, infiltration)

    total_precipitation = sum(precipitation)

    return total_precipitation ≈ sum(regs) + outside + infil        
    
end
