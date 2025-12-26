using Test, SurfaceWaterIntegratedModeling
using LazyArtifacts

include("basicTestFuns.jl")
include("watercourses_test.jl")


surf1_file = joinpath(artifact"swim_testdata", "data", "small", "mini.txt")
surf2_file = joinpath(artifact"swim_testdata", "data", "small", "bay.txt")
surf2_mask = joinpath(artifact"swim_testdata", "data", "small", "baybuildings.txt")
surf3_file = joinpath(artifact"swim_testdata", "data", "synthetic", "synsurf.txt")

grid1 = loadgrid(surf1_file) # mini-grid
grid2 = loadgrid(surf2_file) # bay
mask2 = loadgrid(surf2_mask) .|> Bool; # example building mask for bay grid
grid3 = loadgrid(surf3_file) # synthetic surface


@testset "SpillfieldTests" begin
    @test test_spillfield(grid1, 413, 0, 0)
    @test test_spillfield(grid1, 1039, 0, 0) usediags=false

    @test test_spillfield(grid2, 57180, 0, 0)
    @test test_spillfield(grid2, 58071, 0, 0) usediags=false

    @test test_spillfield(grid3, 3, 0, 0)
    @test test_spillfield(grid3, 3, 0, 0) usediags=false
end

@testset "SpillregionTests" begin
    @test test_spillregions(grid1, -242, 413, 2968)
    @test test_spillregions(grid1, -155, 1039, 1915) usediags=false

    @test test_spillregions(grid2, -398, 896, 98164)
    @test test_spillregions(grid2, -216, 1767, 96504) usediags=false

    @test test_spillregions(grid3, -200, 3, 30737)
    @test test_spillregions(grid3, -200, 3, 31910) usediags=false

end

@testset "UpdateGridTests" begin
    grid_copy = copy(grid1);
    @test test_update_spillfield_changes(grid_copy, Domain2D(80:100, 80:100));
    @test test_update_spillfield_sameanswer(grid_copy, Domain2D(80:100, 80:100));
    @test test_update_regions_changes(grid_copy, Domain2D(80:100, 80:100));
    @test test_update_regions_sameanswer(grid_copy, Domain2D(80:100, 80:100));
    @test test_return_region_reindex(grid_copy, Domain2D(80:100, 80:100));
end

# @testset "ParallelismTests" begin
#     @test test_spillregions(grid1, 0xdff242ae3a512d5e) tiling=(1,1)
#     @test test_spillregions(grid1, 0xdff242ae3a512d5e) tiling=(3,3)
#     @test test_spillfield(grid1, 0xa11c55231672f629) tiling=(3,3)
#     @test test_spillpoint_tiling(grid1) tiling=(3,3)
# end

@testset "MaskingTests" begin
    @test test_buildings_and_flattening(grid2, mask2, (1983, 2005)) usediags=false
    @test test_buildings_and_flattening(grid2, mask2, (1294, 1319)) usediags=true
    @test test_sinks(grid2, [CartesianIndex(119,193), CartesianIndex(180, 193)]) usediags=false
    @test test_sinks(grid2, [CartesianIndex(119,193), CartesianIndex(180, 193)]) usediags=true
end

@testset "Watercourses" begin
    # on synthetic grid
    @test test_watercourses(grid3, 1.0, [], 0.0) # empty traps
    @test test_watercourses(grid3, 1.0, [1, 0, 0, 0], 0.0) # subset of traps filled    
    @test test_watercourses(grid3, 1.0, [1, 1, 0, 0], 0.0) # subset of traps filled
    @test test_watercourses(grid3, 1.0, [1, 1, 1, 0], 0.0) # subset of traps filled
    @test test_watercourses(grid3, 1.0, [1, 1, 1, 1], 0.0) # all traps filled
    @test test_watercourses(grid3, 1.0, [1, 1, 0, 0], 0.3) # include infiltration

    # on 'real' grid
    all_traps = ones(1119);
    some_traps = zeros(1119);
    some_traps[1:1000] .= 1.0;
    @test test_watercourses(grid2, 1.0, [], 0.0) # empty traps
    @test test_watercourses(grid2, 1.0, some_traps, 0.0) # subset of traps filled
    @test test_watercourses(grid2, 1.0, all_traps,   0.0) # all traps filled
    @test test_watercourses(grid2, 1.0, some_traps, 0.3) # include infiltration
end

@testset "Trapping structure" begin
    @test test_spillanalysis(grid1, 461, 165.927629, 65.2978599, 1419, 96)
    @test test_spillanalysis(grid1, 1169, 391.2055799, 167.092600, 3167, 260) usediags=false

    @test test_spillanalysis(grid2, 1119, 71139.310390, 56050.730901, 56155, 446)
    @test test_spillanalysis(grid2, 2254, 111771.868559, 96631.510099, 88101, 974) usediags=false

    @test test_spillanalysis(grid3, 4, 4826.5782538, 835.602984, 2494, 2)
    @test test_spillanalysis(grid3, 4, 4826.5782538, 835.602984, 2494, 2) usediags=false
    
end

@testset "Upstream area" begin
    @test test_upstream_area(grid1, 4757, 22)
    @test test_upstream_area(grid1, 4757, 30879) local_only=false

    @test test_upstream_area(grid2, 37574, 394)
    @test test_upstream_area(grid2, 37574, 394) local_only=false

    @test test_upstream_area(grid3, 19958, 4269)
    @test test_upstream_area(grid3, 19958, 5546) local_only=false
    
end

@testset "Sequencing" begin

    @test test_sequencing(grid1, false, 462, 0.04491857143)
    @test test_sequencing(grid2, false, 1090, 0.227784433)
    @test test_sequencing(grid2, false, 1124, 0.231848557) mask=mask2
    @test test_sequencing(grid3, false, 5, 0.4528216)

    @test test_sequencing(grid1, true, 512, 1.51723914) 
    @test test_sequencing(grid2, true, 1313, 5.0519047) 
    @test test_sequencing(grid2, true, 1353, 5.0519047) mask=mask2
    @test test_sequencing(grid3, true, 8, 4.5241752066)

end
