# Static analysis: input handling of spillanalysis / spillfield.
#
using Test, SurfaceWaterIntegratedModeling

# a small terrain with two pits, sloping down towards one edge
_two_pit_terrain() = (z = [10.0 + 0.5j + 0.3i for i in 1:12, j in 1:15];
                      z[4, 5] -= 3.0; z[8, 10] -= 2.5; z)

@testset "spillanalysis: grid lengths (Lx, Ly)" begin
    z = _two_pit_terrain()
    ts = spillanalysis(z)
    # lengths equal to the cell counts mean square cells, the default
    ts_sq = spillanalysis(z; lengths = size(z))
    @test ts_sq.regions == ts.regions
    @test ts_sq.trapvolumes == ts.trapvolumes
    dir, = spillfield(z; lengths = Float64.(size(z)))
    @test dir == spillfield(z)[1]
    # stretching the cells in one direction changes which way the slope is steepest
    @test spillfield(z; lengths = (12.0, 150.0))[1] != dir
end
