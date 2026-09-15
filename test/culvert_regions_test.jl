# Culverts and spill regions (spillregions.jl): a culvert inlet is excluded from
# its region's bottom cells, so a region draining through a culvert alone has no
# bottom cells left to look up.
#
using Test, SurfaceWaterIntegratedModeling

@testset "culvert inlet as a region's only bottom cell" begin
    # terrain sloping down in +j, with a one-cell pit that is the culvert inlet
    z = [10.0 - j for i in 1:10, j in 1:10]
    z[5, 5] -= 2.0
    inlet, outlet = CartesianIndex(5, 5), CartesianIndex(5, 8)

    field, _ = spillfield(z)
    @test findall(field .== -1) == [inlet]  # the pit is the only bottom cell
    @test z[inlet] > z[outlet]              # sloping, so the culvert is a leak edge

    @test length(spillanalysis(z).trapvolumes) == 1
    # the culvert drains the pit, leaving no trap
    @test length(spillanalysis(z, culverts=[(inlet, outlet)]).trapvolumes) == 0
end
