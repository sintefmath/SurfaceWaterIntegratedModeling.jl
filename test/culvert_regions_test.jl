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

# The same passage can reach us twice (recorded both as a culvert and as a piped
# watercourse), and an inlet may be offered two outlets; one path must survive.
@testset "repeated and rival culverts out of one inlet" begin
    G = SurfaceWaterIntegratedModeling.Graphs
    z = [10.0 - j for i in 1:10, j in 1:10]
    z[5, 5] -= 2.0
    inlet, outlet, rival = CartesianIndex(5, 5), CartesianIndex(5, 8), CartesianIndex(9, 7)
    LI = LinearIndices(size(z))

    drained(cv) = length(spillanalysis(z, culverts=cv).trapvolumes)
    @test drained([(inlet, outlet)]) == 0
    @test drained([(inlet, outlet), (inlet, outlet)]) == 0          # exact repeat
    @test drained([(inlet, outlet), (outlet, inlet)]) == 0          # same pair, reversed

    # a rival culvert out of the same inlet is ignored, the first one stands
    _, g, _ = spillregions(spillfield(z)[1], grid = z,
                           culverts = [(inlet, outlet), (inlet, rival)])
    @test G.has_edge(g, LI[inlet], LI[outlet])
    @test !G.has_edge(g, LI[inlet], LI[rival])

    # a culvert along the cell's own downslope direction must not cancel itself
    plain = [10.0 - j for i in 1:10, j in 1:10]
    src, dst = CartesianIndex(5, 5), CartesianIndex(5, 6)
    _, g2, _ = spillregions(spillfield(plain)[1], grid = plain, culverts = [(src, dst)])
    @test G.has_edge(g2, LI[src], LI[dst])
end
