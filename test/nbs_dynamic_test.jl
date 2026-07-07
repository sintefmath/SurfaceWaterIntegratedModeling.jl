using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

cidx(i, j) = CartesianIndex(i, j)

# ---------------------------------------------------------------------------
@testset "NBS _components join via nbs_links" begin
    # Two independent path->trap pairs with no shared cells and no natural
    # connection: they form two components on their own.
    paths = [SWIM.DynFlowPath([cidx(1, 1)], 1),
             SWIM.DynFlowPath([cidx(5, 5)], 2)]
    traps = [SWIM.DynTrap(10, 0), SWIM.DynTrap(20, 0)]

    no_links = SWIM._components(paths, traps, Tuple{Symbol,Int}[], Tuple{Symbol,Int}[])
    @test length(no_links) == 2

    # An NBS link from trap 1 to trap 2 must fuse them into one component.
    links = [((:trap, 1), (:trap, 2))]
    joined = SWIM._components(paths, traps, Tuple{Symbol,Int}[], Tuple{Symbol,Int}[], links)
    @test length(joined) == 1
end

# ---------------------------------------------------------------------------
@testset "NBS _nbs_elements + _dyn_seeds" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    sz = size(grid); li = LinearIndices(sz); ci = CartesianIndices(sz)
    r = sz[1] ÷ 2; c = sz[2] ÷ 2
    foot = [li[r, c], li[r, c + 1], li[r + 1, c], li[r + 1, c + 1]]
    p = SWIM.NBSPlacement(SWIM.puddle(5.0), foot, [CartesianIndex(0, 0)])
    t = spillanalysis(grid; nbs = [p])

    nbs_objs = SWIM._nbs_elements(t)
    @test length(nbs_objs) == 1
    @test nbs_objs[1].trap_ix == t.regions[foot[1]]
    @test nbs_objs[1].placement_ix == 1

    seeds = SWIM._dyn_seeds(t, Int[], SWIM.DynCulvert[], nbs_objs)
    @test ci[foot[1]] in seeds                    # NBS footprint seeded
    @test nbs_objs[1].outlets[1] in seeds         # resolved outlet seeded
end

# ---------------------------------------------------------------------------
@testset "NBS network build attaches DynNBS" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    sz = size(grid); li = LinearIndices(sz); ci = CartesianIndices(sz)
    r = sz[1] ÷ 2; c = sz[2] ÷ 2
    foot = [li[r, c], li[r, c + 1], li[r + 1, c], li[r + 1, c + 1]]
    p = SWIM.NBSPlacement(SWIM.puddle(5.0), foot, [CartesianIndex(0, 0)])
    t = spillanalysis(grid; nbs = [p])
    nbs_trap = t.regions[foot[1]]

    nbs_objs = SWIM._nbs_elements(t)
    seeds = SWIM._dyn_seeds(t, Int[], SWIM.DynCulvert[], nbs_objs)
    comps = SWIM.setup_network(t, seeds, Int[]; nbs = nbs_objs)

    # the component holding the NBS trap must carry the DynNBS, and the outlet
    # cell must fall inside that same component.
    host = only(filter(net -> nbs_trap in [tr.trap_ix for tr in net.traps], comps))
    @test length(host.nbs) == 1
    @test host.nbs[1].trap_ix == nbs_trap

    occ = Set{CartesianIndex{2}}()
    for pth in host.flow_paths, cc in pth.cells; push!(occ, cc); end
    for tr in host.traps, k in t.footprints[tr.trap_ix]; push!(occ, ci[k]); end
    @test nbs_objs[1].outlets[1] in occ
end
