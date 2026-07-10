using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# A two-cell-tall, two-cell-wide footprint centred in a grid of size `sz`.
function _center_footprint(sz)
    li = LinearIndices(sz)
    r  = sz[1] ÷ 2
    c  = sz[2] ÷ 2
    return [li[r, c], li[r, c + 1], li[r + 1, c], li[r + 1, c + 1]]
end

# Hand-built regions/topography/spillpoints for deterministic outlet-validation
# tests, independent of terrain routing.  Region 1 = the NBS trap (top-left 2x2
# of a 3x3 grid), region 2 = everywhere else.  The trap's spillpoint sits at
# elevation 5.0; its downstream (discharge) cell is (3,3), in region 2.  The two
# regions share no supertrap.
function _synthetic_nbs_case()
    regions = [1 1 2;
               1 1 2;
               2 2 2]
    # topography: trap cells low; a low downstream cell (3,3)=1.0 and a high one (1,3)=9.0
    topo = [0.0 0.0 9.0;
            0.0 0.0 3.0;
            2.0 2.0 1.0]
    li = LinearIndices(size(regions))
    footprint = [li[1, 1], li[1, 2], li[2, 1], li[2, 2]]
    # one Spillpoint per region; trap (region 1) spills at elevation 5.0, from
    # trap-side cell (1,1) to downstream-side cell (3,3).
    sp = [SWIM.Spillpoint{Float64}(2, li[1, 1], li[3, 3], 5.0),
          SWIM.Spillpoint{Float64}(0, li[3, 3], 0, Inf)]
    supertraps_of = [[1], [2]]   # each region is its own top-level trap; no sharing
    return regions, topo, footprint, sp, supertraps_of
end

# ---------------------------------------------------------------------------
@testset "NBS static: end-to-end on a real grid" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    foot = _center_footprint(size(grid))
    # single-layer puddle -> exactly one outlet; leave it unspecified so it is
    # backfilled to the trap's natural spillpoint (no terrain-specific knowledge
    # needed for the test to stay deterministic).
    p = SWIM.DynNBSPlacement(SWIM.puddle(5.0), foot, [CartesianIndex(0, 0)])
    t = spillanalysis(grid; nbs = [p])

    reg = t.regions[foot[1]]
    @test reg > 0
    @test all(t.regions[c] == reg for c in foot)   # footprint is exactly one trap
    @test t.trapvolumes[reg] > 0                    # dug trap has positive volume
    @test length(t.nbs) == 1
    @test t.nbs[1].system.layers[1].A == 4.0        # area overwritten from footprint
    outlet = t.nbs[1].outlets[1]
    @test outlet != CartesianIndex(0, 0)            # outlet was backfilled
    @test t.regions[outlet] != reg                  # and it is outside the NBS region
end

# ---------------------------------------------------------------------------
@testset "NBS _prepare_nbs!" begin
    grid = zeros(4, 4)
    li = LinearIndices(size(grid))

    # area overwrite (footprint = 3 cells -> A = 3.0)
    p = SWIM.DynNBSPlacement(SWIM.puddle(5.0), [li[1, 1], li[1, 2], li[2, 1]],
                          [CartesianIndex(3, 3)])
    SWIM._prepare_nbs!([p], grid)
    @test p.system.layers[1].A == 3.0

    # empty footprint -> error
    @test_throws ErrorException SWIM._prepare_nbs!(
        [SWIM.DynNBSPlacement(SWIM.puddle(5.0), Int[], [CartesianIndex(1, 1)])], grid)

    # out-of-bounds footprint cell -> error
    @test_throws ErrorException SWIM._prepare_nbs!(
        [SWIM.DynNBSPlacement(SWIM.puddle(5.0), [999], [CartesianIndex(1, 1)])], grid)

    # outlet count != layer count -> error (puddle has 1 layer, 2 outlets given)
    @test_throws ErrorException SWIM._prepare_nbs!(
        [SWIM.DynNBSPlacement(SWIM.puddle(5.0), [li[1, 1]],
                           [CartesianIndex(1, 1), CartesianIndex(2, 2)])], grid)
end

# ---------------------------------------------------------------------------
@testset "NBS _dig_nbs_traps!" begin
    grid = Float64[i + j for i in 1:5, j in 1:5]   # min = 2.0 at (1,1)
    li = LinearIndices(size(grid))
    foot = [li[3, 3], li[3, 4]]
    p = SWIM.DynNBSPlacement(SWIM.puddle(5.0), foot, [CartesianIndex(1, 1)])
    level = SWIM._dig_nbs_traps!(grid, p |> x -> [x])
    @test level == minimum(Float64[i + j for i in 1:5, j in 1:5]) - SWIM.NBS_DIG_DROP
    @test all(grid[c] == level for c in foot)      # footprint lowered
    @test grid[li[1, 1]] == 2.0                    # untouched cell unchanged
end

# ---------------------------------------------------------------------------
@testset "NBS _validate_nbs_outlet" begin
    regions, topo, _, sp, _ = _synthetic_nbs_case()
    sp_elev = sp[1].elevation                       # 5.0

    # inside the NBS trap's own region -> error [§25]
    @test_throws ErrorException SWIM._validate_nbs_outlet(
        1, "test", CartesianIndex(1, 1), 1, regions, topo, sp_elev)

    # outside region but NOT strictly below the spillpoint elevation -> error [§15]
    # (1,3) is region 2, topography 9.0 >= 5.0
    @test_throws ErrorException SWIM._validate_nbs_outlet(
        1, "test", CartesianIndex(1, 3), 1, regions, topo, sp_elev)

    # outside region and strictly below the spillpoint -> ok (3,3) region 2, topo 1.0
    @test SWIM._validate_nbs_outlet(
        1, "test", CartesianIndex(3, 3), 1, regions, topo, sp_elev) === nothing

    # out of bounds -> error
    @test_throws ErrorException SWIM._validate_nbs_outlet(
        1, "test", CartesianIndex(9, 9), 1, regions, topo, sp_elev)
end

# ---------------------------------------------------------------------------
@testset "NBS _resolve_nbs! outlet backfill" begin
    regions, topo, footprint, sp, sto = _synthetic_nbs_case()

    # single-layer, unspecified outlet -> backfilled to the downstream discharge cell (3,3)
    p1 = SWIM.DynNBSPlacement(SWIM.puddle(5.0), footprint, [CartesianIndex(0, 0)])
    SWIM._resolve_nbs!([p1], regions, sp, topo, sto)
    @test p1.outlets[1] == CartesianIndex(3, 3)

    # two-layer, both unspecified -> lowermost = discharge cell, upper = lowermost outlet
    twolayer = SWIM.elhadiGreenRoof(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0)
    p2 = SWIM.DynNBSPlacement(twolayer, footprint,
                           [CartesianIndex(0, 0), CartesianIndex(0, 0)])
    SWIM._resolve_nbs!([p2], regions, sp, topo, sto)
    @test p2.outlets[2] == CartesianIndex(3, 3)     # lowermost -> downstream discharge cell
    @test p2.outlets[1] == p2.outlets[2]            # upper -> lowermost outlet

    # specified valid outlet is preserved
    p3 = SWIM.DynNBSPlacement(SWIM.puddle(5.0), footprint, [CartesianIndex(3, 3)])
    SWIM._resolve_nbs!([p3], regions, sp, topo, sto)
    @test p3.outlets[1] == CartesianIndex(3, 3)

    # specified outlet inside the region -> error
    p4 = SWIM.DynNBSPlacement(SWIM.puddle(5.0), footprint, [CartesianIndex(2, 2)])
    @test_throws ErrorException SWIM._resolve_nbs!([p4], regions, sp, topo, sto)
end

# ---------------------------------------------------------------------------
@testset "NBS supertrap cycle guard" begin
    regions, topo, footprint, sp, _ = _synthetic_nbs_case()

    # outlet (3,3) is region 2, strictly below the spillpoint and out of the NBS
    # region — but here regions 1 and 2 share supertrap 3, so it must be rejected.
    shared = [[1, 3], [2, 3]]
    p = SWIM.DynNBSPlacement(SWIM.puddle(5.0), footprint, [CartesianIndex(3, 3)])
    @test_throws ErrorException SWIM._resolve_nbs!([p], regions, sp, topo, shared)

    # the direct helper: shared supertrap -> error; disjoint -> ok
    @test_throws ErrorException SWIM._check_nbs_no_shared_supertrap(
        1, 1, CartesianIndex(3, 3), 1, regions, shared)
    @test SWIM._check_nbs_no_shared_supertrap(
        1, 1, CartesianIndex(3, 3), 1, regions, [[1], [2]]) === nothing
end
