using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# Stage B1: the NBS footprint-as-sink overlay in `watercourses`.  Each footprint
# cell captures the runoff reaching it (rain on the footprint + inflow across its
# boundary) into the owning placement and stops that throughput from continuing
# downstream, so downstream regions no longer double-count it.  See
# agent/NBS_OPTION1_OVERLAY_PLAN.md §1/§2.
# ---------------------------------------------------------------------------
@testset "NBS footprint-as-sink overlay (watercourses)" begin
    # East-rising plane: height = column index, so water flows due west and exits
    # the domain; there are no interior traps, so every cell's rain leaves off-domain.
    N    = 12
    grid = Float64[j for i in 1:N, j in 1:N]
    ts   = spillanalysis(grid)
    LI   = LinearIndices(size(grid))
    nofull = [false]

    # dummy layer model — `watercourses` uses only the footprint, not the system/outlets
    mk(fp) = DynNBSPlacement(puddle(10.0), fp, CartesianIndex{2}[])

    # --- baseline: no NBS -> nothing captured, all rain leaves ------------------
    _, ra0, off0, _, ni0 = watercourses(ts, nofull; precipitation = 1.0, infiltration = 0.0)
    @test isempty(ni0)
    base_down = sum(ra0) + off0
    @test base_down ≈ N * N               # conservation: one unit of rain per cell

    # --- one NBS: a full interior column captures all flow crossing it ----------
    k         = 6
    footprint = Int[LI[i, k] for i in 1:N]
    nbs       = [mk(footprint)]
    _, ra1, off1, _, ni1 =
        watercourses(ts, nofull; precipitation = 1.0, infiltration = 0.0, nbs = nbs)

    @test length(ni1) == 1
    @test ni1[1] > 0                       # the footprint actually captured flow
    # conservation: captured + downstream == baseline downstream
    @test sum(ni1) + sum(ra1) + off1 ≈ base_down
    # no double count: downstream dropped by *exactly* the captured amount
    @test sum(ra1) + off1 ≈ base_down - sum(ni1)

    # --- terrain infiltration under the footprint is ignored (NBS owns it) -------
    infil = zeros(N, N)
    for c in footprint
        infil[c] = 1000.0                  # would swallow everything if it applied
    end
    _, _, _, _, ni_infil =
        watercourses(ts, nofull; precipitation = 1.0, infiltration = infil, nbs = nbs)
    @test ni_infil[1] ≈ ni1[1]

    # --- two disjoint footprints: two tallies, conservation still holds ----------
    nbs2 = [mk(Int[LI[i, 4] for i in 1:N]), mk(Int[LI[i, 9] for i in 1:N])]
    _, ra2, off2, _, ni2 =
        watercourses(ts, nofull; precipitation = 1.0, infiltration = 0.0, nbs = nbs2)
    @test length(ni2) == 2
    @test all(ni2 .> 0)
    @test sum(ni2) + sum(ra2) + off2 ≈ base_down

    # --- overlapping footprints are rejected ------------------------------------
    @test_throws ErrorException watercourses(ts, nofull; nbs = [mk(footprint), mk(footprint)])

    # --- the capture reaches RateInfo via compute_flow (B2 plumbing) ------------
    sg = SWIM.compute_complete_spillgraph(ts, Vector{Bool}(ts.trapvolumes .== 0.0))
    ri = SWIM.compute_flow(sg, 1.0, 0.0, ts, false; nbs = nbs)
    @test length(ri.nbs_inflow) == 1
    @test getnbsinflow(ri, 1) ≈ ni1[1]        # same capture as the direct watercourses call
    # empty by default (no NBS) — existing callers unaffected
    ri0 = SWIM.compute_flow(sg, 1.0, 0.0, ts, false)
    @test isempty(ri0.nbs_inflow)
end

# ---------------------------------------------------------------------------
# Stage B2: the natural (no-NBS) exit split of a footprint, used to spread the
# terrain re-emit of an NBS's top layers over its lower-edge exit boundary.  See
# agent/NBS_OPTION1_OVERLAY_PLAN.md §3.
# ---------------------------------------------------------------------------
@testset "NBS terrain exit weights" begin
    N    = 12
    grid = Float64[j for i in 1:N, j in 1:N]     # east-rising plane: flow due west
    ts   = spillanalysis(grid)
    LI   = LinearIndices(size(grid))
    CI   = CartesianIndices(size(grid))
    footset(fp) = Set(fp)

    # interior column: each cell drains one step west, so every row exits at the
    # column immediately to the west; N equal-weight exits, all outside the footprint.
    col   = Int[LI[i, 6] for i in 1:N]
    w     = SWIM._nbs_exit_weights(ts, col)
    @test sum(x -> x[2], w) ≈ 1.0                 # weights normalised
    @test length(w) == N
    @test all(x -> x[2] ≈ 1 / N, w)               # uniform on a plane
    @test all(x -> x[1] ∉ footset(col), w)        # exits lie outside the footprint
    @test all(x -> CI[x[1]][2] == 5, w)           # all one column to the west

    # a 2-cell-wide block: the upstream column drains through the downstream column,
    # so both exit at the same western boundary — still N equal-weight exits, still
    # summing to one (the internal column is not an exit).
    block = Int[LI[i, j] for i in 1:N for j in (5, 6)]
    w2    = SWIM._nbs_exit_weights(ts, block)
    @test sum(x -> x[2], w2) ≈ 1.0
    @test all(x -> x[1] ∉ footset(block), w2)

    # footprint touching the west domain edge drains off-domain -> sentinel exit 0,
    # which still carries weight so mass is conserved when it is dropped at delivery.
    edge = Int[LI[i, 1] for i in 1:N]
    w3   = SWIM._nbs_exit_weights(ts, edge)
    @test sum(x -> x[2], w3) ≈ 1.0
    @test any(x -> x[1] == 0, w3)                 # off-domain fraction present
end
