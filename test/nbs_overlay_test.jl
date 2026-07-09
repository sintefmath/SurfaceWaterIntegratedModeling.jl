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
    mk(fp) = NBSPlacement(puddle(10.0), fp, CartesianIndex{2}[])

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
end
