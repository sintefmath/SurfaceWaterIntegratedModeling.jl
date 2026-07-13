using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# Gate Phase C: the incremental network driver (build / predict / per-event step /
# finalize).  Uses a self-contained synthetic terrain — a plane rising with column
# index (water flows west, off the j=1 edge) carved with one interior pit — so the
# single trap fills and spills straight out of the domain (spill_path == -1): a clean
# single-basin evolution with no downstream absorption.  (Multi-trap absorption parity
# needs fill_sequence's static loop to keep cur_amounts current; that lands in the
# Phase D/E integration tests.)
# ---------------------------------------------------------------------------
function _pit_terrain(N = 30)
    grid = Float64[j for i in 1:N, j in 1:N]
    for i in 12:18, j in 12:18
        grid[i, j] = 2.0
    end
    return spillanalysis(grid, usediags = true)
end

function _driver_inputs(ts; rain_rate = 1.0e-3)
    filled   = Vector{Bool}(ts.trapvolumes .== 0.0)
    sgraph   = SWIM.compute_complete_spillgraph(ts, filled)
    rain     = fill(rain_rate, size(ts.topography))
    infil    = zeros(size(ts.topography))
    rateinfo = SWIM.compute_flow(sgraph, rain, infil, ts, false)
    zvt      = SWIM._compute_z_vol_tables(ts)
    cur      = fill(SWIM.FilledAmount(0.0, 0.0), numtraps(ts))
    return (; filled, infil, rateinfo, zvt, cur)
end

@testset "network driver: build + predict" begin
    ts = _pit_terrain()
    @test numtraps(ts) == 1
    inp = _driver_inputs(ts)
    driver = SWIM.build_network_driver(ts, [CartesianIndex(15, 15)], SWIM.DynCulvert[],
                                       findall(inp.filled), inp.cur, inp.rateinfo,
                                       inp.infil, inp.zvt, 0.0, Inf)
    @test length(driver.comps) == 1
    @test length(driver.contexts) == 1
    ev = driver.contexts[1].next_event
    @test ev.kind == :fill                       # the pit fills under steady rain
    @test isfinite(ev.time) && ev.time > 0
    # at build the pit is the transitory frontier leaf (no spill path until it fills)
    @test driver.comps[1].traps[1].spill_path == 0
end

@testset "network driver: _driver_next_event" begin
    ts  = _pit_terrain()
    inp = _driver_inputs(ts)
    driver = SWIM.build_network_driver(ts, [CartesianIndex(15, 15)], SWIM.DynCulvert[],
                                       findall(inp.filled), inp.cur, inp.rateinfo,
                                       inp.infil, inp.zvt, 0.0, Inf)
    sel = SWIM._driver_next_event(driver)
    @test sel !== nothing
    @test sel[1] == 1 && sel[2].kind == :fill
end

@testset "network driver: step + finalize (single basin)" begin
    ts  = _pit_terrain()
    inp = _driver_inputs(ts)
    driver = SWIM.build_network_driver(ts, [CartesianIndex(15, 15)], SWIM.DynCulvert[],
                                       findall(inp.filled), inp.cur, inp.rateinfo,
                                       inp.infil, inp.zvt, 0.0, Inf)
    filled_traps = copy(inp.filled)
    lastT = -Inf
    nev   = 0
    while true
        r = SWIM.step_network_driver!(driver, ts, inp.rateinfo, inp.infil, inp.zvt,
                                      filled_traps, inp.cur, Inf)
        r === nothing && break
        nev += 1
        @test r.time >= lastT                    # strictly time-ordered
        lastT = r.time
        @test nev <= 5                           # single basin: no runaway
    end
    @test nev == 1                               # one :fill, then steady out-of-domain spill
    @test filled_traps[driver.comps[1].traps[1].trap_ix]   # trap now full

    SWIM.finalize_network_driver!(driver, inp.cur, ts, inp.infil, inp.zvt, lastT + 1.0e5)
    g = driver.comps[1].traps[1].trap_ix
    @test inp.cur[g].amount ≈ SWIM._own_capacity(ts, g)    # settled at capacity
    @test all(isfinite(ca.amount) for ca in inp.cur)
end

# ---------------------------------------------------------------------------
# Gate Phase D: the driver IS the fill_sequence network path (D2 retired the old one), on
# the real mini.txt terrain.  A `dyn_traps` network spanning a multi-trap chain that
# terminates out of the domain exercises the incremental membership layer end-to-end (build →
# per-event apply → reconcile → finalize) and the full-terminal `-1` spill-path the tracer
# must assign.
# ---------------------------------------------------------------------------
@testset "network driver: through fill_sequence (D1)" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    ts   = spillanalysis(grid, usediags = true)
    nt   = numtraps(ts)
    we   = [SurfaceWaterIntegratedModeling.WeatherEvent(0.0, 1.0e-3)]

    function fold_amounts(seq, n)
        amt = fill(0.0, n)
        for ev in seq
            a = ev.amount
            if eltype(a) <: SWIM.FilledAmount
                for (i, fa) in enumerate(a); amt[i] = fa.amount; end
            else
                for u in a; amt[u.index] = u.value.amount; end
            end
        end
        return amt
    end
    monotone(seq) = all(seq[i].timestamp <= seq[i+1].timestamp + 1e-9 for i in 1:length(seq)-1)
    caps = [SWIM._own_capacity(ts, t) for t in 1:nt]

    # 1. no dyn_traps/culverts → the driver builds an empty network and the run is unchanged
    #    (the plain fill_sequence result); a full snapshot per weather period plus events
    plain = fill_sequence(ts, we)
    @test !isempty(plain)
    @test all(isfinite(a) for a in fold_amounts(plain, nt))

    # 2. a dynamic network at trap 233 (a multi-trap chain terminating out of domain) runs
    #    end-to-end through the driver with all invariants intact
    seq = fill_sequence(ts, we; dyn_traps = [233])
    amt = fold_amounts(seq, nt)
    @test monotone(seq)
    @test all(isfinite, amt)
    @test all(amt[t] <= caps[t] + 1e-6 for t in 1:nt)     # never over capacity
    @test all(amt[t] >= -1e-6 for t in 1:nt)              # never negative
    @test amt[233] ≈ caps[233] rtol=1e-6                   # the seeded trap fills
end

# ---------------------------------------------------------------------------
# Gate Phase E2: an NBS placement carried end-to-end through the driver.  A west-flowing
# plane with a pit near the west edge; a slow-draining `puddle` NBS sits upstream on the
# flow path into the pit.  The overlay must run through fill_sequence (build → per-event
# apply → finalize) and visibly DELAY the downstream pit's fill by holding runoff back.
# ---------------------------------------------------------------------------
@testset "network driver: NBS through fill_sequence (E2)" begin
    N    = 30
    grid = Float64[j for i in 1:N, j in 1:N]              # rising east → flow west
    for i in 12:18, j in 3:6; grid[i, j] = 1.0; end       # pit near the west edge
    ts   = spillanalysis(grid, usediags = true)
    LI   = LinearIndices(size(grid))
    @test numtraps(ts) == 1
    pit  = 1

    # NBS footprint upstream (east) of the pit, on the flow path into it; a slow overflow
    # (kOUT small) so it stores runoff and releases it gradually.
    fp = [LI[CartesianIndex(i, j)] for i in 12:18 for j in 15:20]
    pl = SurfaceWaterIntegratedModeling.DynNBSPlacement(
             SurfaceWaterIntegratedModeling.puddle(50.0; kOUT = 0.01), fp,
             CartesianIndex{2}[])
    w  = [SurfaceWaterIntegratedModeling.WeatherEvent(0.0, 1.0e-3)]

    filltimes(seq) = (d = Dict{Int,Float64}();
        for e in seq[2:end], u in e.filled
            (u isa IncrementalUpdate) && u.value && !haskey(d, u.index) && (d[u.index] = e.timestamp)
        end; d)
    monotone(seq) = all(seq[i].timestamp <= seq[i+1].timestamp + 1e-9 for i in 1:length(seq)-1)

    seq0 = fill_sequence(ts, w)                            # no NBS
    seqN = fill_sequence(ts, w; nbs = [pl])                # NBS on the flow path
    f0, fN = filltimes(seq0), filltimes(seqN)

    @test monotone(seqN)                                  # runs end-to-end, ordered
    @test haskey(f0, pit) && haskey(fN, pit)              # the pit fills either way
    @test fN[pit] > f0[pit] * 1.2                         # the upstream NBS delays it markedly
end
