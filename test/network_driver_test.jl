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

@testset "network driver: _state_for + _driver_next_event" begin
    ts  = _pit_terrain()
    net = setup_network(ts, findall(Vector{Bool}(ts.trapvolumes .== 0.0));
                        dyn_coords = [CartesianIndex(15, 15)])[1]
    g   = net.traps[1].trap_ix
    # committed volume is taken from the map; an absent node falls back to seed0
    st = SWIM._state_for(net, Dict(g => 3.5), _ -> -1.0)
    @test st == [3.5]
    st2 = SWIM._state_for(net, Dict{Int,Float64}(), tix -> 7.0)
    @test st2 == [7.0]

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
