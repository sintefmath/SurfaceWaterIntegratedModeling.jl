using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# network_states_at_timepoints: replay the coupled network from `seq` and read off
# exact network trap volumes + NBS layer volumes at requested timepoints.
# Scenario mirrors the E2 driver test: a west-flowing plane with a pit near the
# west edge and a slow `puddle` NBS on the flow path into it.
# ---------------------------------------------------------------------------
@testset "network_states_at_timepoints: NBS layer + trap volumes" begin
    N    = 30
    grid = Float64[j for i in 1:N, j in 1:N]
    for i in 12:18, j in 3:6; grid[i, j] = 1.0; end
    ts   = spillanalysis(grid, usediags = true)
    LI   = LinearIndices(size(grid))
    @test numtraps(ts) == 1
    pit  = 1

    fp = [LI[CartesianIndex(i, j)] for i in 12:18 for j in 15:20]
    pl = SWIM.DynNBSPlacement(SWIM.puddle(50.0; kOUT = 0.01), fp, CartesianIndex{2}[])
    w  = [SWIM.WeatherEvent(0.0, 1.0e-3)]

    seq = fill_sequence(ts, w; nbs = [pl])
    @test !isempty(seq)

    # pick sample times spanning the run: the event timestamps themselves plus midpoints
    evt_times = [e.timestamp for e in seq]
    tmax = maximum(evt_times) * 1.5 + 1.0
    tps  = sort(unique(vcat(evt_times, [(evt_times[k] + evt_times[k+1]) / 2
                                        for k in 1:length(evt_times)-1], tmax)))
    tps  = filter(t -> t >= seq[1].timestamp, tps)

    res = network_states_at_timepoints(ts, seq, tps; nbs = [pl], verbose = false)
    @test length(res) == length(tps)

    id = pl.id
    for (k, tp) in enumerate(tps)
        r = res[k]
        @test haskey(r.nbs_layers, id)                   # NBS layers reported every timepoint
        layers = r.nbs_layers[id]
        @test length(layers) == 1                        # puddle == single layer
        @test all(isfinite, layers)
        @test layers[1] >= -1e-9                          # non-negative storage
        @test haskey(r.trap_volumes, pit)                # the pit is network-covered
        @test isfinite(r.trap_volumes[pit])
        @test r.trap_volumes[pit] >= -1e-9
    end

    # NBS storage should build up then eventually drain: strictly positive at some interior time
    @test any(res[k].nbs_layers[id][1] > 1e-6 for k in 1:length(tps))

    # fidelity: at each event timestamp, the replayed pit volume matches what `seq` recorded
    # (network-covered trap amount is stored in the event stream)
    function amount_at_time(seq, t, trapix)
        a = 0.0
        for e in seq
            e.timestamp > t + 1e-12 && break
            fld = e.amount
            if eltype(fld) <: SWIM.FilledAmount
                a = fld[trapix].amount
            else
                for u in fld; u.index == trapix && (a = u.value.amount); end
            end
        end
        return a
    end
    for (k, tp) in enumerate(tps)
        any(abs(tp - et) < 1e-9 for et in evt_times) || continue
        expected = amount_at_time(seq, tp, pit)
        @test isapprox(res[k].trap_volumes[pit], expected; atol = 1e-4, rtol = 1e-3)
    end

    # the combiner overlays the network volumes onto the single-trap result and carries layers
    both = all_states_at_timepoints(ts, seq, tps; nbs = [pl], verbose = false)
    plain = trap_states_at_timepoints(ts, seq, tps; verbose = false)
    @test length(both) == length(tps)
    for k in 1:length(tps)
        @test both[k].amounts[pit] == res[k].trap_volumes[pit]   # network volume wins for the pit
        @test both[k].nbs_layers[id] == res[k].nbs_layers[id]    # NBS layers carried through
        @test both[k].filled == plain[k][1]                       # filled status from trap_states
    end
end
