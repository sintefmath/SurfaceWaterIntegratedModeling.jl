using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# NBS as a signed-diff router element (agent/NBS_ROUTING_NODE_PLAN.md).  `watercourses` is
# NBS-oblivious; the NBS captures its live input `I_1` at the footprint-inflow cells and emits
# only the diff `O_1 - O_0` of its live output over the oblivious baseline, attenuated along
# the carrier path's oblivious runoff by `max(V+d,0) - max(V,0)` per cell and delivered to the
# downstream trap inflow.
# ---------------------------------------------------------------------------

@testset "NBS: _attenuate_diff per-cell rule" begin
    P = SWIM._attenuate_diff
    @test P(Float64[], -3.0)        ≈ -3.0    # empty path: unchanged
    @test P([100.0], -3.0)          ≈ -3.0    # throughput-rich cell: full pass
    @test P([2.0], -5.0)            ≈ -2.0    # runoff can't absorb it all: clamps to -runoff
    @test P([-8.0], -5.0)           ≈  0.0    # cell already infiltrating (V<0): diff dies
    @test P([-8.0],  5.0)           ≈  0.0    # positive diff, spare capacity swallows it
    @test P([-3.0],  5.0)           ≈  2.0    # spare capacity refilled first, surplus continues
    @test P([60.0, 20.0], -30.0)    ≈ -20.0   # chain: attenuates to the smaller downstream runoff
end

@testset "NBS: :nbsin intercepts the whole signed flow" begin
    # A footprint-inflow cell intercepts ALL passing network flow into its live input I_1 (nbs_draw),
    # positive or negative — so an upstream NBS storing (a negative diff riding `current`) shows up as
    # a deficit in a downstream NBS's input, not as flow that bypasses it to the trap.
    D = SWIM._path_delivered!
    net = SWIM.DynNetwork()                          # unused (no culvert events)
    run3 = [10.0, 10.0, 10.0]                        # saturated cells: signed flow passes unattenuated
    ev   = [(2, :nbsin, 1)]                          # one :nbsin at cell 2

    nd = [0.0]                                       # positive release from upstream
    d  = D(run3, 5.0, ev, Float64[], Float64[], nothing, net, nothing, nd)
    @test nd[1] ≈ 5.0 && d ≈ 0.0                     # fully captured, nothing passes on

    nd = [0.0]                                       # negative: upstream NBS storing (deficit)
    d  = D(run3, -3.0, ev, Float64[], Float64[], nothing, net, nothing, nd)
    @test nd[1] ≈ -3.0 && d ≈ 0.0                    # deficit captured (not dropped onto the trap)
end

# A west-flowing plane (rising east) with one pit near the west edge; an NBS footprint sits
# upstream on the flow path into the pit.
function _plane_with_pit(N = 30)
    grid = Float64[j for i in 1:N, j in 1:N]
    for i in 12:18, j in 3:6; grid[i, j] = 1.0; end
    return spillanalysis(grid, usediags = true)
end
_upstream_footprint(ts) =
    (LI = LinearIndices(size(ts.topography)); [LI[CartesianIndex(i, j)] for i in 12:18 for j in 15:20])

_filltimes(seq) = (d = Dict{Int,Float64}();
    for e in seq[2:end], u in e.filled
        (u isa IncrementalUpdate) && u.value && !haskey(d, u.index) && (d[u.index] = e.timestamp)
    end; d)
_monotone(seq) = all(seq[i].timestamp <= seq[i+1].timestamp + 1e-9 for i in 1:length(seq)-1)

@testset "NBS: retention delays the downstream pit" begin
    ts = _plane_with_pit()
    @test numtraps(ts) == 1
    fp = _upstream_footprint(ts)
    w  = [SWIM.WeatherEvent(0.0, 1.0e-3)]

    # a slow-draining puddle stores runoff and releases it gradually
    slow = SWIM.DynNBSPlacement(SWIM.puddle(60.0; kOUT = 0.005), fp, CartesianIndex{2}[])
    seq0 = fill_sequence(ts, w)
    seqN = fill_sequence(ts, w; nbs = [slow])

    @test _monotone(seqN)
    @test haskey(_filltimes(seq0), 1) && haskey(_filltimes(seqN), 1)
    @test _filltimes(seqN)[1] > _filltimes(seq0)[1] * 1.3   # markedly later with the NBS
end

@testset "NBS: free pass-through barely changes the pit" begin
    ts = _plane_with_pit()
    fp = _upstream_footprint(ts)
    w  = [SWIM.WeatherEvent(0.0, 1.0e-3)]

    # tiny storage, fast overflow: the NBS re-emits almost everything at once (X → 0),
    # so the downstream pit fills at essentially the no-NBS time.
    thru = SWIM.DynNBSPlacement(SWIM.puddle(1.0e-3; kOUT = 50.0), fp, CartesianIndex{2}[])
    f0 = _filltimes(fill_sequence(ts, w))[1]
    fT = _filltimes(fill_sequence(ts, w; nbs = [thru]))[1]

    @test isapprox(fT, f0; rtol = 0.05)
end
