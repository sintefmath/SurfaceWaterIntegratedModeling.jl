using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# NBS as a signed-diff router element (agent/NBS_INVARIANTS.md).  `watercourses` is
# NBS-oblivious; the NBS captures its live input `Q_in` at the footprint-inflow cells and emits
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
    # A footprint-inflow cell intercepts ALL passing network flow into its live input Q_in (nbs_draw),
    # positive or negative — so an upstream NBS storing (a negative diff riding `current`) shows up as
    # a deficit in a downstream NBS's input, not as flow that bypasses it to the trap.
    D = SWIM._path_delivered!
    net = SWIM.DynNetwork()                          # unused (no culvert events)
    run3 = [10.0, 10.0, 10.0]                        # saturated cells: signed flow passes unattenuated
    ev   = [(2, :nbsin, 1)]                          # one :nbsin at cell 2

    nd = [0.0]                                       # positive release from upstream
    d  = D(SWIM._segment_maps(run3, ev), 5.0, ev, Float64[], Float64[], nothing, net, nothing, nd)
    @test nd[1] ≈ 5.0 && d ≈ 0.0                     # fully captured, nothing passes on

    nd = [0.0]                                       # negative: upstream NBS storing (deficit)
    d  = D(SWIM._segment_maps(run3, ev), -3.0, ev, Float64[], Float64[], nothing, net, nothing, nd)
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

# Build the network rate params for an NBS placement on `ts`, the way fill_sequence does.
function _nbs_rate_params(ts, nbs; rain = 1.0e-3)
    infil = zeros(size(ts.topography))                    # (footprints are forced to zero anyway)
    zvt   = SWIM._compute_z_vol_tables(ts)
    filled = Vector{Bool}(ts.trapvolumes .== 0.0)
    cur    = fill(SWIM.FilledAmount(0.0, 0.0), numtraps(ts))
    sgraph = SWIM.compute_complete_spillgraph(ts, filled)
    rate   = compute_flow(sgraph, rain, infil, ts)
    driver = SWIM.build_network_driver(ts,
                 SWIM._dyn_seeds(ts, CartesianIndex{2}[], SWIM.DynCulvert[]),
                 SWIM.DynCulvert[], findall(filled), cur, rate, infil, zvt, 0.0, Inf;
                 nbs_placements = nbs, precipitation = rain,
                 nbs_state = Dict{Int,Vector{Float64}}())
    ctx = only(c for c in driver.contexts if !isempty(c.net.nbs))
    return ctx.params
end

@testset "NBS: layer cascade conserves mass" begin
    ts  = _plane_with_pit()
    fp  = _upstream_footprint(ts)
    # a puddle that both overflows AND infiltrates to ground, so every mass term is live
    nbs = [SWIM.DynNBSPlacement(SWIM.puddle(1.0; kOUT = 0.5, kINF = 0.1), fp, CartesianIndex{2}[])]
    p   = _nbs_rate_params(ts, nbs)

    nt     = length(p.geom)
    nlayer = p.nbsplan.nlayer_total
    V  = zeros(nt + nlayer)
    for i in (nt + 1):(nt + nlayer); V[i] = 5.0; end      # positive storage → no V<=0 floor fires
    dV = similar(V)
    SWIM.dynNetworkRateFunction!(dV, V, p)

    el = p.nbsplan.elems[1]
    Q_in = el.O_0_total                                    # traps empty ⇒ no spills ⇒ nbs_draw == 0
    @test Q_in > 0.0                                       # the footprint really does receive flow
    O_surface = 0.0; ground = 0.0; sum_dV = 0.0
    for (l, L) in enumerate(el.system.layers)
        S_mm = V[nt + el.state_base + l] * 1000.0 / el.A
        O_surface += SWIM.compute_outflow(L.Kout, L.nout, L.Smax, S_mm) * el.A * 1e-3   # overflow leaves (mm→m³ via A)
        (l == length(el.system.layers)) &&                                             # only the bottom infiltrates to ground
            (ground = SWIM.compute_outflow(L.Kinf, L.ninf, L.Smin, S_mm) * el.A * 1e-3)
        sum_dV += dV[nt + el.state_base + l]
    end
    # storage rate = in − out − to-ground: no mass created or lost in the cascade
    @test isapprox(sum_dV, Q_in - O_surface - ground; atol = 1e-12, rtol = 1e-9)
end

@testset "NBS: a sink inside the footprint is rejected" begin
    # `internal_accumulation_cells` (the plan's ponding endpoints) is footprint ∩ trap_bottoms,
    # and a sink is deliberately not a trap bottom — so a sink inside a footprint would swallow
    # flow the placement never accounts for.  setup_network must refuse it.
    grid = Float64[j for i in 1:30, j in 1:30]
    for i in 12:18, j in 3:6; grid[i, j] = 1.0; end
    inside = CartesianIndex(15, 17)              # sits in _upstream_footprint's 12:18 × 15:20 block

    ts_ok = spillanalysis(grid, usediags = true)
    fp    = _upstream_footprint(ts_ok)
    @test LinearIndices(size(ts_ok.topography))[inside] in fp     # the sink really is in the footprint
    mk(t) = [SWIM.DynNBSPlacement(SWIM.puddle(1.0), fp, CartesianIndex{2}[])]

    # same terrain, same footprint: accepted with no sink, refused once the sink is declared
    @test setup_network(ts_ok, Int[]; nbs = mk(ts_ok)) isa Vector{SWIM.DynNetwork}
    ts_sink = spillanalysis(grid, usediags = true, sinks = [inside])
    @test_throws Exception setup_network(ts_sink, Int[]; nbs = mk(ts_sink))

    # a sink outside every footprint stays legal
    ts_out = spillanalysis(grid, usediags = true, sinks = [CartesianIndex(25, 25)])
    @test setup_network(ts_out, Int[]; nbs = mk(ts_out)) isa Vector{SWIM.DynNetwork}
end

@testset "NBS: a network commit is not altered by the event that triggers it" begin
    # An upstream pit U fills and starts spilling across a green roof whose downpipe drains off the
    # domain, into a downstream pit D in the roof's network.  The event loop updates the runoff for
    # U's spill *before* committing D's network over the interval leading up to it; that commit must
    # still use the rates of the interval (frozen on the context), not the post-event runoff, or the
    # roof's correction removes U's spill from D's inflow over the whole past interval (water lost,
    # D fills late).  Splitting the rain into a second, identical period just before U fills forces
    # the commit without that hazard, so D's fill time must agree.
    N = 40
    grid = Float64[j for i in 1:N, j in 1:N]
    for i in 12:28, j in 3:6;   grid[i, j] = -20.0; end     # downstream pit D (west), deep
    for i in 17:23, j in 30:33; grid[i, j] = 28.0;  end     # upstream pit U (east), shallow
    ts = spillanalysis(grid, usediags = true)
    @test numtraps(ts) == 2
    D, U = argmax(ts.trapvolumes), argmin(ts.trapvolumes)

    LI  = LinearIndices(size(ts.topography))
    fp  = [LI[i, j] for i in 10:30 for j in 15:20]          # across U's spill path into D
    nbs() = [SWIM.DynNBSPlacement(SWIM.elhadiGreenRoof(10.0, 5.0, 0.5, 0.5, 1.0, 1.0, 0.0, 1.0),
                                  fp, [CartesianIndex(5, 21)])]  # downpipe drains west, off domain
    rain = 1.0e-3

    f1 = _filltimes(fill_sequence(ts, [SWIM.WeatherEvent(0.0, rain)]; nbs = nbs()))
    @test f1[U] < f1[D]                                     # U fills first, spilling toward D
    f2 = _filltimes(fill_sequence(ts, [SWIM.WeatherEvent(0.0, rain),
                                       SWIM.WeatherEvent(f1[U] - 1.0e-3, rain)]; nbs = nbs()))
    @test isapprox(f1[U], f2[U]; rtol = 1.0e-9)
    @test isapprox(f1[D], f2[D]; rtol = 1.0e-4)
end

@testset "NBS: overlapping footprints are rejected" begin
    ts = _plane_with_pit()
    fp = _upstream_footprint(ts)
    mk(f) = SWIM.DynNBSPlacement(SWIM.puddle(1.0), f, CartesianIndex{2}[])
    @test_throws Exception setup_network(ts, Int[]; nbs = [mk(fp), mk(fp[1:3])])
    @test setup_network(ts, Int[]; nbs = [mk(fp[1:3]), mk(fp[4:end])]) isa Vector{SWIM.DynNetwork}
end

@testset "NBS: the inflow cell's own residual applies before capture" begin
    # What enters the footprint is what LEAVES the inflow cell (that is what O_0 counts there), so
    # a network flow change is charged by the inflow cell's residual before it is drawn; this keeps
    # Q_in = O_0 + draw >= 0 without clamping.
    D  = SWIM._path_delivered!
    nt = SWIM.DynNetwork()
    ev = [(2, :nbsin, 1)]                            # inflow cell at position 2
    nd = [0.0]
    D(SWIM._segment_maps([10.0, -3.0, 10.0], ev), 5.0, ev, Float64[], Float64[], nothing, nt, nothing, nd)
    @test nd[1] ≈ 2.0                                # 3 of the 5 fill the cell's spare capacity
    nd = [0.0]
    D(SWIM._segment_maps([10.0, 2.0, 10.0], ev), -3.0, ev, Float64[], Float64[], nothing, nt, nothing, nd)
    @test nd[1] ≈ -2.0                               # cannot remove more than the 2 leaving the cell
end

# A context whose network did not change keeps the structural part of its rate params at a touch
# and recomputes only the flow-dependent part (_with_flows); that must equal a full rebuild.
@testset "NBS: reused rate params equal a full rebuild" begin
    ts = _plane_with_pit()
    pl = SWIM.DynNBSPlacement(SWIM.puddle(60.0; kOUT = 0.005), _upstream_footprint(ts),
                              CartesianIndex{2}[])
    SWIM._prepare_nbs!([pl], ts)
    net = only(SWIM.setup_network(ts, Int[]; nbs = [pl]))
    infil = zeros(size(ts.topography))
    params(rain, inflow) = SWIM._build_rate_params(ts, net, infil, inflow;
        runoff = first(watercourses(ts, fill(false, numtraps(ts)); precipitation = rain)),
        precipitation = rain)
    nt = length(net.traps)
    old, new = params(1e-3, fill(0.1, nt)), params(3e-3, fill(0.4, nt))
    R = first(watercourses(ts, fill(false, numtraps(ts)); precipitation = 3e-3))
    reused = SWIM._with_flows(old, ts, infil, fill(0.4, nt); runoff = R, precipitation = 3e-3)
    for f in fieldnames(SWIM.DynNetworkRateParams)
        f in (:geom, :nbsplan) && continue
        @test getfield(reused, f) == getfield(new, f)
    end
    @test reused.geom === old.geom                     # the structural part is shared, not rebuilt
    for (a, b) in zip(reused.nbsplan.elems, new.nbsplan.elems), f in fieldnames(SWIM.NBSElement)
        @test getfield(a, f) == getfield(b, f)
    end
end
