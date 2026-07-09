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

# ---------------------------------------------------------------------------
# Stage B2: the NBS overlay element is carried through setup_network — its
# re-emit targets (terrain exit boundary + piped outlets) pull their downstream
# structure into the network like culvert endpoints, and the DynNBS is assigned to
# the component that owns those targets.  See agent/NBS_OPTION1_OVERLAY_PLAN.md §7a.
# ---------------------------------------------------------------------------
@testset "NBS overlay element in setup_network" begin
    N    = 12
    grid = Float64[j for i in 1:N, j in 1:N]     # east-rising plane, flow due west
    ts   = spillanalysis(grid)
    LI   = LinearIndices(size(grid))
    CI   = CartesianIndices(size(grid))
    full = findall(fill(false, numtraps(ts)))

    foot = Int[LI[i, 6] for i in 1:N]
    nb   = DynNBS(1, foot, 1, CartesianIndex{2}[])
    seed = CI[LI[3, 9]]

    comps = setup_network(ts, [seed], full; nbs = [nb])
    @test sum(length(c.nbs) for c in comps) == 1          # the element is placed exactly once
    host = only(filter(c -> !isempty(c.nbs), comps))
    @test host.nbs[1].placement_ix == 1
    @test !isempty(host.flow_paths)                       # exit boundary pulled in its downstream

    # the exit boundary (column 5, west of the footprint) is materialised as network cells
    exit_cells = Set(c for (c, _w) in SWIM._nbs_exit_weights(ts, foot))
    hostcells  = Set(LI[c] for p in host.flow_paths for c in p.cells)
    @test any(in(hostcells), exit_cells)

    # no NBS -> no element, network unchanged in kind
    comps0 = setup_network(ts, [seed], full)
    @test all(isempty(c.nbs) for c in comps0)
end

# ---------------------------------------------------------------------------
# Stage B2: the NBS layer-storage ODE in the rate function — the layer cascade
# (top layer fed by the static footprint capture, each lower layer by the one
# above's infiltration) and the weighted re-emit that fills the overflow-delivery
# slots.  See agent/NBS_OPTION1_OVERLAY_PLAN.md §3/§6.
# ---------------------------------------------------------------------------
@testset "NBS layer ODE + weighted re-emit (rate function)" begin
    N    = 12
    grid = Float64[j for i in 1:N, j in 1:N]
    ts   = spillanalysis(grid)
    LI   = LinearIndices(size(grid)); CI = CartesianIndices(size(grid))
    full = findall(fill(false, numtraps(ts)))
    foot = Int[LI[i, 6] for i in 1:N]; A = Float64(length(foot))
    hostof(nb) = only(filter(c -> !isempty(c.nbs),
                             setup_network(ts, [CI[LI[3, 9]]], full; nbs = [nb])))

    # --- single-layer puddle: dS/dt = inflow - overflow; delivery == overflow -------
    Smax, K = 10.0, 2.0
    pl  = NBSPlacement(puddle(Smax; kOUT = K, nOUT = 1.0), foot, 1, CartesianIndex{2}[])
    nb  = DynNBS(1, foot, 1, CartesianIndex{2}[])
    net = hostof(nb); nt = length(net.traps)
    I   = 5.0
    p   = SWIM._build_rate_params(ts, net, zeros(N, N), zeros(nt);
                                  nbs_placements = [pl], nbs_inflow = [I])
    @test p.nbsplan.nlayer_total == 1
    @test p.nbsplan.n_slots == N              # one delivery slot per exit cell

    dV = zeros(nt + 1)
    # below threshold: overflow ~ 0, the layer fills at the inflow rate
    Vlo = vcat(zeros(nt), [Smax * A / 1000 * 0.5])     # S_mm = Smax/2 < Smax
    SWIM.dynNetworkRateFunction!(dV, Vlo, p)
    @test dV[nt + 1] ≈ I atol = 1e-6

    # above threshold: dS/dt = inflow - overflow; the weighted delivery sums to overflow
    Vhi  = vcat(zeros(nt), [100.0])
    S_mm = 100.0 * 1000 / A
    qo   = SWIM.compute_outflow(K, 1.0, Smax, S_mm) * 1e-3
    SWIM.dynNetworkRateFunction!(dV, Vhi, p)
    @test dV[nt + 1] ≈ I - qo
    SWIM._nbs_fill_actual!(p.scratch.nbs_actual, Vhi, p, nt)
    @test sum(p.scratch.nbs_actual) ≈ qo               # all exits on-network, weights sum to 1
    @test I ≈ dV[nt + 1] + qo                          # mass: inflow = ΔS/dt + overflow

    # --- two-layer green roof (soil -> drainage, piped outlet): cascade + mass -------
    gr     = elhadiGreenRoof(2.0, 1.0, 3.0, 4.0, 1.0, 1.0, 0.0, 1.0)  # soil, drainage
    outlet = CI[LI[6, 5]]                              # on the exit boundary (network cell)
    pl2    = NBSPlacement(gr, foot, 0, [outlet])       # n_terrain=0: drainage is piped
    nb2    = DynNBS(1, foot, 0, [outlet])
    net2   = hostof(nb2); nt2 = length(net2.traps)
    I2     = 6.0
    p2     = SWIM._build_rate_params(ts, net2, zeros(N, N), zeros(nt2);
                                     nbs_placements = [pl2], nbs_inflow = [I2])
    @test p2.nbsplan.nlayer_total == 2

    V2  = vcat(zeros(nt2), [50.0, 30.0])               # soil, drainage storage (m^3)
    dV2 = zeros(nt2 + 2)
    SWIM.dynNetworkRateFunction!(dV2, V2, p2)
    L   = gr.layers
    S1, S2 = 50.0 * 1000 / A, 30.0 * 1000 / A
    qi1 = SWIM.compute_outflow(L[1].Kinf, L[1].ninf, L[1].Smin, S1) * 1e-3  # soil -> drainage
    qo2 = SWIM.compute_outflow(L[2].Kout, L[2].nout, L[2].Smax, S2) * 1e-3  # drainage -> outlet
    @test dV2[nt2 + 1] ≈ I2 - qi1                      # soil: inflow - infiltration down
    @test dV2[nt2 + 2] ≈ qi1 - qo2                     # drainage: fed by soil, drains via pipe
    SWIM._nbs_fill_actual!(p2.scratch.nbs_actual, V2, p2, nt2)
    @test sum(p2.scratch.nbs_actual) ≈ qo2             # piped outlet on-network
    @test I2 ≈ (dV2[nt2 + 1] + dV2[nt2 + 2]) + qo2     # mass (drainage Kinf=0, no ground loss)
end

# ---------------------------------------------------------------------------
# Stage B2: an end-to-end solveDynNetwork! over a network carrying an NBS overlay
# element — its layer storage is appended to the ODE state and evolved to steady
# state, where inflow balances outflow.  Exercises the state-vector extension and
# the NBS-aware steady-state callback (here on a trap-free NBS-only component).
# ---------------------------------------------------------------------------
@testset "NBS overlay element solved to steady state" begin
    N    = 12
    grid = Float64[j for i in 1:N, j in 1:N]
    ts   = spillanalysis(grid)
    LI   = LinearIndices(size(grid)); CI = CartesianIndices(size(grid))
    full = findall(fill(false, numtraps(ts)))
    foot = Int[LI[i, 6] for i in 1:N]; A = Float64(length(foot))

    Smax, K, I = 10.0, 2.0, 5.0
    pl  = NBSPlacement(puddle(Smax; kOUT = K, nOUT = 1.0), foot, 1, CartesianIndex{2}[])
    nb  = DynNBS(1, foot, 1, CartesianIndex{2}[])
    net = only(filter(c -> !isempty(c.nbs),
                      setup_network(ts, [CI[LI[3, 9]]], full; nbs = [nb])))
    nt  = length(net.traps)

    state = vcat(zeros(nt), [1.0])                      # layer starts nearly empty
    res   = solveDynNetwork!(state, ts, net, zeros(N, N), zeros(nt);
                             nbs_placements = [pl], nbs_inflow = [I])
    @test res.kind == :none                            # settled to steady state (no topology event)

    # steady state: dS/dt = 0  =>  overflow == inflow.  For a linear puddle
    # (K*(S_mm - Smax)) that pins S_mm = Smax + I*1000/K.
    S_mm = state[nt + 1] * 1000 / A
    qo   = SWIM.compute_outflow(K, 1.0, Smax, S_mm) * 1e-3
    @test qo ≈ I atol = 1e-2                            # outflow balances the captured inflow
    @test S_mm ≈ Smax + I * 1000 / K rtol = 1e-3        # settled at the analytic steady storage
end

# ---------------------------------------------------------------------------
# Stage B2: NBS→NBS coupling.  When element A's overflow re-emits onto element B's
# footprint, A and B must land in one network component and B must capture A's
# (variable, storage-driven) overflow as extra layer-1 inflow — with no mass leak.
# The coupling is storage-driven, so it is a one-way forcing read from a state
# snapshot: no algebraic loop, no topological ordering.  Plan §4/§7a, decision 3.
# ---------------------------------------------------------------------------
@testset "NBS→NBS coupling (variable-rate, mass-conserving)" begin
    N    = 12
    grid = Float64[j for i in 1:N, j in 1:N]           # east-rising plane, flow due west
    ts   = spillanalysis(grid)
    LI   = LinearIndices(size(grid)); CI = CartesianIndices(size(grid))
    full = findall(fill(false, numtraps(ts)))
    A    = Float64(N)

    # A at column 6 re-emits west onto B at column 5 (A's exit boundary == B's footprint).
    colA = Int[LI[i, 6] for i in 1:N]
    colB = Int[LI[i, 5] for i in 1:N]
    plA  = NBSPlacement(puddle(10.0; kOUT = 2.0), colA, 1, CartesianIndex{2}[])
    plB  = NBSPlacement(puddle(10.0; kOUT = 3.0), colB, 1, CartesianIndex{2}[])
    dA   = DynNBS(1, colA, 1, CartesianIndex{2}[])
    dB   = DynNBS(2, colB, 1, CartesianIndex{2}[])

    comps = setup_network(ts, [CI[LI[3, 9]]], full; nbs = [dA, dB])
    net   = only(filter(c -> !isempty(c.nbs), comps))
    @test length(net.nbs) == 2                          # A and B in ONE component (coupled)
    nt    = length(net.traps)

    IA, IB = 4.0, 5.0
    p    = SWIM._build_rate_params(ts, net, zeros(N, N), zeros(nt);
                                   nbs_placements = [plA, plB], nbs_inflow = [IA, IB])
    plan = p.nbsplan
    kA   = findfirst(==(1), plan.placement_ix)          # local index of A (placement 1)
    kB   = findfirst(==(2), plan.placement_ix)          # local index of B (placement 2)
    @test isempty(plan.nbs_into[kA])                    # A has no upstream feeder
    @test length(plan.nbs_into[kB]) == N                # B captures A's N exit slots

    V   = vcat(zeros(nt), zeros(2)); V[nt + plan.state_base[kA] + 1] = 80.0
    V[nt + plan.state_base[kB] + 1] = 60.0
    dV  = zeros(nt + 2)
    SWIM.dynNetworkRateFunction!(dV, V, p)
    SA  = V[nt + plan.state_base[kA] + 1]; SB = V[nt + plan.state_base[kB] + 1]
    qoA = SWIM.compute_outflow(2.0, 1.0, 10.0, SA * 1000 / A) * 1e-3
    qoB = SWIM.compute_outflow(3.0, 1.0, 10.0, SB * 1000 / A) * 1e-3
    dSA = dV[nt + plan.state_base[kA] + 1]; dSB = dV[nt + plan.state_base[kB] + 1]

    @test dSA ≈ IA - qoA                                # A: static in, overflow out
    @test dSB ≈ IB + qoA - qoB                          # B: static in + A's overflow, own overflow out
    @test IA + IB ≈ dSA + dSB + qoB                     # mass: nothing created or lost (no leak)
end

# ---------------------------------------------------------------------------
# Stage B2: end-to-end through fill_sequence.  The NBS placement is carried into
# the dynamic-network context (element built from the placement, layer storage
# persisted across events, footprint capture fed from rateinfo), so it measurably
# changes the fill sequence — a downstream trap fills later because the NBS retains
# part of its inflow.  Also checks the no-NBS path stays inert.
# ---------------------------------------------------------------------------
@testset "NBS through fill_sequence (context integration)" begin
    N    = 15
    grid = Float64[(i - 8)^2 + (j - 8)^2 for i in 1:N, j in 1:N]   # bowl -> one central trap
    ts   = spillanalysis(grid)
    LI   = LinearIndices(size(grid))
    weather = [WeatherEvent(0.0, 1.0)]

    # trap -> the time it first becomes filled, reconstructed from the event log
    filltimes(seq) = (d = Dict{Int,Float64}();
        for e in seq[2:end], u in e.filled
            u.value && !haskey(d, u.index) && (d[u.index] = e.timestamp)
        end; d)

    # a north-rim footprint that drains south into the bowl (re-emits toward the trap)
    foot = Int[LI[i, j] for i in 3:5 for j in 4:11]
    pl   = NBSPlacement(puddle(2000.0; kOUT = 2.0), foot, 1, CartesianIndex{2}[])

    ftNo  = filltimes(fill_sequence(ts, weather; dyn_traps = [1]))
    ftYes = filltimes(fill_sequence(ts, weather; dyn_traps = [1], nbs = [pl]))

    @test haskey(ftNo, 1) && haskey(ftYes, 1)          # the trap fills in both runs
    @test ftYes[1] > ftNo[1] + 1.0                     # the NBS retains inflow -> fills later

    # an empty NBS vector leaves the run byte-identical (the NBS path stays inert)
    @test [e.timestamp for e in fill_sequence(ts, weather; dyn_traps = [1], nbs = NBSPlacement[])] ==
          [e.timestamp for e in fill_sequence(ts, weather; dyn_traps = [1])]
end
