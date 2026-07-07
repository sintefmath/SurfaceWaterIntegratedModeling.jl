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

# ---------------------------------------------------------------------------
# Build a network on the mini grid with a single-layer puddle NBS whose outlet is
# backfilled to its natural downstream cell.  Returns (tstruct, net, nt, nbs_local).
function _mini_nbs_net(; Smax = 5.0, kOUT = 1.0)
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    sz = size(grid); li = LinearIndices(sz)
    r = sz[1] ÷ 2; c = sz[2] ÷ 2
    foot = [li[r, c], li[r, c + 1], li[r + 1, c], li[r + 1, c + 1]]
    p = SWIM.NBSPlacement(SWIM.puddle(Smax; kOUT = kOUT, nOUT = 1.0), foot, [CartesianIndex(0, 0)])
    t = spillanalysis(grid; nbs = [p])
    nbs_objs = SWIM._nbs_elements(t)
    seeds = SWIM._dyn_seeds(t, Int[], SWIM.DynCulvert[], nbs_objs)
    net = only(filter(n -> !isempty(n.nbs), SWIM.setup_network(t, seeds, Int[]; nbs = nbs_objs)))
    return t, net, size(grid)
end

# ---------------------------------------------------------------------------
@testset "NBS rate function conserves (charged layer, no inflow)" begin
    t, net, gsz = _mini_nbs_net(; Smax = 0.0, kOUT = 1.0)  # Smax=0 -> any storage overflows
    nt = length(net.traps)
    p = SWIM._build_rate_params(t, net, zeros(gsz), zeros(nt))
    nl = SWIM._nbs_state_count(p)
    @test nl == 1

    V = zeros(nt + nl)
    V[nt + 1] = 4.0                    # charge the single layer (A=4 -> S=1000 mm)
    dV = similar(V)
    SWIM.dynNetworkRateFunction!(dV, V, p, 0.0)

    # layer sheds Qout = 1.0 m^3/time; that exact amount reappears downstream; the
    # NBS trap node is frozen; nothing is created or destroyed.
    @test isapprox(dV[nt + 1], -1.0; atol = 1e-9)      # layer outflow
    @test isapprox(sum(dV), 0.0; atol = 1e-9)          # global balance (no infil, on-network delivery)
    @test dV[p.nbsplan.trap_local[1]] == 0.0           # NBS trap frozen
end

# ---------------------------------------------------------------------------
@testset "NBS solve conserves mass over a window" begin
    t, net, gsz = _mini_nbs_net(; Smax = 5.0, kOUT = 1.0)
    nt = length(net.traps)
    infil = zeros(gsz)
    p = SWIM._build_rate_params(t, net, infil, zeros(nt))
    nl = SWIM._nbs_state_count(p)
    nbs_local = p.nbsplan.trap_local[1]

    Q = 0.1
    inflow = zeros(nt); inflow[nbs_local] = Q       # constant inflow into the NBS only
    state = zeros(nt + nl)
    T = 1.0
    res = solveDynNetwork!(state, t, net, infil, inflow; tmax = T)

    @test res.kind == :none                          # window elapsed, no topology event
    @test state[nt + 1] > 0.0                         # the layer filled
    # infiltration is zero and (for this grid/window) nothing exits the domain, so all
    # injected water is accounted for in the final trap volumes + layer storage.
    @test isapprox(sum(state), Q * T; atol = 1e-9)
end

# ---------------------------------------------------------------------------
@testset "NBS state persistence round-trip" begin
    t, net, gsz = _mini_nbs_net(; Smax = 5.0)
    nt  = length(net.traps)
    pix = net.nbs[1].placement_ix

    # read: a populated store restores the layer block; an empty store gives zeros.
    @test SWIM._nbs_layer_block(net, t, Dict(pix => [2.5])) == [2.5]
    @test SWIM._nbs_layer_block(net, t, Dict{Int,Vector{Float64}}()) == [0.0]

    # write-back: a context's evolved layer state lands in the store, keyed by placement.
    gix   = Int[tr.trap_ix for tr in net.traps]
    state = vcat(zeros(nt), [7.0])
    ctx   = SWIM.DynNetworkContext(net, state, gix, Set{Int}(), CartesianIndex{2}[],
                                   0.0, Float64[], (; time = Inf, trap = 0, kind = :none))
    store = Dict{Int,Vector{Float64}}()
    SWIM._store_nbs_state!(store, ctx, t)
    @test store[pix] == [7.0]
end

# ---------------------------------------------------------------------------
@testset "NBS fill_sequence runs over two weather periods" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    sz = size(grid); li = LinearIndices(sz)
    r = sz[1] ÷ 2; c = sz[2] ÷ 2
    foot = [li[r, c], li[r, c + 1], li[r + 1, c], li[r + 1, c + 1]]
    p = SWIM.NBSPlacement(SWIM.puddle(5.0; kOUT = 1.0), foot, [CartesianIndex(0, 0)])
    t = spillanalysis(grid; nbs = [p])

    # rain then dry: the NBS fills in period 1 and drains its carried-over storage in
    # period 2, exercising the cross-period nbs_state persistence path.
    we = [WeatherEvent(0.0, 1.0), WeatherEvent(5.0, 0.0)]
    seq = fill_sequence(t, we)
    @test !isempty(seq)
    @test issorted([e.timestamp for e in seq])
end

# ---------------------------------------------------------------------------
@testset "NBS outlet routed as a culvert-style delivery" begin
    t, net, gsz = _mini_nbs_net(; Smax = 5.0)
    p = SWIM._build_rate_params(t, net, zeros(gsz), zeros(length(net.traps)))
    plan = p.nbsplan
    @test plan.n_slots == 1                       # one on-network outlet -> one delivery slot
    # the outlet is delivered at exactly one place: a path cell position (:nbsout event)
    # or straight into a trap — never head-injected, never both.
    npath = sum(length(e) for e in plan.nbs_path_events)
    ntrap = sum(length(o) for o in plan.nbs_trap_outlets)
    @test npath + ntrap == 1
    @test plan.deliver_slot[1][1] == 1            # the single layer maps to slot 1
end

# ---------------------------------------------------------------------------
# Build a network on the mini grid for a given multi-layer NBSSystem, all outlets
# unspecified (backfilled).  Returns (tstruct, net, gridsize).
function _mini_nbs_net_system(system)
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    sz = size(grid); li = LinearIndices(sz)
    r = sz[1] ÷ 2; c = sz[2] ÷ 2
    foot = [li[r, c], li[r, c + 1], li[r + 1, c], li[r + 1, c + 1]]
    outlets = fill(CartesianIndex(0, 0), length(system.layers))
    p = SWIM.NBSPlacement(system, foot, outlets)
    t = spillanalysis(grid; nbs = [p])
    nbs_objs = SWIM._nbs_elements(t)
    seeds = SWIM._dyn_seeds(t, Int[], SWIM.DynCulvert[], nbs_objs)
    net = only(filter(n -> !isempty(n.nbs), SWIM.setup_network(t, seeds, Int[]; nbs = nbs_objs)))
    return t, net, sz
end

@testset "NBS multi-layer (mantilla, 3 layers) conserves mass" begin
    # surface -> soil -> drainage; EVC=0 (ET deferred/ignored), bottom Kinf=0 (no ground loss)
    sys = SWIM.mantillaRRmodel(5.0, 5.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0)
    t, net, sz = _mini_nbs_net_system(sys)
    nt = length(net.traps)
    infil = zeros(sz)
    p = SWIM._build_rate_params(t, net, infil, zeros(nt))
    nl = SWIM._nbs_state_count(p)
    @test nl == 3                                 # three layer states appended
    @test p.nbsplan.n_slots == 3                  # three outlet-delivery slots

    Q = 0.1; nbs_local = p.nbsplan.trap_local[1]
    inflow = zeros(nt); inflow[nbs_local] = Q
    state = zeros(nt + nl); T = 1.0
    res = solveDynNetwork!(state, t, net, infil, inflow; tmax = T)

    @test isapprox(sum(state), Q * T; atol = 1e-8)     # nothing created/destroyed
    @test any(state[(nt + 1):(nt + nl)] .> 0)          # the cascade filled some layer
end

@testset "NBS multi-layer (elhadi green roof, 2 layers) conserves mass" begin
    sys = SWIM.elhadiGreenRoof(2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0)
    t, net, sz = _mini_nbs_net_system(sys)
    nt = length(net.traps)
    infil = zeros(sz)
    p = SWIM._build_rate_params(t, net, infil, zeros(nt))
    nl = SWIM._nbs_state_count(p)
    @test nl == 2

    Q = 0.1; nbs_local = p.nbsplan.trap_local[1]
    inflow = zeros(nt); inflow[nbs_local] = Q
    state = zeros(nt + nl); T = 1.0
    solveDynNetwork!(state, t, net, infil, inflow; tmax = T)
    @test isapprox(sum(state), Q * T; atol = 1e-8)
end

@testset "NBS multi-layer fill_sequence runs over two periods" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    sz = size(grid); li = LinearIndices(sz)
    r = sz[1] ÷ 2; c = sz[2] ÷ 2
    foot = [li[r, c], li[r, c + 1], li[r + 1, c], li[r + 1, c + 1]]
    sys = SWIM.mantillaRRmodel(5.0, 5.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0)
    p = SWIM.NBSPlacement(sys, foot, fill(CartesianIndex(0, 0), 3))
    t = spillanalysis(grid; nbs = [p])
    seq = fill_sequence(t, [WeatherEvent(0.0, 1.0), WeatherEvent(5.0, 0.0)])
    @test !isempty(seq)
    @test issorted([e.timestamp for e in seq])
end

# ---------------------------------------------------------------------------
@testset "NBS and culvert coexist in one build" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    sz = size(grid); li = LinearIndices(sz); ci = CartesianIndices(sz)
    r = sz[1] ÷ 2; c = sz[2] ÷ 2
    foot = [li[r, c], li[r, c + 1], li[r + 1, c], li[r + 1, c + 1]]
    p = SWIM.NBSPlacement(SWIM.puddle(5.0; kOUT = 1.0), foot, [CartesianIndex(0, 0)])
    t = spillanalysis(grid; nbs = [p])
    nbs_objs = SWIM._nbs_elements(t)
    nbs_trap = t.regions[foot[1]]

    # a downhill culvert clear of the NBS region (its footprint is the global minimum,
    # so a culvert into it would loop — the acyclicity guard would rightly reject that).
    idxs = findall(i -> t.regions[i] != nbs_trap, eachindex(t.regions))
    hi = ci[idxs[argmax(t.topography[idxs])]]
    lo = ci[idxs[argmin(t.topography[idxs])]]
    cv = SWIM.DynCulvert(hi, lo, 0.5, 0.6, 0.5, 0.02, 1.7)

    seeds = SWIM._dyn_seeds(t, Int[], [cv], nbs_objs)
    comps = SWIM.setup_network(t, seeds, Int[]; culverts = [cv], nbs = nbs_objs)
    @test sum(length(n.culverts) for n in comps) == 1     # culvert included
    @test sum(length(n.nbs) for n in comps) == 1          # NBS included
    # every component solves without error (mixed culvert + NBS build)
    for n in comps
        nt = length(n.traps)
        pp = SWIM._build_rate_params(t, n, zeros(sz), zeros(nt))
        st = zeros(nt + SWIM._nbs_state_count(pp))
        @test (solveDynNetwork!(st, t, n, zeros(sz), zeros(nt); tmax = 0.5); true)
    end
end
