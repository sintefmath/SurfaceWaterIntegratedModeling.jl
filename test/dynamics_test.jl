using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

c(i, j) = CartesianIndex(i, j)

# every internal reference (trap/path index) is 0 or a valid local index
function valid_network(net)
    np, nt = length(net.flow_paths), length(net.traps)
    ok = true
    for p in net.flow_paths
        ok &= 0 <= p.target_trap <= nt
        ok &= all(1 <= m <= np for (m, _) in p.merges)
    end
    for t in net.traps
        ok &= 0 <= t.spill_path <= np
    end
    return ok
end

# no grid cell is shared by two flow paths across the given networks
function disjoint_cells(nets)
    seen = Set{CartesianIndex{2}}()
    for net in nets, p in net.flow_paths, cell in p.cells
        cell in seen && return false
        push!(seen, cell)
    end
    return true
end

# ---------------------------------------------------------------------------
# Tier 2: functions that only touch a couple of tstruct fields (mockable)
# ---------------------------------------------------------------------------

@testset "_build_network" begin
    mock = (topography = zeros(3, 3),)  # only .topography is used (for CartesianIndices)

    # normal: a path flowing into a trap that spills out of the domain
    net = SWIM._build_network([[1, 2, 3]], [100], false, mock)
    @test length(net.flow_paths) == 1
    @test length(net.traps) == 1
    @test net.flow_paths[1].target_trap == 1          # path -> trap
    @test net.traps[1].spill_path == 0                # trap exits domain
    @test net.traps[1].trap_ix == 100
    @test net.flow_paths[1].cells == [c(1,1), c(2,1), c(3,1)]
    @test eltype(net.flow_paths[1].cells) == CartesianIndex{2}

    # start point inside the first full trap: chain begins with a trap that spills
    # into the first path (regression test for the inverted-wiring bug)
    net = SWIM._build_network([[1, 2, 3]], [100, 200], true, mock)
    @test length(net.flow_paths) == 1
    @test length(net.traps) == 2
    @test net.traps[1].spill_path == 1                # trap A spills into the path
    @test net.flow_paths[1].target_trap == 2          # path flows into trap B
    @test net.traps[2].spill_path == 0                # trap B exits domain

    # start inside a trap that spills straight out of the domain: lone trap, no path
    net = SWIM._build_network(Vector{Int}[], [100], true, mock)
    @test isempty(net.flow_paths)
    @test length(net.traps) == 1
    @test net.traps[1].spill_path == 0
end

@testset "_unfilled_trap_at" begin
    # cells 1,2 -> region 1 (supertraps 3,5,7); cell 3 -> region 2; cell 4 -> out
    mock = (regions = [1 2; 1 -1], supertraps_of = [[3, 5, 7], [4]])

    @test SWIM._unfilled_trap_at(mock, 1, [5, 7]) == 3      # uppermost unfilled
    @test SWIM._unfilled_trap_at(mock, 1, [3, 5, 7]) == 0   # all full
    @test SWIM._unfilled_trap_at(mock, 4, Int[]) == 0       # drains out of domain
end

# ---------------------------------------------------------------------------
# Tier 1: pure structural functions (no tstruct at all)
# ---------------------------------------------------------------------------

@testset "_combine_networks: index remapping" begin
    nA = DynNetwork([DynFlowPath([c(1,1)], 1)], [DynTrap(100, 0)], DynCulvert[])
    nB = DynNetwork([DynFlowPath([c(5,5)], 1)], [DynTrap(200, 1)], DynCulvert[])

    paths, traps, _ = SWIM._combine_networks([nA, nB])

    @test length(paths) == 2 && length(traps) == 2
    # network A's references keep their values (zero offset)
    @test paths[1].target_trap == 1
    @test traps[1].spill_path == 0           # 0 stays 0
    # network B's references are shifted by the per-type offsets
    @test paths[2].target_trap == 2
    @test traps[2].spill_path == 2
end

@testset "_path_components" begin
    # path1 -> trap1 -> path2 (connected); path3 standalone
    paths = [DynFlowPath([c(1,1)], 1),
             DynFlowPath([c(2,2)], 0),
             DynFlowPath([c(3,3)], 0)]
    traps = [DynTrap(100, 2)]               # trap1 spills into path2
    comps = SWIM._path_components(paths, traps)
    @test length(comps) == 2
    @test Set(map(Set, comps)) == Set([Set([1, 2]), Set([3])])

    # connection purely via a merge (junction position = 1 since path1 has 1 cell)
    paths = [DynFlowPath([c(1,1)], 0, Int[], Int[], [(2, 1)]),   # path1 is main; trib=path2 at pos 1
             DynFlowPath([c(2,2)], 0)]
    comps = SWIM._path_components(paths, DynTrap[])
    @test length(comps) == 1
    @test Set(comps[1]) == Set([1, 2])
end

@testset "_resolve_cell_overlaps!: truncation and merge" begin
    # path2 shares cell (1,2) with path1 -> truncated there, becomes a tributary
    paths = [DynFlowPath([c(1,1), c(1,2), c(1,3)], 1),
             DynFlowPath([c(2,1), c(1,2), c(1,3)], 2)]
    SWIM._resolve_cell_overlaps!(paths, DynCulvert[])

    @test paths[1].merges == [(2, 2)]        # path2 is trib of path1; junction at cell index 2 (c(1,2))
    @test paths[2].cells == [c(2,1)]        # truncated before the shared cell
    @test paths[2].target_trap == 0         # truncated path no longer targets a trap
    @test disjoint_cells([DynNetwork(paths, DynTrap[], DynCulvert[])])
end

@testset "_topological_order" begin
    # chain p1 -> t1 -> p2, fed in shuffled order
    paths = [DynFlowPath([c(1,1)], 1),   # global 1 = p1
             DynFlowPath([c(2,2)], 0)]   # global 2 = p2
    traps = [DynTrap(100, 2)]            # global 1 = t1, spills into p2
    sp, st = SWIM._topological_order([2, 1], [1], paths, traps)
    @test sp == [1, 2]                   # upstream (p1) before downstream (p2)
    @test st == [1]
end

@testset "_merge_networks: disjoint networks pass through" begin
    nA = DynNetwork([DynFlowPath([c(1,1), c(1,2)], 1)], [DynTrap(100, 0)], DynCulvert[])
    nB = DynNetwork([DynFlowPath([c(5,5), c(5,6)], 1)], [DynTrap(200, 0)], DynCulvert[])
    out = SWIM._merge_networks([nA, nB])

    @test length(out) == 2
    @test all(valid_network, out)
    @test disjoint_cells(out)
    @test Set(t.trap_ix for net in out for t in net.traps) == Set([100, 200])
end

@testset "_merge_networks: overlapping networks collapse" begin
    # both paths share (2,2),(3,3) and (deterministically) reach the same trap 100
    nA = DynNetwork([DynFlowPath([c(1,1), c(2,2), c(3,3)], 1)], [DynTrap(100, 0)], DynCulvert[])
    nB = DynNetwork([DynFlowPath([c(9,9), c(2,2), c(3,3)], 1)], [DynTrap(100, 0)], DynCulvert[])
    out = SWIM._merge_networks([nA, nB])

    @test length(out) == 1
    net = out[1]
    @test valid_network(net)
    @test disjoint_cells(out)
    @test length(net.flow_paths) == 2       # primary + truncated tributary
    @test length(net.traps) == 1            # duplicate trap removed
    @test net.traps[1].trap_ix == 100
    # exactly one tributary relationship, pointing at a valid path
    merges = reduce(vcat, [p.merges for p in net.flow_paths])
    @test length(merges) == 1

    @test isempty(SWIM._merge_networks(DynNetwork[]))
end

# ---------------------------------------------------------------------------
# Tier 3: integration tests against a real TrapStructure (mini.txt artifact)
# ---------------------------------------------------------------------------

# number of traps that no flow path leads into; such a "source" trap sits at the
# top of the chain, which happens only when the start point is inside that trap
source_traps(net) = count(ti -> all(p -> p.target_trap != ti, net.flow_paths),
                          1:length(net.traps))

@testset "setup_network on mini.txt" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    ts = spillanalysis(grid, usediags=true)
    allfull = collect(1:numtraps(ts))   # valid (downward-closed) fill: all traps full

    # a high start whose chain threads many traps -> a single nontrivial network
    long = setup_network(ts, [CartesianIndex(7, 119)], allfull)
    @test length(long) == 1
    @test valid_network(long[1])
    @test disjoint_cells(long)
    @test length(long[1].traps) > 20            # long spill path (58 traps in practice)

    # start inside a full trap -> chain begins with a source trap (no path feeds it)
    @test source_traps(long[1]) == 1
    # start on a slope -> every trap is fed by a path, so there is no source trap
    slope = setup_network(ts, [CartesianIndex(1, 1)], allfull)
    @test source_traps(slope[1]) == 0

    # two starts whose paths converge -> a single merged network
    merged = setup_network(ts, [CartesianIndex(193, 8), CartesianIndex(190, 9)], allfull)
    @test length(merged) == 1
    @test valid_network(merged[1])
    @test disjoint_cells(merged)

    # two starts that drain independently -> two disjoint networks
    separate = setup_network(ts, [CartesianIndex(195, 7), CartesianIndex(179, 37)], allfull)
    @test length(separate) == 2
    @test all(valid_network, separate)
    @test disjoint_cells(separate)
end

# ---------------------------------------------------------------------------
# networksolver.jl
# ---------------------------------------------------------------------------

@testset "networksolver: flow routing (pure topology)" begin
    T(ix, sp) = DynTrap(ix, sp)
    rf = SWIM._route_flow

    # path_cell_infil is per-cell; for single-cell paths it is a 1-element vector.
    # chain  A(full) -> path -> B(leaf)
    net = DynNetwork([DynFlowPath([c(1,1)], 2)], [T(10,1), T(20,0)], DynCulvert[])
    @test rf(net, [10.0,0.0], Bool[true,false], [0.0,0.0], [[0.0]]) ≈ [10.0,10.0]   # no loss
    @test rf(net, [10.0,0.0], Bool[true,false], [0.0,0.0], [[3.0]]) ≈ [10.0,7.0]    # path loss
    @test rf(net, [10.0,0.0], Bool[true,false], [4.0,0.0], [[0.0]]) ≈ [10.0,6.0]    # footprint loss
    @test rf(net, [10.0,0.0], Bool[false,false], [0.0,0.0], [[0.0]]) ≈ [10.0,0.0]   # leaf doesn't spill
    @test rf(net, [5.0,0.0],  Bool[true,false], [0.0,0.0], [[8.0]]) ≈ [5.0,0.0]     # loss floored at 0

    # path_inflow: direct inflow onto a leaf path reaches the downstream trap (after path loss)
    @test rf(net, [0.0,0.0], Bool[false,false], [0.0,0.0], [[0.0]];
             path_inflow=[5.0]) ≈ [0.0, 5.0]                                         # no path loss
    @test rf(net, [0.0,0.0], Bool[false,false], [0.0,0.0], [[3.0]];
             path_inflow=[5.0]) ≈ [0.0, 2.0]                                         # path loss applied
    @test rf(net, [0.0,0.0], Bool[false,false], [0.0,0.0], [[8.0]];
             path_inflow=[5.0]) ≈ [0.0, 0.0]                                         # loss exceeds inflow

    # tributary merge: path1(A)->trap2(C,leaf), path2(B) merges into A at pos 1 (only cell)
    net2 = DynNetwork([DynFlowPath([c(1,1)], 2, Int[], Int[], [(2, 1)]),
                       DynFlowPath([c(5,5)], 0)],
                      [T(10,1), T(20,0), T(30,2)], DynCulvert[])
    # B: max(20-2,0)=18 joins A at pos 1 (no pre-junc infil); A: max((10+18)-1,0)=27 -> C
    @test rf(net2, [10.0,0.0,20.0], Bool[true,false,true], [0.0,0.0,0.0], [[1.0],[2.0]]) ≈ [10.0,27.0,20.0]
    # path_inflow onto tributary path2: flows through path2 then path1 infiltration
    @test rf(net2, [0.0,0.0,0.0], Bool[false,false,false], [0.0,0.0,0.0], [[1.0],[2.0]];
             path_inflow=[0.0,7.0]) ≈ [0.0, max(7.0-2.0-1.0, 0.0), 0.0]             # 7-2(p2)-1(p1)=4

    # merge-fix: trib joins a 2-cell main path at junction 2 (not at head)
    # main path cells: [c(1,1), c(1,2)], infil=[1.0, 3.0]; trib at junction 2 → post-infil=3.0
    # head=0.5 (path_inflow[1]), trib delivers 5.0 (path_inflow[2], trib infil=0.0)
    # correct:  max(0.5-1.0,0)=0  + 5.0 = 5.0;  max(5.0-3.0,0)=2.0
    # old approx would give: max(0.5+5.0-4.0,0)=1.5  (wrong)
    net_mf = DynNetwork([DynFlowPath([c(1,1), c(1,2)], 2, Int[], Int[], [(2, 2)]),
                         DynFlowPath([c(5,5)], 0)],
                        [T(10,1), T(20,0)], DynCulvert[])
    @test rf(net_mf, [0.0,0.0], Bool[false,false], [0.0,0.0], [[1.0,3.0],[0.0]];
             path_inflow=[0.5, 5.0]) ≈ [0.0, 2.0]

    # spill exits the domain
    net3 = DynNetwork([DynFlowPath([c(1,1)], 0)], [T(10,1)], DynCulvert[])
    @test rf(net3, [10.0], Bool[true], [0.0], [[2.0]]) ≈ [10.0]
end

@testset "networksolver: geometry / rate / solve on mini.txt" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    ts   = spillanalysis(grid, usediags=true)
    net  = setup_network(ts, [CartesianIndex(7, 119)], collect(1:numtraps(ts)))[1]
    nt   = length(net.traps)

    @testset "geometry helpers (wetted-area infiltration)" begin
        infil = fill(0.5, size(ts.topography))
        geom  = SWIM._build_trap_geometry(ts, net, infil)
        @test length(geom) == nt
        # zvt caching: pre-computed tables must give identical geometry
        zvt_cached = SWIM._compute_z_vol_tables(ts)
        geom2 = SWIM._build_trap_geometry(ts, net, infil; zvt=zvt_cached)
        @test length(geom2) == nt
        for (g, g2) in zip(geom, geom2)
            @test g.capacity == g2.capacity && g.zmin == g2.zmin
        end
        for g in geom
            g.capacity <= 0 && continue
            nfp = length(g.footprint)
            # empty: level at bottom; full: level Inf and whole footprint wetted
            @test SWIM.water_level(g, 0.0) == g.zmin
            @test isinf(SWIM.water_level(g, g.capacity))
            @test SWIM.wetted_infiltration(g, g.capacity) ≈ 0.5 * nfp
            # half-full: level within the trap, infiltration monotone in fill
            Vh = g.capacity / 2
            @test g.zmin <= SWIM.water_level(g, Vh) <= ts.spillpoints[g.trap_ix].elevation + 1e-9
            @test SWIM.wetted_infiltration(g, 0.0) <=
                  SWIM.wetted_infiltration(g, Vh)  <=
                  SWIM.wetted_infiltration(g, g.capacity) + 1e-9
            # table-based level agrees with the Roots-based waterheight
            @test isapprox(SWIM.water_level(g, Vh), SWIM.waterheight(ts, g.trap_ix, Vh);
                           rtol=1e-3, atol=1e-3)
        end
    end

    @testset "rate function dV/dt" begin
        # zero infiltration: full traps are steady pass-throughs, leaf accumulates
        pA = SWIM._build_rate_params(ts, net, zeros(size(ts.topography)), fill(1.0, nt))
        caps = [g.capacity for g in pA.geom]
        VA = copy(caps); VA[end] = 0.0
        dV = similar(VA); dynNetworkRateFunction!(dV, VA, pA, 0.0)
        spilling = [VA[i] >= caps[i] for i in 1:nt]
        @test all(abs.(dV[findall(spilling)]) .< 1e-9)   # full traps dV ~ 0
        @test dV[end] > 0                                 # leaf fills

        # starved (no inflow, only infiltration): every full trap drains
        pC = SWIM._build_rate_params(ts, net, fill(0.1, size(ts.topography)), zeros(nt))
        VC = [g.capacity for g in pC.geom]
        dVC = similar(VC); dynNetworkRateFunction!(dVC, VC, pC, 0.0)
        @test all(dVC .< 0)
        @test dVC ≈ -pC.footprint_infil
    end

    @testset "solveDynNetwork events" begin
        caps = [g.capacity for g in SWIM._build_trap_geometry(ts, net, zeros(size(ts.topography)))]
        leaf = nt

        # FILL: zero infiltration, uniform inflow, leaf empty -> fills and spills
        V0 = copy(caps); V0[leaf] = 0.0
        res = solveDynNetwork(ts, net, zeros(size(ts.topography)), fill(1.0, nt), V0)
        @test res.kind == :fill
        @test res.trap == net.traps[leaf].trap_ix
        @test isfinite(res.time) && res.time > 0
        @test isapprox(res.state[leaf], caps[leaf]; rtol=1e-3)
        @test length(res.state) == nt

        # zvt caching: pre-computed tables give the same fill result
        zvt_cached = SWIM._compute_z_vol_tables(ts)
        res_zvt = solveDynNetwork(ts, net, zeros(size(ts.topography)), fill(1.0, nt), V0;
                                  zvt=zvt_cached)
        @test res_zvt.kind == res.kind && res_zvt.trap == res.trap
        @test isapprox(res_zvt.time, res.time; rtol=1e-6)

        # path_inflow: leaf fed only via path inflow (no trap inflow), still fills
        V0_pi = copy(caps); V0_pi[leaf] = 0.0
        path_qi = zeros(length(net.flow_paths))
        # find the spill path of the last full trap, which feeds into the leaf
        spill_p = net.traps[leaf-1].spill_path
        spill_p > 0 && (path_qi[spill_p] = 1.0)
        res_pi = solveDynNetwork(ts, net, zeros(size(ts.topography)), zeros(nt), V0_pi;
                                 path_inflow=path_qi)
        @test res_pi.kind == :fill && res_pi.trap == net.traps[leaf].trap_ix
        @test isfinite(res_pi.time) && res_pi.time > 0

        # STAGNATION -> Inf: only the leaf is fed, infiltration only on its footprint,
        # so it reaches a steady level below capacity (no further event)
        infil = zeros(size(ts.topography))
        infil[ts.footprints[net.traps[leaf].trap_ix]] .= 0.3
        inflow = zeros(nt); inflow[leaf] = 0.5
        res2 = solveDynNetwork(ts, net, infil, inflow, V0)
        @test res2.kind == :none
        @test res2.trap == 0
        @test res2.time == Inf
        @test res2.state[leaf] < caps[leaf]

        # UNSPILL: all traps full, infiltration exceeds inflow -> immediately draining
        VC = copy(caps)
        res_us = solveDynNetwork(ts, net, fill(0.1, size(ts.topography)), zeros(nt), VC)
        @test res_us.kind == :unspill
        @test res_us.time == 0.0
        @test res_us.trap > 0
        @test all(res_us.state .<= caps)               # no state above capacity
        @test any(res_us.state .< caps)                # at least one trap below full

        # NO EVOLVING: every trap already full, zero infiltration -> steady state
        res3 = solveDynNetwork(ts, net, zeros(size(ts.topography)), zeros(nt), copy(caps))
        @test res3.time == Inf && res3.trap == 0 && res3.kind == :none
    end
end
