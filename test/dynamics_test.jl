using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

c(i, j) = CartesianIndex(i, j)

# a DynCulvert with arbitrary (unused-for-topology) hydraulic parameters
cvlt(inlet, outlet) = DynCulvert(inlet, outlet, 0.5, 0.6, 0.5, 0.02, 1.7)

# every internal reference (trap/path/culvert index) is 0 or a valid local index
function valid_network(net)
    np, nt, nc = length(net.flow_paths), length(net.traps), length(net.culverts)
    ok = true
    for p in net.flow_paths
        ok &= 0 <= p.target_trap <= nt
        ok &= all(1 <= m <= np for (m, _) in p.merges)
        ok &= all(1 <= ci <= nc for ci in p.culvert_inlets)
        ok &= all(1 <= ci <= nc for ci in p.culvert_outlets)
    end
    for t in net.traps
        ok &= 0 <= t.spill_path <= np
        ok &= all(1 <= ci <= nc for ci in t.culvert_inlets)
        ok &= all(1 <= ci <= nc for ci in t.culvert_outlets)
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

@testset "_combine_subnets: index remapping" begin
    nA = DynNetwork([DynFlowPath([c(1,1)], 1)], [DynTrap(100, 0)], DynCulvert[])
    nB = DynNetwork([DynFlowPath([c(5,5)], 1)], [DynTrap(200, 1)], DynCulvert[])

    paths, traps = SWIM._combine_subnets([nA, nB])

    @test length(paths) == 2 && length(traps) == 2
    # network A's references keep their values (zero offset)
    @test paths[1].target_trap == 1
    @test traps[1].spill_path == 0           # 0 stays 0
    # network B's references are shifted by the per-type offsets
    @test paths[2].target_trap == 2
    @test traps[2].spill_path == 2
end

@testset "_components" begin
    nocv = (Tuple{Symbol,Int}[], Tuple{Symbol,Int}[])
    pathsets(comps) = Set(Set(p) for (p, _) in comps)

    # path1 -> trap1 -> path2 (connected); path3 standalone
    paths = [DynFlowPath([c(1,1)], 1),
             DynFlowPath([c(2,2)], 0),
             DynFlowPath([c(3,3)], 0)]
    traps = [DynTrap(100, 2)]               # trap1 spills into path2
    comps = SWIM._components(paths, traps, nocv...)
    @test length(comps) == 2
    @test pathsets(comps) == Set([Set([1, 2]), Set([3])])
    # trap1 lives in the {1,2} component
    comp12 = comps[findfirst(((p, _),) -> Set(p) == Set([1, 2]), comps)]
    @test comp12[2] == [1]

    # connection purely via a merge (junction position = 1 since path1 has 1 cell)
    paths = [DynFlowPath([c(1,1)], 0, Int[], Int[], [(2, 1)]),   # path1 is main; trib=path2 at pos 1
             DynFlowPath([c(2,2)], 0)]
    comps = SWIM._components(paths, DynTrap[], nocv...)
    @test length(comps) == 1
    @test Set(comps[1][1]) == Set([1, 2])

    # a culvert links two otherwise-disjoint single-path components into one
    paths = [DynFlowPath([c(1,1)], 0), DynFlowPath([c(9,9)], 0)]
    inlet_owner  = [(:path, 1)]
    outlet_owner = [(:path, 2)]
    comps = SWIM._components(paths, DynTrap[], inlet_owner, outlet_owner)
    @test length(comps) == 1
    @test Set(comps[1][1]) == Set([1, 2])

    # a lone trap (no path) survives as its own component
    comps = SWIM._components(DynFlowPath[], [DynTrap(100, 0)], nocv...)
    @test length(comps) == 1
    @test comps[1] == (Int[], [1])
end

@testset "_culvert_owners: trap footprint wins over flow path" begin
    # 3x3 grid; trap 1's footprint is the single cell (1,1) (linear index 1).  A path
    # runs through (1,1) and (2,2), so (1,1) is claimed by both -> the trap must win.
    mock  = (topography = zeros(3, 3), footprints = [[1]])
    paths = [DynFlowPath([c(1,1), c(2,2)], 0)]
    traps = [DynTrap(1, 0)]
    cvs   = [cvlt(c(1,1), c(2,2)),    # inlet on the shared cell; outlet path-only
             cvlt(c(3,3), c(1,1))]    # inlet off-network; outlet on the shared cell
    io, oo = SWIM._culvert_owners(mock, paths, traps, cvs)
    @test io[1] == (:trap, 1)        # shared cell resolves to the trap, not the path
    @test oo[1] == (:path, 1)
    @test io[2] == (:none, 0)        # (3,3) belongs to nothing
    @test oo[2] == (:trap, 1)
    # no culverts -> empty owner vectors (tstruct is allowed to be `nothing`)
    @test SWIM._culvert_owners(nothing, paths, traps, DynCulvert[]) ==
          (Tuple{Symbol,Int}[], Tuple{Symbol,Int}[])
end

@testset "_resolve_cell_overlaps!: truncation and merge" begin
    # path2 shares cell (1,2) with path1 -> truncated there, becomes a tributary
    paths = [DynFlowPath([c(1,1), c(1,2), c(1,3)], 1),
             DynFlowPath([c(2,1), c(1,2), c(1,3)], 2)]
    SWIM._resolve_cell_overlaps!(paths)

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

@testset "setup_network culverts on mini.txt" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    ts = spillanalysis(grid, usediags=true)
    allfull = collect(1:numtraps(ts))

    # Inclusion / assignment: a culvert between two footprint cells of two in-network
    # traps is included and registered on the inlet trap and the outlet trap.
    # On the :long network (start (7,119)): inlet cell (7,119) -> trap 233 (the start
    # trap), outlet cell (199,4) -> trap 13 (far downstream).
    @testset "inclusion and assignment (trap <-> trap)" begin
        inlet, outlet = CartesianIndex(7, 119), CartesianIndex(199, 4)
        out = setup_network(ts, [CartesianIndex(7, 119)], allfull; culverts=[cvlt(inlet, outlet)])
        @test length(out) == 1
        net = out[1]
        @test valid_network(net)
        @test disjoint_cells(out)
        @test length(net.culverts) == 1
        @test net.culverts[1].inlet == inlet && net.culverts[1].outlet == outlet
        ti_in  = findfirst(t -> t.trap_ix == 233, net.traps)
        ti_out = findfirst(t -> t.trap_ix == 13,  net.traps)
        @test ti_in !== nothing && ti_out !== nothing
        @test net.traps[ti_in].culvert_inlets   == [1]
        @test net.traps[ti_out].culvert_outlets == [1]
    end

    # Cross-network merge: two starts give two disjoint networks; a culvert from a
    # trap in one to a trap in the other merges them into a single network.
    @testset "cross-network merge" begin
        separate = setup_network(ts, [CartesianIndex(195, 7), CartesianIndex(179, 37)], allfull)
        @test length(separate) == 2          # precondition: disjoint without the culvert
        # trap 71 (footprint (179,37)) is in one net, trap 22 (footprint (196,6)) the other
        out = setup_network(ts, [CartesianIndex(195, 7), CartesianIndex(179, 37)], allfull;
                            culverts=[cvlt(CartesianIndex(179, 37), CartesianIndex(196, 6))])
        @test length(out) == 1
        net = out[1]
        @test valid_network(net)
        @test disjoint_cells(out)
        @test length(net.culverts) == 1
        @test 71 in (t.trap_ix for t in net.traps)
        @test 22 in (t.trap_ix for t in net.traps)
    end

    # Non-inclusion: a culvert whose inlet and outlet are both off-network must not be
    # added; the network is unchanged and has no culverts.
    @testset "non-inclusion (both endpoints off-network)" begin
        base = setup_network(ts, [CartesianIndex(179, 37)], allfull)
        out  = setup_network(ts, [CartesianIndex(179, 37)], allfull;
                            culverts=[cvlt(CartesianIndex(1, 1), CartesianIndex(10, 10))])
        @test length(out) == length(base) == 1
        @test isempty(out[1].culverts)
        @test length(out[1].traps)      == length(base[1].traps)
        @test length(out[1].flow_paths) == length(base[1].flow_paths)
    end

    # Outlet in terrain -> downstream expansion: an included culvert whose outlet lands
    # on a bare-terrain cell traces a fresh downstream chain that joins the network,
    # growing the trap/path counts.  Here inlet (179,37) is in the network and outlet
    # (8,119) is a slope cell of the unrelated long chain.
    @testset "terrain outlet expands the network" begin
        base = setup_network(ts, [CartesianIndex(179, 37)], allfull)[1]
        out  = setup_network(ts, [CartesianIndex(179, 37)], allfull;
                            culverts=[cvlt(CartesianIndex(179, 37), CartesianIndex(8, 119))])
        @test length(out) == 1
        net = out[1]
        @test valid_network(net)
        @test disjoint_cells(out)
        @test length(net.culverts) == 1
        @test length(net.traps)      > length(base.traps)
        @test length(net.flow_paths) > length(base.flow_paths)
    end

    # Inlet/outlet on a flow path (not a trap): the culvert is registered on the
    # owning flow paths' culvert_inlets / culvert_outlets, and on no trap.
    @testset "inlet/outlet on flow paths" begin
        inlet, outlet = CartesianIndex(182, 34), CartesianIndex(190, 31)
        out = setup_network(ts, [CartesianIndex(179, 37)], allfull; culverts=[cvlt(inlet, outlet)])
        @test length(out) == 1
        net = out[1]
        @test valid_network(net)
        @test length(net.culverts) == 1
        pin  = findfirst(p -> 1 in p.culvert_inlets,  net.flow_paths)
        pout = findfirst(p -> 1 in p.culvert_outlets, net.flow_paths)
        @test pin  !== nothing && inlet  in net.flow_paths[pin].cells
        @test pout !== nothing && outlet in net.flow_paths[pout].cells
        @test all(isempty(t.culvert_inlets) && isempty(t.culvert_outlets) for t in net.traps)
    end

    # Several culverts in one network: each is included and registered exactly once,
    # with the local culvert indices forming a clean 1:n set.  On the :long network,
    # cv1 (233 -> 13) and cv2 (444 -> 233) connect three in-network traps.
    @testset "multiple culverts in one network" begin
        out = setup_network(ts, [CartesianIndex(7, 119)], allfull;
                            culverts=[cvlt(CartesianIndex(7, 119), CartesianIndex(199, 4)),
                                      cvlt(CartesianIndex(115, 68), CartesianIndex(7, 119))])
        @test length(out) == 1
        net = out[1]
        @test valid_network(net)
        @test length(net.culverts) == 2
        # every culvert is registered once on the inlet side and once on the outlet side
        in_regs  = reduce(vcat, [t.culvert_inlets for t in net.traps];
                          init=Int[]) ∪ reduce(vcat, [p.culvert_inlets for p in net.flow_paths]; init=Int[])
        out_regs = reduce(vcat, [t.culvert_outlets for t in net.traps];
                          init=Int[]) ∪ reduce(vcat, [p.culvert_outlets for p in net.flow_paths]; init=Int[])
        @test sort(in_regs)  == [1, 2]
        @test sort(out_regs) == [1, 2]
    end

    # Fix-point inclusion: a culvert is pulled in only because an *earlier* culvert's
    # expansion grew the network to touch it.  cvA (179,37 -> terrain 8,119) drags in
    # the long chain; cvB (115,68 -> 199,4) has *both* endpoints on that chain, so it
    # only becomes reachable after cvA's expansion (it is excluded on its own).
    @testset "fix-point inclusion (chained culverts)" begin
        cvA = cvlt(CartesianIndex(179, 37), CartesianIndex(8, 119))
        cvB = cvlt(CartesianIndex(115, 68), CartesianIndex(199, 4))
        # cvB alone touches nothing in the (179,37) network -> not included
        only_b = setup_network(ts, [CartesianIndex(179, 37)], allfull; culverts=[cvB])
        @test sum(length(n.culverts) for n in only_b) == 0
        # with cvA present, cvA's expansion brings cvB into reach -> both included
        out = setup_network(ts, [CartesianIndex(179, 37)], allfull; culverts=[cvA, cvB])
        @test length(out) == 1
        @test valid_network(out[1])
        @test length(out[1].culverts) == 2
    end

    # Terrain inlet -> expansion (mirror of the terrain-outlet case): an included
    # culvert whose *inlet* is on bare terrain traces its host chain into the network.
    @testset "terrain inlet expands the network" begin
        base = setup_network(ts, [CartesianIndex(179, 37)], allfull)[1]
        out  = setup_network(ts, [CartesianIndex(179, 37)], allfull;
                            culverts=[cvlt(CartesianIndex(8, 119), CartesianIndex(179, 37))])
        @test length(out) == 1
        net = out[1]
        @test valid_network(net)
        @test disjoint_cells(out)
        @test length(net.culverts) == 1
        @test length(net.traps)      > length(base.traps)
        @test length(net.flow_paths) > length(base.flow_paths)
    end
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

# ---------------------------------------------------------------------------
# Culvert hydraulics: constructor + culvert_rate
# ---------------------------------------------------------------------------

# culvert_rate / the convenience constructor only read `tstruct.topography`, so a
# bare NamedTuple stands in for a full TrapStructure here.
mock_ts(M) = (; topography = M)

@testset "DynCulvert convenience constructor" begin
    # flat terrain, inlet c(1,1) -> outlet c(1,2): horiz = 1 m, drop = 0 -> L = 1 m.
    ts = mock_ts(zeros(2, 2))
    cv = DynCulvert(ts, c(1, 1), c(1, 2); r = 0.5)   # n=0.013, D=1, L=1, R=D/4=0.25
    @test cv.r  == 0.5
    @test cv.Cd == 0.6 && cv.Ke == 0.5 && cv.Cw == 1.7   # SI defaults
    @test cv.Kf ≈ 19.63 * 0.013^2 * 1.0 / (0.25)^(4/3) rtol = 1e-12
    # overriding a coefficient propagates
    @test DynCulvert(ts, c(1, 1), c(1, 2); r = 0.5, Cw = 2.0).Cw == 2.0
end

@testset "culvert_rate" begin
    A = pi * 0.5^2                       # r = 0.5 -> D = 1, A = pi/4
    # steep drop: inlet invert 10 m above outlet, so outlet control never binds
    steep = mock_ts([10.0 0.0])          # z[1,1]=10, z[1,2]=0
    cv = DynCulvert(c(1, 1), c(1, 2), 0.5, 0.6, 0.5, 1.0, 1.7)   # raw: Kf = 1.0

    # weir regime: inlet not submerged -> Q = Cw * D * H^1.5, inlet control governs
    qw = culvert_rate(cv, steep; inlet_submerged = false, inlet_head = 0.5,
                                 outlet_submerged = false, outlet_head = 0.0)
    @test qw ≈ 1.7 * 1.0 * 0.5^1.5 rtol = 1e-6

    # orifice regime: inlet submerged -> Q = Cd * A * sqrt(2 g H)
    qo = culvert_rate(cv, steep; inlet_submerged = true, inlet_head = 2.0,
                                 outlet_submerged = false, outlet_head = 0.0)
    @test qo ≈ 0.6 * A * sqrt(2 * 9.81 * 2.0) rtol = 1e-6

    # both submerged on flat terrain with small head difference -> outlet control
    # is the bottleneck and governs via the min().
    flat = mock_ts([0.0 0.0])
    qc = culvert_rate(cv, flat; inlet_submerged = true, inlet_head = 3.0,
                                outlet_submerged = true, outlet_head = 2.0)
    @test qc ≈ A * sqrt(2 * 9.81 * 1.0) / sqrt(1 + 0.5 + 1.0) rtol = 1e-6
    # and it is indeed the more restrictive of the two
    @test qc < 0.6 * A * sqrt(2 * 9.81 * 3.0)

    # free-outfall branch runs and stays bounded by inlet control
    qf = culvert_rate(cv, flat; inlet_submerged = true, inlet_head = 1.0,
                                outlet_submerged = false, outlet_head = 0.0)
    @test 0.0 <= qf <= 0.6 * A * sqrt(2 * 9.81 * 1.0) + 1e-9

    # zero head -> zero flow
    @test culvert_rate(cv, flat; inlet_submerged = false, inlet_head = 0.0,
                                 outlet_submerged = false, outlet_head = 0.0) == 0.0
end

@testset "culvert_rate reverse flow" begin
    A = pi * 0.5^2
    flat = mock_ts([0.0 0.0])
    cv = DynCulvert(c(1, 1), c(1, 2), 0.5, 0.6, 0.5, 1.0, 1.7)   # Kf = 1.0

    # outlet pool higher than inlet pool: drowned.  Default (downhill-only) -> 0.
    @test culvert_rate(cv, flat; inlet_submerged = true, inlet_head = 1.0,
                                 outlet_submerged = true, outlet_head = 3.0) == 0.0

    # same conditions with allow_reverse -> negative flow (outlet -> inlet),
    # governed here by outlet control on the reverse driving head dH = 3 - 1 = 2.
    qr = culvert_rate(cv, flat; inlet_submerged = true, inlet_head = 1.0,
                                outlet_submerged = true, outlet_head = 3.0,
                                allow_reverse = true)
    @test qr < 0
    @test qr ≈ -A * sqrt(2 * 9.81 * 2.0) / sqrt(1 + 0.5 + 1.0) rtol = 1e-6

    # symmetry: swapping the two heads flips the sign, same magnitude
    qf = culvert_rate(cv, flat; inlet_submerged = true, inlet_head = 3.0,
                                outlet_submerged = true, outlet_head = 1.0,
                                allow_reverse = true)
    @test qf ≈ -qr rtol = 1e-6

    # a forward-flow case is unchanged by allow_reverse
    @test culvert_rate(cv, flat; inlet_submerged = true, inlet_head = 3.0,
                                 outlet_submerged = true, outlet_head = 1.0) ≈ qf rtol = 1e-6
end
