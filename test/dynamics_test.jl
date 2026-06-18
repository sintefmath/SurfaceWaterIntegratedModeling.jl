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
        ok &= all(1 .<= p.merges .<= np)
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

    # connection purely via a merge
    paths = [DynFlowPath([c(1,1)], 0),
             DynFlowPath([c(2,2)], 0, Int[], Int[], [1])]  # path2 merges into path1
    comps = SWIM._path_components(paths, DynTrap[])
    @test length(comps) == 1
    @test Set(comps[1]) == Set([1, 2])
end

@testset "_resolve_cell_overlaps!: truncation and merge" begin
    # path2 shares cell (1,2) with path1 -> truncated there, becomes a tributary
    paths = [DynFlowPath([c(1,1), c(1,2), c(1,3)], 1),
             DynFlowPath([c(2,1), c(1,2), c(1,3)], 2)]
    SWIM._resolve_cell_overlaps!(paths, DynCulvert[])

    @test paths[1].merges == [2]            # path2 registered as tributary of path1
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
