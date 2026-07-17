# Dynamic-network membership / live-lifecycle tests (network_updating.jl):
# init_in_counts!, detach_spill! (+ compaction / re-root), grow_spill!, fusion, and the
# set-level entry points apply_fill! / apply_unfill! / apply_empty!.
#
using Test, SurfaceWaterIntegratedModeling, LazyArtifacts, Random

const DM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
ci(i, j) = CartesianIndex(i, j)

# a full flow path with explicit departure point
fp(cells, dep, tt; cin = Tuple{Int,Int}[], cout = Tuple{Int,Int}[],
                    nbs_in = Tuple{Int,Int}[], nbs = Tuple{Int,Int}[], mg = Tuple{Int,Int}[]) =
    DynFlowPath(cells, dep, tt, cin, cout, nbs_in, nbs, mg)

# in_count of every trap equals a from-scratch recompute
counts_consistent(net) =
    [t.in_count for t in net.traps] ==
    [t.in_count for t in DM.init_in_counts!(deepcopy(net)).traps]
counts_consistent(comps::Vector) = all(counts_consistent, comps)

all_trap_ix(comps) = sort([t.trap_ix for net in comps for t in net.traps])
partition(comps)   = Set(Set(t.trap_ix for t in net.traps) for net in comps if !isempty(net.traps))
no_dup(comps)      = (a = all_trap_ix(comps); length(a) == length(unique(a)))

# a mini.txt spillanalysis + a spread of dyn_coord seeds, shared across the terrain tests
const _mini_tstruct = spillanalysis(loadgrid(joinpath(artifact"swim_testdata", "data", "small",
                                                      "mini.txt")); usediags = true)
const _mini_seeds = let CIx = CartesianIndices(_mini_tstruct.topography),
                        N = length(_mini_tstruct.topography)
    [CIx[i] for i in 1:(N ÷ 80):N]
end
_mini_build(full) = DM.setup_network(_mini_tstruct, full; dyn_coords = _mini_seeds)

# ---------------------------------------------------------------------------
@testset "reachability primitives (hand-built nets)" begin
    # chain: seed -> t1 -> t2 -> t3 ; detaching t1 cascades away t2, t3
    paths = [fp([ci(1,1)], ci(1,1), 1), fp([ci(2,2)], ci(2,2), 2), fp([ci(3,3)], ci(3,3), 3)]
    traps = [DynTrap(101,2,Int[],Int[]), DynTrap(102,3,Int[],Int[]), DynTrap(103,0,Int[],Int[])]
    net = DynNetwork(paths, traps, DynCulvert[])
    DM.init_in_counts!(net)
    @test [t.in_count for t in net.traps] == [1, 1, 1]

    detached = detach_spill!(net, 1)
    @test detached == [102, 103]                       # returns stable trap_ix
    @test length(net.flow_paths) == 1                  # compacted
    @test length(net.traps) == 1 && net.traps[1].trap_ix == 101
    @test net.traps[1].spill_path == 0 && net.flow_paths[1].target_trap == 1

    # re-root: a live tributary keeps the target fed; survivor emanates from its source
    p2 = [fp([ci(5,5), ci(6,6)], ci(5,5), 2; mg = [(2,2)]), fp([ci(9,9)], ci(9,9), 0),
          fp([ci(1,1)], ci(1,1), 1), fp([ci(8,8)], ci(8,8), 3)]
    t2 = [DynTrap(201,1,Int[],Int[]), DynTrap(202,0,Int[],Int[]), DynTrap(203,2,Int[],Int[])]
    net2 = DynNetwork(p2, t2, DynCulvert[])
    DM.init_in_counts!(net2)
    @test detach_spill!(net2, 1) == Int[]              # trap1 unspills; tributary keeps t2 fed
    surv = net2.flow_paths[findfirst(p -> p.target_trap == 2, net2.flow_paths)]
    @test surv.cells == [ci(9,9), ci(6,6)] && surv.departure_point == ci(9,9)
    tC = net2.traps[findfirst(t -> t.trap_ix == 203, net2.traps)]
    @test tC.spill_path > 0 && net2.flow_paths[tC.spill_path] === surv   # no redirect needed
end

# ---------------------------------------------------------------------------
@testset "seed order is immaterial (no downstream sort)" begin
    full = Set(1:2:length(_mini_tstruct.spillpoints))
    ref  = partition(_mini_build(full))
    for ss in (reverse(_mini_seeds), shuffle(MersenneTwister(1), _mini_seeds),
               shuffle(MersenneTwister(2), _mini_seeds))
        comps = DM.setup_network(_mini_tstruct, full; dyn_coords = ss)
        @test partition(comps) == ref
    end
end

# ---------------------------------------------------------------------------
@testset "grow_spill! on mini.txt" begin
    grew = 0
    for lin in 1:length(_mini_tstruct.topography)
        comps = try
            DM.setup_network(_mini_tstruct, Int[]; dyn_coords = [CartesianIndices(_mini_tstruct.topography)[lin]])
        catch; continue; end
        for net in comps
            li = findfirst(t -> t.spill_path == 0, net.traps)
            li === nothing && continue
            sp = _mini_tstruct.spillpoints[net.traps[li].trap_ix]
            sp.downstream_region_cell == sp.current_region_cell && continue
            g = deepcopy(net)
            added = DM.grow_spill!(g, _mini_tstruct, Set([net.traps[li].trap_ix]), li)
            @test counts_consistent(g)
            @test no_dup([g])
            @test g.traps[li].spill_path != 0
            isempty(added) || (grew += 1)
            grew >= 20 && break
        end
        grew >= 20 && break
    end
    @test grew > 0
end

# ---------------------------------------------------------------------------
@testset "fusion regrow round-trip" begin
    full  = Set(1:2:length(_mini_tstruct.spillpoints))
    comps = _mini_build(full)
    a = findfirst(net -> length(net.traps) >= 2, comps)
    before = all_trap_ix(comps)
    DM._fuse_components!(comps, [a], _mini_tstruct, full)        # single-component regrow
    @test all_trap_ix(comps) == before                          # exact coverage preserved
    @test counts_consistent(comps)
end

# ---------------------------------------------------------------------------
@testset "composition harness (incremental vs fresh)" begin
    nt = length(_mini_tstruct.spillpoints)
    # incremental comps must have the same dynamic trap set as fresh and be a coarsening of it
    # (incremental never under-merges; fission is deferred so it may over-merge)
    function matches_fresh(comps, fresh)
        inc = Dict{Int,Int}()
        for (k, net) in enumerate(comps), t in net.traps; inc[t.trap_ix] = k; end
        Set(keys(inc)) == Set(t.trap_ix for net in fresh for t in net.traps) || return false
        for fc in fresh
            isempty(fc.traps) && continue
            length(unique(inc[t.trap_ix] for t in fc.traps)) == 1 || return false
        end
        counts_consistent(comps)
    end

    for run in 1:4
        rng = MersenneTwister(run)
        cur = Set(1:2:nt)
        comps = _mini_build(cur)
        for _ in 1:50
            fillc   = [t.trap_ix for net in comps for t in net.traps if t.spill_path == 0]
            unfillc = [t.trap_ix for net in comps for t in net.traps if t.spill_path != 0]
            dofill = isempty(unfillc) ? true : isempty(fillc) ? false : rand(rng, Bool)
            if dofill
                tix = rand(rng, fillc); push!(cur, tix); apply_fill!(comps, _mini_tstruct, cur, tix)
            else
                tix = rand(rng, unfillc); delete!(cur, tix); apply_unfill!(comps, _mini_tstruct, tix)
            end
            @test matches_fresh(comps, _mini_build(cur))
        end
    end
end

# ---------------------------------------------------------------------------
@testset "apply_empty! de-subsumption cycle" begin
    nt  = length(_mini_tstruct.spillpoints)
    cur = Set(1:2:nt)
    comps = _mini_build(cur)
    # find a not-full node whose fill completes its parent's subtraps (subsumption)
    A = 0
    for net in comps, t in net.traps
        if t.spill_path == 0 && DM._fill_subsumes(_mini_tstruct, t.trap_ix, union(cur, Set([t.trap_ix])))
            A = t.trap_ix; break
        end
    end
    @test A > 0
    P = parentof(_mini_tstruct, A)

    push!(cur, A); apply_fill!(comps, _mini_tstruct, cur, A)     # subsume {A, siblings} -> P
    @test DM._locate_trap(comps, A) === nothing                 # A no longer a node
    @test DM._locate_trap(comps, P) !== nothing                 # parent is the node

    delete!(cur, A); apply_empty!(comps, _mini_tstruct, cur, P)  # A drains below rim -> de-subsume
    @test DM._locate_trap(comps, A) !== nothing                 # A restored
    @test DM._locate_trap(comps, P) === nothing
    @test partition(comps) == partition(_mini_build(cur))       # matches a fresh build
end
