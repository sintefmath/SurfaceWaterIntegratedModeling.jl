using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# These tests were written against the old positional `setup_network(ts, dyn_coords, full)`;
# the retired old path is gone (gate D2).  `mk_network` forwards to the new keyword-form
# `setup_network` (build_network.jl) so the solver/geometry/culvert coverage stays exercised.
mk_network(ts, coords, full; kw...) = setup_network(ts, full; dyn_coords = coords, kw...)

# A source flow path (no cells) for hand-built solver fixtures.  The primary DynFlowPath
# constructor needs an explicit departure_point since there is no `first(cells)` to take.
srcpath(target) = DynFlowPath(CartesianIndex{2}[], CartesianIndex(1, 1), target,
                              Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[],
                              Tuple{Int,Int}[])

# Tolerance for the network-ODE fill/drain times reproducing the analytic (plain) path.
# The dynamic solver runs at abstol=1e-6 (m^3, ~mL) / reltol=1e-4 — accuracy calibrated to
# the physical need (millilitres, milliseconds), not to machine precision.  The resulting
# fill-time drift vs the analytic path is a few ×1e-5 (tens of microseconds) on the coupled
# full-coverage cases, far under a millisecond.  This is the physical bound; a regression
# that shifts a fill time by more than this is genuine, not solver noise.
const PARITY_TOL = 1e-4

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
        ok &= all(1 <= ci <= nc && 1 <= pos <= length(p.cells) for (ci, pos) in p.culvert_inlets)
        ok &= all(1 <= ci <= nc && 1 <= pos <= length(p.cells) for (ci, pos) in p.culvert_outlets)
    end
    for t in net.traps
        ok &= -1 <= t.spill_path <= np      # -1 = full trap spilling out of domain
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
# Tier 3: integration tests against a real TrapStructure (mini.txt artifact)
# ---------------------------------------------------------------------------
@testset "out-of-domain full trap uses spill_path == -1 sentinel" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    ts   = spillanalysis(grid, usediags=true)

    # With every trap full, the chain from (7,119) ends on a full trap that spills
    # straight out of the domain: it carries the out-of-domain sentinel spill_path == -1,
    # while every upstream full trap has a real in-network spill path (> 0).
    net = mk_network(ts, [CartesianIndex(7, 119)], collect(1:numtraps(ts)))[1]
    @test net.traps[end].spill_path == -1
    @test all(t.spill_path > 0 for t in net.traps[1:end-1])
    @test valid_network(net)

    # solveDynNetwork! must ACCEPT a FULL trap carrying the -1 sentinel (it satisfies
    # the FULL contract spill_path != 0).  All traps full with infiltration > inflow
    # drains immediately -> :unspill via the t=0 fast path, not a validator throw.
    caps = [g.capacity for g in SWIM._build_trap_geometry(ts, net, zeros(size(ts.topography)))]
    res  = solveDynNetwork!(copy(caps), ts, net,
                            fill(0.1, size(ts.topography)), zeros(length(net.traps)))
    @test res.kind == :unspill
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
        out = mk_network(ts, [CartesianIndex(7, 119)], allfull; culverts=[cvlt(inlet, outlet)])
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
        separate = mk_network(ts, [CartesianIndex(195, 7), CartesianIndex(179, 37)], allfull)
        @test length(separate) == 2          # precondition: disjoint without the culvert
        # trap 71 (footprint (179,37)) is in one net, trap 22 (footprint (196,6)) the other
        out = mk_network(ts, [CartesianIndex(195, 7), CartesianIndex(179, 37)], allfull;
                            culverts=[cvlt(CartesianIndex(179, 37), CartesianIndex(196, 6))])
        @test length(out) == 1
        net = out[1]
        @test valid_network(net)
        @test disjoint_cells(out)
        @test length(net.culverts) == 1
        @test 71 in (t.trap_ix for t in net.traps)
        @test 22 in (t.trap_ix for t in net.traps)
    end

    # Outlet in terrain -> downstream expansion: an included culvert whose outlet lands
    # on a bare-terrain cell traces a fresh downstream chain that joins the network,
    # growing the trap/path counts.  Here inlet (179,37) is in the network and outlet
    # (8,119) is a slope cell of the unrelated long chain.
    @testset "terrain outlet expands the network" begin
        base = mk_network(ts, [CartesianIndex(179, 37)], allfull)[1]
        out  = mk_network(ts, [CartesianIndex(179, 37)], allfull;
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
        out = mk_network(ts, [CartesianIndex(179, 37)], allfull; culverts=[cvlt(inlet, outlet)])
        @test length(out) == 1
        net = out[1]
        @test valid_network(net)
        @test length(net.culverts) == 1
        # culvert 1 is registered on a path's inlet/outlet list, and the stored cell
        # position points to the actual inlet/outlet cell.
        function find_pos(listf)
            for p in net.flow_paths, (ci, pos) in listf(p)
                ci == 1 && return (p, pos)
            end
            return (nothing, 0)
        end
        pin_p,  ipos = find_pos(p -> p.culvert_inlets)
        pout_p, opos = find_pos(p -> p.culvert_outlets)
        @test pin_p  !== nothing && pin_p.cells[ipos]  == inlet
        @test pout_p !== nothing && pout_p.cells[opos] == outlet
        @test all(isempty(t.culvert_inlets) && isempty(t.culvert_outlets) for t in net.traps)
    end

    # Several culverts in one network: each is included and registered exactly once,
    # with the local culvert indices forming a clean 1:n set.  On the :long network,
    # cv1 (233 -> 13) and cv2 (233 -> 444) connect three in-network traps.  Both run
    # downhill from the start trap 233, so the network stays acyclic and orderable.
    @testset "multiple culverts in one network" begin
        out = mk_network(ts, [CartesianIndex(7, 119)], allfull;
                            culverts=[cvlt(CartesianIndex(7, 119), CartesianIndex(199, 4)),
                                      cvlt(CartesianIndex(7, 119), CartesianIndex(115, 68))])
        @test length(out) == 1
        net = out[1]
        @test valid_network(net)
        @test length(net.culverts) == 2
        # every culvert is registered once on the inlet side and once on the outlet
        # side (traps store bare indices, paths store (index, position) tuples)
        in_regs  = reduce(vcat, [t.culvert_inlets for t in net.traps]; init=Int[]) ∪
                   [ci for p in net.flow_paths for (ci, _) in p.culvert_inlets]
        out_regs = reduce(vcat, [t.culvert_outlets for t in net.traps]; init=Int[]) ∪
                   [ci for p in net.flow_paths for (ci, _) in p.culvert_outlets]
        @test sort(in_regs)  == [1, 2]
        @test sort(out_regs) == [1, 2]
    end

    # An uphill / reverse culvert (inlet downstream of its outlet) makes the network
    # cyclic, so it cannot be ordered upstream-to-downstream; construction rejects it.
    @testset "uphill culvert is rejected" begin
        # (115,68)=trap 444 lies downstream of (7,119)=trap 233 on the :long chain,
        # so a culvert 444 -> 233 runs against terrain flow.
        @test_throws Exception setup_network(
            ts, [CartesianIndex(7, 119)], allfull;
            culverts=[cvlt(CartesianIndex(115, 68), CartesianIndex(7, 119))])
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
    net2 = DynNetwork([DynFlowPath([c(1,1)], 2, Tuple{Int,Int}[], Tuple{Int,Int}[],
                                   Tuple{Int,Int}[], [(2, 1)]),
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
    net_mf = DynNetwork([DynFlowPath([c(1,1), c(1,2)], 2, Tuple{Int,Int}[], Tuple{Int,Int}[],
                                     Tuple{Int,Int}[], [(2, 2)]),
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
    # This testset models a chain of full traps feeding one evolving frontier leaf
    # (`net.traps[end]`).  Build with every trap full EXCEPT that leaf, so the leaf is
    # the terminal TRANSITORY node (spill_path == 0).  An all-full build would instead
    # mark the leaf as a full trap spilling out of the domain (spill_path == -1);
    # making it evolve afterwards would then violate the three-state contract.
    all_traps   = collect(1:numtraps(ts))
    leaf_trapix = mk_network(ts, [CartesianIndex(7, 119)], all_traps)[1].traps[end].trap_ix
    net  = mk_network(ts, [CartesianIndex(7, 119)], setdiff(all_traps, [leaf_trapix]))[1]
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
            # full: wetted infiltration = 0.5 over the PONDING cells only.  Cells at or
            # above the spillpoint never pond as part of the trap, so they carry no
            # infiltration (see `_ponding_mask`); the full-capacity loss is therefore
            # <= the whole-footprint value and equals 0.5 times the ponding-cell count.
            nland = count(<(ts.spillpoints[g.trap_ix].elevation), g.bottom)
            @test SWIM.wetted_infiltration(g, g.capacity) ≈ 0.5 * nland
            @test SWIM.wetted_infiltration(g, g.capacity) <= 0.5 * nfp + 1e-12
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

    @testset "solveDynNetwork! events" begin
        caps = [g.capacity for g in SWIM._build_trap_geometry(ts, net, zeros(size(ts.topography)))]
        leaf = nt

        # FILL: zero infiltration, uniform inflow, leaf empty -> fills and spills.
        # state is mutated in place; no res.state field.
        state = copy(caps); state[leaf] = 0.0
        res = solveDynNetwork!(state, ts, net, zeros(size(ts.topography)), fill(1.0, nt))
        @test res.kind == :fill
        @test res.trap == net.traps[leaf].trap_ix
        @test isfinite(res.time) && res.time > 0
        @test isapprox(state[leaf], caps[leaf]; rtol=1e-3)   # state updated in place
        @test length(state) == nt

        # zvt caching: pre-computed tables give the same fill result
        zvt_cached = SWIM._compute_z_vol_tables(ts)
        state_zvt = copy(caps); state_zvt[leaf] = 0.0
        res_zvt = solveDynNetwork!(state_zvt, ts, net, zeros(size(ts.topography)), fill(1.0, nt);
                                   zvt=zvt_cached)
        @test res_zvt.kind == res.kind && res_zvt.trap == res.trap
        @test isapprox(res_zvt.time, res.time; rtol=1e-6)

        # path_inflow: leaf fed only via path inflow (no trap inflow), still fills
        state_pi = copy(caps); state_pi[leaf] = 0.0
        path_qi  = zeros(length(net.flow_paths))
        spill_p  = net.traps[leaf-1].spill_path
        spill_p > 0 && (path_qi[spill_p] = 1.0)
        res_pi = solveDynNetwork!(state_pi, ts, net, zeros(size(ts.topography)), zeros(nt);
                                  path_inflow=path_qi)
        @test res_pi.kind == :fill && res_pi.trap == net.traps[leaf].trap_ix
        @test isfinite(res_pi.time) && res_pi.time > 0

        # STAGNATION -> Inf: single isolated trap (57) with 11 cells, inflow=0.3,
        # infil=0.1/cell → max-infil=1.1 > inflow, equilibrium at ~3 cells wet (V_eq
        # well below capacity=2.3126).  The steady-state DiscreteCallback must fire
        # and return :none with the state settled at the sub-capacity equilibrium.
        # (Trap 13, the :long-chain leaf, cannot be used here because its second cell's
        # bottom equals the spillpoint, giving no sub-capacity zone where infil > inflow.)
        net_s   = DynNetwork([srcpath(0)],
                             [DynTrap(57, 0)], DynCulvert[])
        infil_s = zeros(size(ts.topography))
        infil_s[ts.footprints[57]] .= 0.1
        state_s = [0.0]
        caps_s  = SWIM._build_trap_geometry(ts, net_s, infil_s)[1].capacity
        res2 = solveDynNetwork!(state_s, ts, net_s, infil_s, [0.3])
        @test res2.kind == :none
        @test res2.trap == 0
        @test res2.time == Inf
        @test 0.0 < state_s[1] < caps_s   # settled below capacity

        # UNSPILL: upstream FULL traps have infiltration exceeding inflow -> immediately
        # draining.  The leaf is left EMPTY (0.0); it does not participate in the check.
        # The t=0 fast path fires and state is NOT modified (caller's job to decrement).
        state_us = copy(caps); state_us[leaf] = 0.0
        res_us = solveDynNetwork!(state_us, ts, net, fill(0.1, size(ts.topography)), zeros(nt))
        @test res_us.kind == :unspill
        @test res_us.time == 0.0
        @test res_us.trap > 0
        # state unchanged for t=0 unspill: caller must set state[trap] = prevfloat(C)
        @test state_us == [i == leaf ? 0.0 : caps[i] for i in 1:nt]

        # STEADY STATE: leaf TRANSITORY at prevfloat(C), zero inflow/infiltration ->
        # evolving but rate == 0 everywhere -> :none.
        state_ss = copy(caps); state_ss[leaf] = prevfloat(caps[leaf])
        res3 = solveDynNetwork!(state_ss, ts, net, zeros(size(ts.topography)), zeros(nt))
        @test res3.time == Inf && res3.trap == 0 && res3.kind == :none

        # VALIDATOR: FULL trap with spill_path == 0 is rejected
        state_bad = copy(caps)   # leaf is FULL but spill_path == 0 in this network
        @test_throws ErrorException solveDynNetwork!(
            state_bad, ts, net, zeros(size(ts.topography)), zeros(nt))
    end

    @testset "solveDynNetwork!: steady-state callback waits for all traps to settle" begin
        # Two independent leaf traps (57: 11-cell footprint, 313: 10-cell footprint) in a
        # synthetic DynNetwork.  Both start EMPTY with positive inflow; both settle at
        # sub-capacity equilibrium (inflow = infiltration).  Trap 313 settles first.
        # The :steadystate DiscreteCallback fires only when max(|dV/dt|) < abstol
        # for BOTH traps, not as soon as one trap's rate crosses zero.
        net_two = DynNetwork(
            [srcpath(1), srcpath(2)],
            [DynTrap(57, 0), DynTrap(313, 0)],
            DynCulvert[])
        ifil_ms   = fill(0.1, size(ts.topography))   # 0.1/cell → max 1.1 and 1.0
        inflow_ms = [0.7, 0.5]
        state_ms  = [0.0, 0.0]

        res_ms = solveDynNetwork!(state_ms, ts, net_two, ifil_ms, inflow_ms)
        @test res_ms.kind == :none

        # Both traps must be at equilibrium: rate ≈ 0 at the returned state.
        # If only trap 313 had stagnated and integration terminated early,
        # trap 57 would still be filling (du[1] >> 0).
        p_ms  = SWIM._build_rate_params(ts, net_two, ifil_ms, inflow_ms)
        du_ms = zeros(2)
        dynNetworkRateFunction!(du_ms, state_ms, p_ms, 0.0)
        @test abs(du_ms[1]) < 1e-5
        @test abs(du_ms[2]) < 1e-5
        @test 0.0 < state_ms[1] < p_ms.geom[1].capacity
        @test 0.0 < state_ms[2] < p_ms.geom[2].capacity
    end

    @testset "solveDynNetwork!: :empty registered only for parent traps" begin
        # _event_conditions must add :fill for ALL evolving traps and :empty only for
        # those whose trap_ix > numregions (parent/merged traps).  Leaf traps at V=0
        # are at their physical floor — no topology changes there.
        # Steady state is detected by the single global :steadystate entry (trap==0).
        nreg = numregions(ts)
        p_ec = SWIM._build_rate_params(ts, net, zeros(size(ts.topography)), fill(1.0, nt))

        # Pick the first parent trap and the leaf as the evolving set (bypass validator).
        pi = findfirst(i -> net.traps[i].trap_ix > nreg, 1:nt)
        li = nt
        @assert pi !== nothing "test requires a parent trap in the :long chain"

        conds = SWIM._event_conditions(p_ec, [pi, li], nreg)

        @test  any(c -> c.trap == li && c.kind == :fill,  conds)
        @test !any(c -> c.trap == li && c.kind == :empty, conds)  # leaf: no :empty

        @test  any(c -> c.trap == pi && c.kind == :fill,  conds)
        @test  any(c -> c.trap == pi && c.kind == :empty, conds)  # parent: has :empty

        # no stagnation entries — steady-state is handled by a separate DiscreteCallback
        @test !any(c -> c.kind == :stagnation  || c.kind == :steadystate, conds)

        # :unspill registered for ALL non-evolving traps, including zero-initial-inflow ones.
        # LeftRootFind prevents degenerate rootfinding on zero-starting conditions so
        # culvert-driven networks where inflow later goes 0 → positive → negative are covered.
        @test count(c -> c.kind == :unspill, conds) == nt - 2
    end

    @testset "solveDynNetwork!: two-call cascade (caller rebuilds network after fill)" begin
        # Demonstrate the correct caller protocol for a fill event:
        #   1. Solve with leaf 233 → :fill event.
        #   2. Push 233 to full_traps, rebuild; next leaf is 220.
        #   3. Solve again → :fill for 220.
        start  = CartesianIndex(7, 119)
        ifil_c = zeros(size(ts.topography))

        ft_c = Int[]
        n1_c = mk_network(ts, [start], ft_c)[1]
        @test n1_c.traps[1].trap_ix == 233 && n1_c.traps[1].spill_path == 0

        g1_c  = SWIM._build_trap_geometry(ts, n1_c, ifil_c)
        C_233 = g1_c[1].capacity
        s1_c  = [0.0]
        res1_c = solveDynNetwork!(s1_c, ts, n1_c, ifil_c, [1.0])
        @test res1_c.kind == :fill && res1_c.trap == 233
        @test isfinite(res1_c.time) && res1_c.time > 0
        @test isapprox(s1_c[1], C_233; rtol=1e-3)

        push!(ft_c, 233)
        n2_c = mk_network(ts, [start], ft_c)[1]
        @test length(n2_c.traps) == 2
        @test n2_c.traps[2].trap_ix == 220 && n2_c.traps[2].spill_path == 0

        g2_c  = SWIM._build_trap_geometry(ts, n2_c, ifil_c)
        C_220 = g2_c[2].capacity
        s2_c  = [C_233, 0.0]
        res2_c = solveDynNetwork!(s2_c, ts, n2_c, ifil_c, [1.0, 1.0])
        @test res2_c.kind == :fill && res2_c.trap == 220
        @test isapprox(s2_c[2], C_220; rtol=1e-3)
    end

    @testset "solveDynNetwork!: prevfloat(C) handoff after :unspill" begin
        # After :unspill the caller must:
        #   (a) leave state unchanged (guaranteed by the t=0 fast path)
        #   (b) set state[trap] = prevfloat(C)
        #   (c) remove the trap from full_traps and rebuild the network
        # The second call must enter the ODE normally (t=0 fast path must NOT re-fire
        # since the trap is now TRANSITORY, not FULL).
        start  = CartesianIndex(7, 119)
        ifil_u = fill(0.1, size(ts.topography))  # infiltration > inflow → unspill

        # 2-trap network: trap 233 FULL, trap 220 EMPTY leaf; zero inflow triggers :unspill.
        ft_u = [233]
        n2_u = mk_network(ts, [start], ft_u)[1]
        g2_u = SWIM._build_trap_geometry(ts, n2_u, ifil_u)
        C_u  = g2_u[1].capacity

        s_us = [C_u, 0.0]
        res_us = solveDynNetwork!(s_us, ts, n2_u, ifil_u, zeros(2))
        @test res_us.kind == :unspill && res_us.time == 0.0
        @test s_us == [C_u, 0.0]            # state must be unchanged for t=0 :unspill

        # Caller sets prevfloat(C), removes trap 233 from full_traps, rebuilds.
        filter!(!=(233), ft_u)
        n1_u = mk_network(ts, [start], ft_u)[1]
        s_pf = [prevfloat(C_u)]             # trap 233 just below capacity (TRANSITORY)

        # Second call: positive inflow, zero infiltration → trap refills to capacity.
        # If the t=0 fast path re-fired (V == C check), the call would return :unspill
        # immediately.  The correct result is :fill because V = prevfloat(C) < C.
        res_pf = solveDynNetwork!(s_pf, ts, n1_u, zeros(size(ts.topography)), [1.0])
        @test res_pf.kind == :fill
        @test res_pf.trap == 233
        @test isapprox(s_pf[1], C_u; rtol=1e-3)
    end

    @testset "solveDynNetwork!: tmax bounds the integration" begin
        # Single leaf trap (233), empty, inflow 1.0, zero infiltration -> dV/dt = 1,
        # so it fills exactly at t_e = capacity.
        start = CartesianIndex(7, 119)
        net   = mk_network(ts, [start], Int[])[1]
        @test net.traps[1].trap_ix == 233
        infil = zeros(size(ts.topography))
        C     = SWIM._build_trap_geometry(ts, net, infil)[1].capacity

        # unbounded (tmax = Inf default): :fill at t_e ≈ C, state clamped to C
        s = [0.0]
        res = solveDynNetwork!(s, ts, net, infil, [1.0])
        @test res.kind == :fill && res.trap == 233
        t_e = res.time
        @test isapprox(t_e, C; rtol=1e-3) && isapprox(s[1], C; rtol=1e-3)

        # tmax short of the event: no event; return (tmax, 0, :none) with state at tmax.
        # With dV/dt = 1 the volume at tmax is exactly tmax.
        tm  = t_e / 2
        s2  = [0.0]
        res2 = solveDynNetwork!(s2, ts, net, infil, [1.0]; tmax = tm)
        @test res2.kind == :none && res2.trap == 0
        @test res2.time == tm
        @test isapprox(s2[1], tm; rtol=1e-3) && 0.0 < s2[1] < C

        # tmax past the event: the event still fires (state clamped to C).
        s3  = [0.0]
        res3 = solveDynNetwork!(s3, ts, net, infil, [1.0]; tmax = t_e * 1.5)
        @test res3.kind == :fill && res3.trap == 233
        @test isapprox(res3.time, t_e; rtol=1e-3) && isapprox(s3[1], C; rtol=1e-3)

        # tmax == 0: nothing advances, state unchanged, no event.
        s0  = [0.3]
        res0 = solveDynNetwork!(s0, ts, net, infil, [1.0]; tmax = 0.0)
        @test res0.kind == :none && res0.time == 0.0 && s0 == [0.3]

        # Steady state below a finite tmax must still report Inf (not tmax): trap 57
        # settles at a sub-capacity equilibrium (inflow 0.3 < max infiltration 1.1).
        net_s   = DynNetwork([srcpath(0)], [DynTrap(57, 0)], DynCulvert[])
        infil_s = zeros(size(ts.topography)); infil_s[ts.footprints[57]] .= 0.1
        s_s     = [0.0]
        res_s   = solveDynNetwork!(s_s, ts, net_s, infil_s, [0.3]; tmax = 1e6)
        @test res_s.kind == :none && res_s.time == Inf
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

@testset "culvert_rate: partial-depth tailwater & reverse" begin
    A = pi * 0.5^2
    flat = mock_ts([0.0 0.0])
    cv = DynCulvert(c(1, 1), c(1, 2), 0.5, 0.6, 0.5, 1.0, 1.7)   # D=1, Kf=1.0
    # outlet-control rate for a real tailwater giving dH = 0.10 m
    q_dH(dH) = A * sqrt(2 * 9.81 * dH) / sqrt(1 + 0.5 + 1.0)

    # NEITHER end submerged, but both hold water (weir regime both ways).  The
    # downstream pool -- not a free outfall -- sets the tailwater, so the driving
    # head is the real surface difference (0.95 - 0.85 = 0.10), and allow_reverse
    # must NOT change a forward-dominant result (it used to: Q_fwd - Q_rev bug).
    nr_off = culvert_rate(cv, flat; inlet_submerged = false, inlet_head = 0.95,
                                    outlet_submerged = false, outlet_head = 0.85)
    nr_on  = culvert_rate(cv, flat; inlet_submerged = false, inlet_head = 0.95,
                                    outlet_submerged = false, outlet_head = 0.85,
                                    allow_reverse = true)
    @test nr_off ≈ q_dH(0.10) rtol = 1e-6     # real tailwater, not free outfall
    @test nr_on  ≈ nr_off     rtol = 1e-6     # inlet higher -> no reverse triggered

    # outlet pool higher, neither submerged: downhill-only drowns (0), reverse
    # gives the genuine outlet->inlet flow as a negative, symmetric to nr_off.
    @test culvert_rate(cv, flat; inlet_submerged = false, inlet_head = 0.85,
                                 outlet_submerged = false, outlet_head = 0.95) == 0.0
    rev = culvert_rate(cv, flat; inlet_submerged = false, inlet_head = 0.85,
                                 outlet_submerged = false, outlet_head = 0.95,
                                 allow_reverse = true)
    @test rev ≈ -q_dH(0.10) rtol = 1e-6
    @test rev ≈ -nr_off     rtol = 1e-6
end

@testset "_route_flow: trap-to-trap culvert (mass conservation)" begin
    inlet, outlet = c(1, 1), c(1, 2)
    cv  = DynCulvert(inlet, outlet, 0.5, 0.6, 0.5, 1.0, 1.7)    # D = 1
    t1  = DynTrap(101, 0, [1], Int[])     # trap 1 hosts the culvert inlet (drained)
    t2  = DynTrap(102, 0, Int[], [1])     # trap 2 hosts the culvert outlet (filled)
    net = DynNetwork(DynFlowPath[], [t1, t2], [cv])
    ts  = mock_ts([0.0 0.0])              # both culvert inverts at elevation 0
    plan = SWIM._build_culvert_plan(net, ts)
    rf(levels) = SWIM._route_flow(net, [0.0, 0.0], [false, false], [0.0, 0.0],
                                  Vector{Float64}[]; cvplan = plan, trap_level = levels)

    # trap 1 higher -> culvert flows 1 -> 2; drawn at inlet == delivered at outlet
    Q = culvert_rate(cv, ts; inlet_submerged = true, inlet_head = 3.0,
                     outlet_submerged = false, outlet_head = 0.0)
    @test Q > 0
    inflow = rf([3.0, 0.0])
    @test inflow[1] ≈ -Q rtol = 1e-12              # trap 1 drained
    @test inflow[2] ≈  Q rtol = 1e-12              # trap 2 filled
    @test inflow[1] + inflow[2] ≈ 0.0 atol = 1e-12 # mass conserved (no external flux)

    # trap 2 higher -> uphill; reverse disallowed -> no flow
    @test rf([0.0, 3.0]) == [0.0, 0.0]
    # equal surfaces -> no driving head -> no flow
    @test rf([1.0, 1.0]) == [0.0, 0.0]
    # without a culvert plan the culvert is ignored entirely
    @test SWIM._route_flow(net, [5.0, 7.0], [false, false], [0.0, 0.0],
                           Vector{Float64}[]) == [5.0, 7.0]
end

@testset "_route_flow: flow-path culvert endpoint (capped by passing flow)" begin
    # path p1 (3 cells, exits domain); a culvert inlet sits at cell position 2 and
    # delivers into trap t1.  Path endpoints use head = D, not submerged.
    pcells = [c(1,1), c(1,2), c(1,3)]
    cv  = DynCulvert(pcells[2], c(5,5), 0.5, 0.6, 0.5, 1.0, 1.7)        # D = 1
    p1  = DynFlowPath(pcells, 0, [(1, 2)], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[])
    t1  = DynTrap(101, 0, Int[], [1])               # trap hosts the culvert outlet
    net = DynNetwork([p1], [t1], [cv])
    ts  = mock_ts(zeros(6, 6))                       # all inverts at 0
    plan = SWIM._build_culvert_plan(net, ts)
    cellinfil = [[0.0, 0.0, 0.0]]                    # no path infiltration

    # capacity at the path endpoint (inlet head = D, outlet dry)
    Q = culvert_rate(cv, ts; inlet_submerged = false, inlet_head = 1.0,
                     outlet_submerged = false, outlet_head = 0.0)
    @test Q > 0
    rf(F) = SWIM._route_flow(net, [0.0], [false], [0.0], cellinfil;
                             path_inflow = [F], cvplan = plan, trap_level = [0.0])

    # plenty of flow -> culvert abstracts its full capacity into the trap
    @test rf(10.0)[1] ≈ Q rtol = 1e-12
    # little flow -> abstraction is capped at the passing flow (mass-conserving)
    @test rf(0.3 * Q)[1] ≈ 0.3 * Q rtol = 1e-12
end

@testset "solveDynNetwork: culvert drain triggers :unspill" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    ts = spillanalysis(grid, usediags=true)
    allfull = collect(1:SWIM.numtraps(ts))
    net = mk_network(ts, [CartesianIndex(7, 119)], allfull;
                        culverts=[cvlt(CartesianIndex(7, 119), CartesianIndex(199, 4))])[1]
    nt = length(net.traps)
    ti_in = findfirst(t -> t.culvert_inlets == [1], net.traps)
    @test ti_in !== nothing

    # Upstream FULL traps + culvert draining one of them, leaf EMPTY.
    # The culvert makes the inlet trap the only one with a negative net inflow -> unspill.
    zer   = zeros(size(ts.topography))
    p     = SWIM._build_rate_params(ts, net, zer, zeros(nt))
    state = [g.capacity for g in p.geom]
    state[end] = 0.0   # leaf must be TRANSITORY/EMPTY per contract
    res = solveDynNetwork!(state, ts, net, zer, zeros(nt))
    @test res.kind == :unspill
    @test res.trap == net.traps[ti_in].trap_ix
end

@testset "fill_sequence dyn_traps parity (slice 3)" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    ts   = spillanalysis(grid, usediags=true)
    weather = [WeatherEvent(0.0, 1.0)]
    # trap -> time it first becomes filled, reconstructed from the event log
    filltimes(seq) = (d = Dict{Int,Float64}();
        for e in seq[2:end], u in e.filled
            u.value && !haskey(d, u.index) && (d[u.index] = e.timestamp)
        end; d)
    maxΔ(a, b) = maximum(abs(a[t] - b[t]) for t in keys(a))

    seqP = fill_sequence(ts, weather)
    ftP  = filltimes(seqP)

    # no dyn_traps must be byte-identical to the default (network path is inert)
    @test [e.timestamp for e in fill_sequence(ts, weather; dyn_traps=Int[])] ==
          [e.timestamp for e in seqP]

    # single isolated network reproduces plain fill_sequence essentially exactly
    ft1 = filltimes(fill_sequence(ts, weather; dyn_traps=[233]))
    @test Set(keys(ft1)) == Set(keys(ftP))
    @test maxΔ(ftP, ft1) < PARITY_TOL

    # FULL coverage: every trap solved as a dynamic network must reproduce plain to a
    # tight tolerance (ODE vs analytic; same event count, same set of filled traps)
    seqF = fill_sequence(ts, weather; dyn_traps=collect(1:numtraps(ts)))
    ftF  = filltimes(seqF)
    @test length(seqF) == length(seqP)
    @test Set(keys(ftF)) == Set(keys(ftP))
    @test maxΔ(ftP, ftF) < PARITY_TOL

    # MIXED coverage (subset networked) — same traps fill; timings agree to floating-point
    # precision.  Boundary traps newly absorbed by the network are projected to cur_time
    # before initialising their ODE state, so no accumulation-lag is introduced.
    ftM = filltimes(fill_sequence(ts, weather; dyn_traps=[233, 220]))
    @test Set(keys(ftM)) == Set(keys(ftP))
    @test maxΔ(ftP, ftM) < PARITY_TOL

    # Subtrap-hierarchy seed (trap 414 seeds a network that later absorbs sibling trap 18).
    ft414 = filltimes(fill_sequence(ts, weather; dyn_traps=[414]))
    @test Set(keys(ft414)) == Set(keys(ftP))
    @test maxΔ(ftP, ft414) < PARITY_TOL

    @testset "drought drain ordering and parity (staletime regression)" begin
        infil2 = fill(0.05, size(ts.topography))

        # Run rain until every trap is filled, then switch to drought.
        T_fill   = maximum(values(filltimes(
                       fill_sequence(ts, [WeatherEvent(0.0, 1.0)]; infiltration=infil2)))) + 1.0
        weather2 = [WeatherEvent(0.0, 1.0), WeatherEvent(T_fill, 0.0)]

        # Weather-period boundary events carry e.filled::Vector{Bool} (full snapshot),
        # not Vector{IncrementalUpdate{Bool}}.  Only incremental events record the
        # trap-by-trap transitions we care about here.
        inc_events(seq) = Iterators.filter(
            e -> e.filled isa Vector{IncrementalUpdate{Bool}}, seq)
        function draintimes(seq)
            d = Dict{Int,Float64}()
            for e in inc_events(seq), u in e.filled
                !u.value && !haskey(d, u.index) && (d[u.index] = e.timestamp)
            end; d
        end
        function filltimes_r(seq)   # filltimes robust to multi-weather sequences
            d = Dict{Int,Float64}()
            for e in inc_events(seq), u in e.filled
                u.value && !haskey(d, u.index) && (d[u.index] = e.timestamp)
            end; d
        end

        seqP = fill_sequence(ts, weather2; infiltration=infil2)
        dtP  = draintimes(seqP)

        # (1) Every trap that filled must eventually drain.
        # Regression: old min_net_inflow formula = getinflow − (getsmax − getsmin) gave
        # zero drain rate during uniform-infiltration drought (smax==smin), producing
        # Inf estimate for traps inside a filled parent → no unspill event ever fired.
        @test Set(keys(filltimes_r(seqP))) ⊆ Set(keys(dtP))

        # (2) Topological drain ordering: a child drains no earlier than its parent.
        # Regression: the estimate anchor used cur_amounts[parent].time (drought start)
        # instead of earliest_changetime (when the parent actually fires :unspill).
        # For a trap submerged by a grandparent, these differ; using the stale anchor
        # placed the child's drain estimate before the parent had even started draining.
        nreg = numregions(ts); nt = numtraps(ts)
        for t in (nreg+1):nt
            haskey(dtP, t) || continue
            for child in SWIM.subtrapsof(ts, t)
                haskey(dtP, child) && @test dtP[child] >= dtP[t]
            end
        end

        # (3) Dynamic path drain times match plain path.
        dtM    = draintimes(fill_sequence(ts, weather2; infiltration=infil2,
                                          dyn_traps=[233, 220]))
        common = intersect(keys(dtP), keys(dtM))
        @test !isempty(common)
        @test maximum(abs(dtP[t] - dtM[t]) for t in common) < PARITY_TOL
    end

    @testset "weather-boundary handoff of partially-filled network traps (§10)" begin
        # A weather boundary that changes nothing physically (identical rain rate)
        # must reproduce the single-period result — but only if each network trap's
        # volume is carried across the boundary from its multi-trap ODE state, not a
        # single-trap constant-rate projection (§10).  The projection omits in-network
        # upstream spill (network traps read only terrain inflow, §4), so a partially
        # filled network trap that is fed by an upstream network trap would be handed
        # the wrong boundary volume.  FULL coverage with the boundary placed mid-run
        # exercises this for every partially filled trap at once.
        inc_only(seq) = Iterators.filter(
            e -> e.filled isa Vector{IncrementalUpdate{Bool}}, seq)
        function ft_multi(seq)
            d = Dict{Int,Float64}()
            for e in inc_only(seq), u in e.filled
                u.value && !haskey(d, u.index) && (d[u.index] = e.timestamp)
            end; d
        end
        t_mid = 0.02                     # within (0, 0.0449): the single-period fill span
        w2    = [WeatherEvent(0.0, 1.0), WeatherEvent(t_mid, 1.0)]

        # sanity: the plain path with the (inert) boundary already reproduces ftP
        @test maxΔ(ftP, ft_multi(fill_sequence(ts, w2))) < 1e-9

        # full-coverage dynamic path across the boundary must match to ODE tolerance
        ftF2 = ft_multi(fill_sequence(ts, w2; dyn_traps=collect(1:numtraps(ts))))
        @test Set(keys(ftF2)) == Set(keys(ftP))
        @test maxΔ(ftP, ftF2) < PARITY_TOL
    end
end

@testset "culverts through fill_sequence (end-to-end structural invariants)" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    ts   = spillanalysis(grid, usediags=true)

    # first fill time per trap; robust to weather-boundary full snapshots (Vector{Bool})
    function filltimes(seq)
        d = Dict{Int,Float64}()
        for e in seq[2:end], u in e.filled
            (u isa IncrementalUpdate) && u.value && !haskey(d, u.index) &&
                (d[u.index] = e.timestamp)
        end
        return d
    end
    monotone(seq)     = all(seq[i].timestamp <= seq[i+1].timestamp for i in 1:length(seq)-1)
    total_stored(seq) = sum(fa.amount for fa in amount_at(seq))

    # trap<->trap culvert on the :long chain: inlet (7,119) -> trap 233, outlet
    # (199,4) -> trap 13.  Larger bore (r=1.5) so its effect is well above ODE noise.
    cv = DynCulvert(CartesianIndex(7, 119), CartesianIndex(199, 4), 1.5, 0.6, 0.5, 0.02, 1.7)

    # (1) an empty culvert list is byte-identical to the plain path (culverts inert)
    w1 = [WeatherEvent(0.0, 1.0)]
    @test [e.timestamp for e in fill_sequence(ts, w1; culverts=DynCulvert[])] ==
          [e.timestamp for e in fill_sequence(ts, w1)]

    # (2) rain-then-off; the culvert must complete the full event loop with strictly
    # ordered events (a three-state / spill_path contract violation would throw).
    w  = [WeatherEvent(0.0, 0.05), WeatherEvent(0.3, 0.0)]
    s0 = fill_sequence(ts, w)                    # baseline
    sC = fill_sequence(ts, w; culverts=[cv])     # with culvert
    @test monotone(s0) && monotone(sC)

    # (3) directional behaviour: the network seeds BOTH culvert endpoints, so the inlet trap
    # (233) is an evolving node the culvert draws from — not just filled statically.  With this
    # large bore the culvert bleeds 233 faster than the rain fills it, so 233 never reaches
    # capacity, and it delivers that water to the outlet trap (13), which fills far EARLIER.
    ft0 = filltimes(s0); ftC = filltimes(sC)
    @test haskey(ft0, 233)                     # 233 fills without the culvert
    @test !haskey(ftC, 233)                    # but the culvert keeps it drained below capacity
    @test ftC[13] < ft0[13] - 1e-3             # outlet strongly accelerated

    # (4) mass conservation.  Exact drawn == delivered at the routing layer is covered by the
    # _route_flow tests.  Here the culvert's outlet (trap 13) spills straight out of the domain,
    # so the culvert only moves water OUT faster: the culvert run stores no MORE than the plain
    # run (it never creates water), and the shortfall is bounded by the drained inlet's capacity.
    @test total_stored(sC) <= total_stored(s0) + 1e-6
    @test total_stored(s0) - total_stored(sC) <= SWIM._own_capacity(ts, 233) + 1e-6

    # (5) a different topology (terrain-outlet expansion) also survives the pipeline:
    # inlet (179,37) in a trap, outlet (8,119) on bare terrain traces a fresh chain.
    cv2 = DynCulvert(CartesianIndex(179, 37), CartesianIndex(8, 119), 1.0, 0.6, 0.5, 0.02, 1.7)
    sX  = fill_sequence(ts, w1; culverts=[cv2])
    @test monotone(sX) && length(sX) > 0
end

@testset "solveDynNetwork!: parent at its floor drains/exposes children (:empty)" begin
    grid = loadgrid(joinpath(artifact"swim_testdata", "data", "small", "mini.txt"))
    ts   = spillanalysis(grid, usediags=true)
    CI   = CartesianIndices(size(ts.topography))
    # parent 414 (children 9, 18) as a single subsuming node
    net  = mk_network(ts, [CI[ts.footprints[414][1]]], [9, 18])[1]
    @test [t.trap_ix for t in net.traps] == [414]
    C    = SWIM._build_trap_geometry(ts, net, zeros(size(ts.topography)))[1].capacity
    infil = zeros(size(ts.topography)); infil[ts.footprints[414]] .= 1.0

    # (a) draining from just-below-full down to zero fires :empty (clamped to 0)
    s = [prevfloat(C)]
    ra = solveDynNetwork!(s, ts, net, infil, [0.0])
    @test ra.kind == :empty && ra.trap == 414 && s[1] == 0.0

    # (b) a parent already at its floor (V == 0) with negative net inflow exposes its
    # children immediately: :empty at t = 0, not steady state
    s2 = [0.0]
    rb = solveDynNetwork!(s2, ts, net, infil, [0.0])
    @test rb.kind == :empty && rb.trap == 414 && rb.time == 0.0

    # (c) a parent at its floor with positive inflow fills instead
    s3 = [0.0]
    @test solveDynNetwork!(s3, ts, net, zeros(size(ts.topography)), [1.0]).kind == :fill
end
