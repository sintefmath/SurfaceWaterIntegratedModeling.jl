using Test, SurfaceWaterIntegratedModeling, LazyArtifacts

const SWIM = SurfaceWaterIntegratedModeling

# ---------------------------------------------------------------------------
# network_watercourses: flow-intensity field at one timepoint, accounting for the
# dynamic network.  NBS scenario mirrors the E2 driver test — a west-flowing plane
# with a pit near the west edge and a slow `puddle` NBS on the flow path into it.
# ---------------------------------------------------------------------------
@testset "network_watercourses: NBS retention field" begin
    N    = 30
    grid = Float64[j for i in 1:N, j in 1:N]
    for i in 12:18, j in 3:6; grid[i, j] = 1.0; end
    ts   = spillanalysis(grid, usediags = true)
    LI   = LinearIndices(size(grid))
    fp   = [LI[CartesianIndex(i, j)] for i in 12:18 for j in 15:20]
    pl   = SWIM.DynNBSPlacement(SWIM.puddle(50.0; kOUT = 0.01), fp, CartesianIndex{2}[])
    w    = [SWIM.WeatherEvent(0.0, 1.0e-3)]

    seq  = fill_sequence(ts, w; nbs = [pl])
    tp   = seq[end].timestamp * 0.5 + 1.0            # mid-run: the NBS is holding water

    allst = all_states_at_timepoints(ts, seq, [tp]; nbs = [pl], verbose = false)[1]
    netst = network_states_at_timepoints(ts, seq, [tp]; nbs = [pl], verbose = false)[1]
    @test netst.nbs_layers[pl.id][1] > 1e-3          # storage is nonzero at this instant

    prec = 1.0e-3
    ro, = SWIM.watercourses(ts, allst.filled; precipitation = prec, infiltration = 0.0)
    rn  = network_watercourses(ts, allst.filled, netst.trap_volumes, netst.nbs_layers;
                               nbs = [pl], precipitation = prec, infiltration = 0.0)

    @test size(rn) == size(grid)
    @test all(isnan, rn[fp])                          # footprint: undefined surface flow
    offfoot = setdiff(1:length(rn), fp)
    @test all(isfinite, rn[offfoot])                  # everywhere else is defined
    @test all(!isnan, ro[fp])                         # the oblivious field does route the footprint

    # the NBS holds runoff back, so downstream of it the network field carries LESS flow than the
    # oblivious one, and never more
    probe = LI[CartesianIndex(15, 10)]                # on the flow path, west (downstream) of the NBS
    @test rn[probe] < ro[probe]
    @test sum(x -> max(x, 0.0), rn[offfoot]) < sum(x -> max(x, 0.0), ro[offfoot])
end

# ---------------------------------------------------------------------------
# Culvert: diverts flow from its inlet to its outlet.  A west-flowing plane; a culvert
# carries flow from a cell on the flow path to another cell downstream.  The field must
# conserve mass — the total leaving the domain is unchanged by an internal diversion.
# ---------------------------------------------------------------------------
@testset "network_watercourses: culvert diverts, conserves mass" begin
    N    = 20
    grid = Float64[j for i in 1:N, j in 1:N]          # rising east → flow west, off the j=1 edge
    ts   = spillanalysis(grid, usediags = true)
    @test numtraps(ts) == 0                            # pure slope, no traps
    prec = 1.0e-3

    inlet  = CartesianIndex(5, 15)
    outlet = CartesianIndex(15, 8)                     # a lower (further west) cell
    cv     = DynCulvert(ts, inlet, outlet; r = 0.5)

    full_traps  = Vector{Bool}(ts.trapvolumes .== 0.0)
    trap_volumes = Dict{Int,Float64}()
    nbs_layers   = Dict{Int,Vector{Float64}}()

    ro, _, off_ro, _ = SWIM.watercourses(ts, full_traps; precipitation = prec, infiltration = 0.0)
    rn = network_watercourses(ts, full_traps, trap_volumes, nbs_layers;
                              culverts = [cv], precipitation = prec, infiltration = 0.0)

    @test all(isfinite, rn)
    li = LinearIndices(size(grid))
    # flow is drawn just downstream of the inlet and re-appears at the outlet
    @test rn[li[inlet]] < ro[li[inlet]]                # inlet cell loses the drawn flow
    # mass conservation: total off-domain outflow is unchanged (diversion is internal).  Recompute
    # the off-domain total for the network field by summing what leaves the west edge.
    edge = [li[CartesianIndex(i, 1)] for i in 1:N]
    # every cell drains west to j=1; the domain-exit total is the sum over the west edge column's
    # accumulated flow.  It must match the oblivious total to within rounding.
    @test isapprox(sum(max.(rn[edge], 0.0)), sum(max.(ro[edge], 0.0)); rtol = 1e-9)
end
