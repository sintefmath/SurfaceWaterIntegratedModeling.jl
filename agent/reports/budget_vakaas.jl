# Global water budget on the vakaas webapp case over intervals [a, b] within one weather period:
# S(b) - S(a) = ∫ (rain - outflow - infiltration) dt, with every term from the model's own outputs
# and the rates integrated by composite Simpson's rule, sampled every DT or closer.  The rates
# have kinks (a trap running empty, a wetted cell drying), so a coarse rule over-estimates them:
# three points over 50 min of the dry period gave 8 %, while the volumes match an exact per-trap
# integration to 0.05 %.
#   S        all trap volumes plus NBS layer volumes (all_states_at_timepoints).  Each sampled
#            volume carries the sampler's projection error (an ODE solve per trap, ~0.3 %), which
#            over the terrain's ~18,000 draining traps makes S good to a few m3: hence intervals,
#            not a derivative
#   outflow  the network-aware field (network_watercourses) at the cells where water ends outside
#            any trap: no downstream cell and no trap region, or a sink
#   infil    on a dry cell what the field used (capacity - max(-field, 0)); on a wet cell (submerged
#            by a full trap, or wetted by a non-full trap's water) its whole capacity, i.e. the
#            field's part plus the ponding loss; plus each NBS stack's bottom percolation
# Run from SUrbArea/WebInterface, for K green roofs (0 = base SWIM only):
#   SWIM_INFIL_MMH=20 taskset -c 0,2,4,6,8,10,12,14 julia --project=. <path>/budget_vakaas.jl [K ...]
include("vakaas_webapp_setup.jl")
const G = SWIM.Graphs
# BUDGET_DT and BUDGET_INTERVAL ("a,b") override the sampling and the intervals, to check one closer
const DT = parse(Float64, get(ENV, "BUDGET_DT", "2.0"))  # min, at most between rate samples
intervals = haskey(ENV, "BUDGET_INTERVAL") ? [Tuple(parse.(Float64, split(ENV["BUDGET_INTERVAL"], ",")))] :
            [(0.5, 10.0), (10.0, 20.0), (20.0, 29.5), (30.5, 50.0), (50.0, 70.0), (70.0, 90.0),
             (90.0, 119.5), (120.5, 150.0), (150.0, 200.0), (200.0, 300.0)]

rain_at(t) = w[findlast(e -> e.timestamp <= t, w)].rain_rate
terminal = [c for c in eachindex(ts.topography) if isempty(G.outneighbors(ts.flowgraph, c))]
sinkset  = Set(LinearIndices(size(ts.topography))[s] for s in ts.sinks)
outcells = [c for c in terminal if ts.regions[c] <= 0 || c in sinkset]
@printf("terminal cells %d, of which outflow cells %d (sinks %d)\n",
        length(terminal), length(outcells), count(in(sinkset), outcells))
zvt = SWIM._compute_z_vol_tables(ts)
own(t) = ts.trapvolumes[t] - ts.subvolumes[t]

# the level of trap t's water at own volume V, and its footprint bottoms (a parent's raised to its
# subtraps' spillpoint, as the models do)
function bottoms(t)
    b = Float64.(ts.topography[ts.footprints[t]])
    t > SWIM.numregions(ts) && (b = max.(b, ts.spillpoints[SWIM.subtrapsof(ts, t)[1]].elevation))
    return b
end
level(t, V, b) = SWIM.water_level(SWIM.TrapGeometry(t, own(t), minimum(b), zvt[t][2], zvt[t][1],
                                                    Float64[], Int[], Float64[]), V)

# cells under standing water: the footprints of the topmost full traps, and the wetted part of
# every other trap holding water; cells at or above a trap's spillpoint never pond.  A trap emptied
# by the model can keep a rounding residue (<= 1e-12 m3, ~9,000 traps on vakaas); at any V > 0 its
# level is its bottom, which on flattened terrain wets a whole flat patch, so such traps count as
# empty
const EMPTY = 1e-9
function wet_cells(filled, amounts)
    wet = falses(length(ts.topography))
    for t in eachindex(filled)
        topfull = filled[t] && !any(filled[G.outneighbors(ts.agglomerations, t)])
        (topfull || !filled[t] && amounts[t] > EMPTY) || continue
        sp, b = ts.spillpoints[t].elevation, bottoms(t)
        z = topfull ? sp : level(t, amounts[t], b)
        for (k, c) in enumerate(ts.footprints[t])
            ts.topography[c] < sp && b[k] <= z && (wet[c] = true)
        end
    end
    return wet
end

nbs_percolation(nbs, layers) = sum((last(SWIM._layer_fluxes(nb.system.layers[end], layers[nb.id][end],
                                                            Float64(length(nb.footprint))))
                                    for nb in nbs); init = 0.0)

# the rates at one sampled state: rain, outflow and infiltration (dry cells, wet cells, NBS)
function rates(s, t, nbs)
    r     = rain_at(t)
    field = SWIM.network_watercourses(ts, s.filled, Dict(enumerate(s.amounts)), s.nbs_layers;
                                      nbs = nbs, precipitation = r, infiltration = infil, zvt = zvt)
    wet   = wet_cells(s.filled, s.amounts)
    dry   = sum(isnan(field[c]) || wet[c] ? 0.0 : infil[c] - max(-field[c], 0.0) for c in eachindex(field))
    return (rain = r * length(field), out = sum(max(field[c], 0.0) for c in outcells if !isnan(field[c])),
            dry, wet = sum(infil[c] for c in eachindex(field) if wet[c]), nbs = nbs_percolation(nbs, s.nbs_layers))
end

function budget(K)
    nbs = placements(K)
    seq = SWIM.fill_sequence(ts, w; infiltration = infil, nbs = nbs)
    nbs = SWIM._prepare_nbs!(nbs, ts)                   # ids, as fill_sequence assigned them
    grid(a, b) = range(a, b, length = 2 * ceil(Int, (b - a) / (2DT)) + 1)
    tps = sort(unique(vcat([collect(grid(a, b)) for (a, b) in intervals]...)))
    st  = Dict(zip(tps, SWIM.all_states_at_timepoints(ts, seq, tps; infiltration = infil, nbs = nbs,
                                                      verbose = false)))
    S(t) = sum(st[t].amounts) + sum(sum, values(st[t].nbs_layers); init = 0.0)
    R    = Dict(t => rates(st[t], t, nbs) for t in tps)
    function simpson(f, a, b)                           # composite, on grid(a, b)
        x = grid(a, b)
        return step(x) / 3 * sum((k == 1 || k == length(x) ? 1 : iseven(k) ? 4 : 2) * f(R[x[k]])
                                 for k in eachindex(x))
    end
    println("\n==== K = $K: volumes over each interval, m3")
    println("   interval  |      rain   outflow  infil(dry       wet      NBS) |        dS | residual    rel")
    for (a, b) in intervals
        rain, out = simpson(x -> x.rain, a, b), simpson(x -> x.out, a, b)
        dry, wet, nb = simpson(x -> x.dry, a, b), simpson(x -> x.wet, a, b), simpson(x -> x.nbs, a, b)
        dS  = S(b) - S(a)
        res = rain - out - dry - wet - nb - dS
        @printf("%5.1f-%5.1f | %9.0f %9.0f %9.1f %9.1f %8.2f | %9.1f | %8.1f %7.1e\n", a, b, rain, out,
                dry, wet, nb, dS, res, abs(res) / max(rain, dry + wet + nb))
    end
end

foreach(budget, isempty(ARGS) ? [0, 10] : parse.(Int, ARGS))
