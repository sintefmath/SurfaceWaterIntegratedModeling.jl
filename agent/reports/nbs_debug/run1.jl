# First pass: build both synthetic setups, run fill_sequence, sample gauges and the
# flow field over time, and flag suspicious behaviour.
include("setup.jl")
using Printf

dtm = load_dtm_flat()
gs = size(dtm)
println("grid ", gs, " z range ", extrema(dtm))
t0 = time()
tstruct = spillanalysis(dtm; sinks = fill(false, gs...))
println("spillanalysis: ", round(time() - t0, digits=2), " s; traps=", numtraps(tstruct),
        " regions=", numregions(tstruct))

const RAIN = [(120.0, 1000.0)]   # mm/h, minutes (matches gauge CSV levels: 0.0333 l/s per cell)
const TMAX = 1000.0

function analyse(name, file)
    println("\n================ ", name)
    feats = load_features(file)
    nbs = build_placements(feats, gs)
    gp  = gauges(feats, gs)
    for n in nbs
        println("NBS id=", n.id, " cells=", length(n.footprint), " outlets=", n.outlets,
                " bbox rows ", extrema(CartesianIndices(gs)[n.footprint] .|> x->x[1]),
                " cols ", extrema(CartesianIndices(gs)[n.footprint] .|> x->x[2]))
    end
    println("gauges: ", gp)
    we = weather(RAIN)
    t0 = time()
    seq = fill_sequence(tstruct, we; nbs = nbs)
    println("fill_sequence: ", round(time() - t0, digits=2), " s, events=", length(seq))
    seq0 = fill_sequence(tstruct, we)
    println("no-NBS events=", length(seq0))
    return (; nbs, gp, seq, seq0, we)
end

r1 = analyse("setup 1 (one roof)",  joinpath(SURB, "debug/synthetic_userdefined.geojson"))
r2 = analyse("setup 2 (three roofs)", joinpath(SURB, "debug/synthetic_userdefined_2.geojson"))

JLD2.@save "run1.jld2" r1 r2 tstruct dtm
