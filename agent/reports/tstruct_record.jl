# Build the TrapStructure for a set of cases and record a fingerprint of every field (graphs as
# adjacency lists), with spillanalysis timings, to check that a change leaves the static analysis
# bit-identical.
#   taskset -c 0,2,4,6,8,10,12,14 julia --project=. agent/reports/tstruct_record.jl OUT [FILES...]
# FILES: .jld2 files with either a bare DTM (key `dtm`, e.g. the webapp's *_dtm_flattened.jld2) or
# a webapp model's spillanalysis inputs (save_vakaas_inputs.jl); each becomes a case.
using SurfaceWaterIntegratedModeling, LazyArtifacts, JLD2, Printf
const SWIM = SurfaceWaterIntegratedModeling

# one fingerprint per field (graphs by their adjacency lists): hashes, not the values, so a
# record of large terrains stays small enough to load two of them for comparing
function fields(ts)
    # structs without their own hash (Spillpoint) hash by objectid, which differs between runs
    content(v) = v isa SWIM.Graphs.SimpleDiGraph ? (v.fadjlist, v.badjlist) :
                 v isa AbstractVector{<:SWIM.Spillpoint} ? [ntuple(i -> getfield(s, i), 4) for s in v] : v
    fingerprint(v) = (hash(content(v)), applicable(length, v) ? length(v) : 0)
    return Dict(string(f) => fingerprint(v)
                for f in fieldnames(typeof(ts)) for v in (getfield(ts, f),))
end

function case(name, args...; kw...)
    spillanalysis(args...; kw...)                                   # compile / warm up
    t = minimum(@elapsed(spillanalysis(args...; kw...)) for _ in 1:2)
    ts = spillanalysis(args...; kw...)
    @printf("%-22s %7.2f s  %6d traps\n", name, t, numtraps(ts)); flush(stdout)
    return name => (t, fields(ts))
end

data = joinpath(artifact"swim_testdata", "data")
mini, bay = SWIM.loadgrid(joinpath(data, "small/mini.txt")), SWIM.loadgrid(joinpath(data, "small/bay.txt"))
mask = SWIM.loadgrid(joinpath(data, "small/baybuildings.txt")) .|> Bool
syn  = SWIM.loadgrid(joinpath(data, "synthetic/synsurf.txt"))
results = Dict{String,Any}(case("mini", mini), case("mini no diags", mini; usediags = false),
                           case("bay", bay), case("bay mask", bay; clip_mask = mask),
                           case("synsurf", syn))
for f in ARGS[2:end]
    d, name = JLD2.load(f), replace(basename(f), ".jld2" => "")
    push!(results, haskey(d, "dtm") ? case(name, d["dtm"]) :
          case(name, d["grid"]; waterbodies = d["waterbodies"], sinks = d["sinks"],
               culverts = d["culverts"]))
end
JLD2.jldsave(ARGS[1]; results)
println("written $(ARGS[1])")
