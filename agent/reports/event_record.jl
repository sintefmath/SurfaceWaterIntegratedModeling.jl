# Record fill_sequence's complete event history (time, and which traps changed and how) for a set
# of cases, with timings, to check that a change leaves the results bit-identical.
#   taskset -c 0,2,4,6,8,10,12,14 julia --project=. agent/reports/event_record.jl OUT [CASE]
# OUT: a .jld2 file; CASE: a saved vakaas case (save_vakaas_case.jl), included when given.
using SurfaceWaterIntegratedModeling, LazyArtifacts, JLD2, Printf
const SWIM = SurfaceWaterIntegratedModeling

record(seq) = [(ev.timestamp, [(u.index, u.value) for u in filled_updates_of(ev)]) for ev in seq]
filled_updates_of(ev) = eltype(ev.filled) <: SWIM.IncrementalUpdate ? ev.filled :
                        [SWIM.IncrementalUpdate(i, v) for (i, v) in enumerate(ev.filled)]

function run(name, ts, w; kw...)
    fill_sequence(ts, w; kw...)                                       # compile / warm up
    stats = @timed fill_sequence(ts, w; kw...)
    @printf("%-28s %7.2f s  %6d events\n", name, stats.time, length(stats.value))
    flush(stdout)
    return name => (stats.time, record(stats.value))
end

out  = ARGS[1]
data = joinpath(artifact"swim_testdata", "data")
grids = Dict("mini" => "small/mini.txt", "bay" => "small/bay.txt", "synsurf" => "synthetic/synsurf.txt")
results = Dict{String,Any}()
for (name, f) in sort(collect(grids))
    ts = spillanalysis(SWIM.loadgrid(joinpath(data, f)))
    infil = 0.5 * ones(size(ts.topography))
    push!(results, run("$name", ts, [WeatherEvent(0.0, 1.0)]))
    push!(results, run("$name infiltration", ts, [WeatherEvent(0.0, 1.0), WeatherEvent(1.0, 0.0)];
                       infiltration = infil))
end
if length(ARGS) >= 2
    d = JLD2.load(ARGS[2])
    ts, infil = d["ts"], d["infil"]
    roof() = SWIM.elhadiGreenRoof(10.0, 5.0, 0.5, 0.5, 1.0, 1.0, 0.0, 1.0)
    nbs = [SWIM.DynNBSPlacement(roof(), fp, 2, CartesianIndex{2}[]) for fp in d["footprints"]]
    rain = 120.0 / (1000 * 60)
    w = [WeatherEvent(0.0, rain), WeatherEvent(30.0, 0.0), WeatherEvent(120.0, 0.0)]
    push!(results, run("vakaas no NBS", ts, w; infiltration = infil))
    push!(results, run("vakaas 10 roofs", ts, w; infiltration = infil, nbs = nbs))
end
JLD2.jldsave(out; results)
println("written $out")
