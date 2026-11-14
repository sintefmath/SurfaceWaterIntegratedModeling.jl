# Profile SWIM on the SUrbArea "vakaas" terrain (3162 x 3689 = 11.7M cells), stage by stage:
# time, allocations and peak memory (RSS) per phase, plus the resident size of the result.
#
# Stage 1 (this script): bare flattened DTM, as cached by the webapp - spillanalysis, then
# fill_sequence for a rain event, then a few timepoint queries.
#
# Run from the repo root (pinned to the P-cores):
#   taskset -c 0,2,4,6,8,10,12,14 julia --project=. agent/reports/perf_vakaas.jl
# The SUrbArea checkout is assumed next to this repository (override with SURBAREA_DIR).
using SurfaceWaterIntegratedModeling, JLD2, Printf
const SWIM = SurfaceWaterIntegratedModeling

const SURB = get(ENV, "SURBAREA_DIR", normpath(joinpath(@__DIR__, "..", "..", "..", "SUrbArea")))
const DTM  = joinpath(SURB, "WebInterface", "model_data", "geometry", "vakaas_dtm_flattened.jld2")

gb(x) = x / 2^30

# resident bytes, counted directly (Base.summarysize tracks every object it visits in an IdDict,
# which takes very long for millions of small vectors)
_bytes(x::AbstractArray{<:Number}) = sizeof(x)
_bytes(x::AbstractArray) = sizeof(x) + sum(_bytes, x; init = 0)
_bytes(x::SWIM.Graphs.SimpleDiGraph) = _bytes(x.fadjlist) + _bytes(x.badjlist)
_bytes(x) = Base.summarysize(x)
function phase(name, f)
    GC.gc()
    rss0 = Sys.maxrss()
    stats = @timed f()
    @printf("%-28s %8.1f s  alloc %7.2f GB  gc %5.1f%%  peak RSS %6.2f GB (+%.2f)\n",
            name, stats.time, gb(stats.bytes), 100 * stats.gctime / max(stats.time, eps()),
            gb(Sys.maxrss()), gb(Sys.maxrss() - rss0))
    flush(stdout)
    return stats.value
end

dtm = phase("load DTM", () -> JLD2.load(DTM, "dtm"))
println("grid ", size(dtm), " = ", round(length(dtm) / 1e6, digits = 2), "M cells")

ts = phase("spillanalysis", () -> spillanalysis(dtm))
@printf("traps %d, regions %d, hierarchy depth %d\n", numtraps(ts), numregions(ts),
        maximum(length, ts.supertraps_of))
flush(stdout)
println("resident size of TrapStructure fields (GB):")
for f in fieldnames(typeof(ts))
    @printf("  %-22s %7.3f\n", f, gb(_bytes(getfield(ts, f))))
end
flush(stdout)

# rain: 120 mm/h for 30 min, then dry for 90 min (rates in m/min, as the webapp uses)
rain = 120.0 / (1000 * 60)
w = [WeatherEvent(0.0, rain), WeatherEvent(30.0, 0.0), WeatherEvent(120.0, 0.0)]
seq = phase("fill_sequence", () -> fill_sequence(ts, w; verbose = false))
@printf("events %d, filled at end %d\n", length(seq), count(filled_at(seq, length(seq))))
nfull = count(ev -> !(eltype(ev.runoff) <: SWIM.IncrementalUpdate), seq)
@printf("full snapshots in the sequence: %d of %d (each holds an %.2f GB runoff grid)\n",
        nfull, length(seq), gb(sizeof(dtm)))
flush(stdout)

tps = collect(10.0:10.0:120.0)
phase("trap_states (12 timepoints)", () -> trap_states_at_timepoints(ts, seq, tps; verbose = false))
phase("watercourses (one call)", () -> watercourses(ts, filled_at(seq, length(seq)); precipitation = rain))
