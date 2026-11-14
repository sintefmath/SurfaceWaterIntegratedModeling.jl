# Run fill_sequence on a saved vakaas case (save_vakaas_case.jl) with its K green roofs, in SWIM's
# own environment, and report time, events and allocation.  Prints fill_sequence's progress lines
# (verbose), so watch_progress.py can watch the run.  Run from the repo root, pinned:
#   taskset -c 0,2,4,6,8,10,12,14 julia --project=. agent/reports/run_vakaas_case.jl [CASE]
using SurfaceWaterIntegratedModeling, JLD2
const SWIM = SurfaceWaterIntegratedModeling

d = JLD2.load(isempty(ARGS) ? "/tmp/vakaas_case.jld2" : ARGS[1])
ts, infil = d["ts"], d["infil"]
roof() = SWIM.elhadiGreenRoof(10.0, 5.0, 0.5, 0.5, 1.0, 1.0, 0.0, 1.0)
nbs = [SWIM.DynNBSPlacement(roof(), fp, 2, CartesianIndex{2}[]) for fp in d["footprints"]]
rain = 120.0 / (1000 * 60)
w = [WeatherEvent(0.0, rain), WeatherEvent(30.0, 0.0), WeatherEvent(120.0, 0.0)]
fill_sequence(ts, w; infiltration = infil, nbs = nbs[1:1])          # compile first
stats = @timed fill_sequence(ts, w; infiltration = infil, nbs = nbs, verbose = true)
println("RESULT $(length(nbs)) roofs: ", round(stats.time, digits = 1), " s, events ",
        length(stats.value), ", alloc ", round(stats.bytes / 2^30, digits = 2), " GB")
