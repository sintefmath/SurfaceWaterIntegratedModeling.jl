# Profile SWIM on the SUrbArea "vakaas" terrain as the webapp sets it up (stages 2 and 3 of the
# large-case profile; stage 1, the bare terrain, is perf_vakaas.jl).
#
# Stage 2: the webapp's own layer cache and trap structure (buildings raised, waterbodies,
#          manhole and courtyard sinks, culverts, permeable infiltration), then fill_sequence.
# Stage 3: the same with green roofs on the K largest buildings, for several K.
#
# Loads the webapp headless (as its tests do) and calls its _compute_layer_cache /
# _compute_tstruct, so the case is built exactly as in the app.  Run from SUrbArea/WebInterface:
#   taskset -c 0,2,4,6,8,10,12,14 julia --project=. <path>/perf_vakaas_webapp.jl [K1 K2 ...]
ENV["GENIE_ENV"] = "test"
push!(LOAD_PATH, abspath("src"))
using Genie
Genie.loadapp()
using Printf
const EC   = Main.WebInterface.ExplorersController
const SWIM = EC.SWIM
const JLD2 = EC.JLD2

gb(x) = x / 2^30
function phase(name, f)
    GC.gc()
    rss0 = Sys.maxrss()
    stats = @timed f()
    @printf("%-34s %8.1f s  alloc %7.2f GB  gc %5.1f%%  peak RSS %6.2f GB (+%.2f)\n",
            name, stats.time, gb(stats.bytes), 100 * stats.gctime / max(stats.time, eps()),
            gb(Sys.maxrss()), gb(Sys.maxrss() - rss0))
    flush(stdout)
    return stats.value
end

# ---- the workspace the app would have after loading the model --------------------------------
const MODEL = "vakas"                                  # URL name; files are prefixed "vakaas"
ws = Dict{String,Any}("model_name" => MODEL)
ws["bbox"], = EC.boundingbox_and_downsampling_suggestion(MODEL)
ws["resolution"] = size(JLD2.load(joinpath("model_data", "geometry", "vakaas_dtm_flattened.jld2"), "dtm"))
ws["jsondata"] = phase("load FKB json", () -> EC.load_model_jsondata(MODEL))
for (k, v) in ws["jsondata"]; println("  ", k, ": ", length(v)); end

selection = Dict{String,Any}(
    "buildings"           => EC.BUILDING_LAYER_NAMES,
    "landcover"           => EC.WATER_BODY_NAMES,
    "drainage"            => EC.DRAINAGE_NAMES,           # culverts and manholes
    "barriers"            => String[],
    "permeable"           => EC.PERMEABLE_NAMES,
    "permeable_rates"     => Dict(n => 20.0 for n in EC.PERMEABLE_NAMES),   # mm/h
    "infiltration_groups" => String[],
    "effective_culverts"  => EC.EFFECTIVE_CULVERT_NAMES)

cache = phase("layer cache (rasterize FKB)", () -> EC._compute_layer_cache(ws, selection))
println("  layers: ", join(sort(collect(string.(keys(cache)))), ", "))
ts = phase("spillanalysis (webapp setup)", () -> EC._compute_tstruct(ws, selection, cache))
@printf("traps %d, regions %d, hierarchy depth %d, sinks %d\n", SWIM.numtraps(ts),
        SWIM.numregions(ts), maximum(length, ts.supertraps_of), length(ts.sinks))
flush(stdout)

infil = get(cache, "infiltration_rate_grid", zeros(ws["resolution"]))
@printf("infiltrating cells: %d\n", count(>(0), infil))
rain = 120.0 / (1000 * 60)                             # 120 mm/h, in m/min as the app uses
w = [SWIM.WeatherEvent(0.0, rain), SWIM.WeatherEvent(30.0, 0.0), SWIM.WeatherEvent(120.0, 0.0)]
seq = phase("fill_sequence (no NBS)", () -> SWIM.fill_sequence(ts, w; infiltration = infil))
@printf("events %d\n", length(seq))
flush(stdout)

# ---- stage 3: green roofs on the K largest buildings -----------------------------------------
ids = cache["building_ids"]
LI  = LinearIndices(size(ids))
cells_of = Dict{Int,Vector{Int}}()
for c in eachindex(ids)
    ids[c] > 0 && push!(get!(cells_of, ids[c], Int[]), LI[c])
end
largest = sort(collect(keys(cells_of)); by = b -> -length(cells_of[b]))
roof() = SWIM.elhadiGreenRoof(10.0, 5.0, 0.5, 0.5, 1.0, 1.0, 0.0, 1.0)
# both layers overflow onto the terrain (n_terrain = 2), so no drawn pipe outlet is needed
placements(K) = [SWIM.DynNBSPlacement(roof(), cells_of[b], 2, CartesianIndex{2}[]) for b in largest[1:K]]

for K in (isempty(ARGS) ? [1, 10] : parse.(Int, ARGS))
    nbs = placements(K)
    @printf("K = %d green roofs, %d footprint cells\n", K, sum(length(p.footprint) for p in nbs))
    # SWIM_VERBOSE=1 prints fill_sequence's progress lines, for the progress monitor
    s = phase("fill_sequence (K=$K NBS)", () -> SWIM.fill_sequence(ts, w; infiltration = infil, nbs = nbs,
                                                                  verbose = get(ENV, "SWIM_VERBOSE", "0") == "1"))
    @printf("events %d\n", length(s))
    flush(stdout)
end
