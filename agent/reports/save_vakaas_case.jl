# Build the vakaas webapp case (as perf_vakaas_webapp.jl) with K green roofs and save it, so
# debugging runs can load it in SWIM's own environment without starting the webapp.
# Run from SUrbArea/WebInterface, pinned:
#   taskset -c 0,2,4,6,8,10,12,14 julia --project=. <path>/save_vakaas_case.jl [K] [OUT]
ENV["GENIE_ENV"] = "test"
push!(LOAD_PATH, abspath("src"))
using Genie
Genie.loadapp()
const EC   = Main.WebInterface.ExplorersController
const SWIM = EC.SWIM
const JLD2 = EC.JLD2

K   = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
OUT = length(ARGS) >= 2 ? ARGS[2] : "/tmp/vakaas_case.jld2"

ws = Dict{String,Any}("model_name" => "vakas")
ws["bbox"], = EC.boundingbox_and_downsampling_suggestion("vakas")
ws["resolution"] = size(JLD2.load(joinpath("model_data", "geometry", "vakaas_dtm_flattened.jld2"), "dtm"))
ws["jsondata"] = EC.load_model_jsondata("vakas")
selection = Dict{String,Any}(
    "buildings" => EC.BUILDING_LAYER_NAMES, "landcover" => EC.WATER_BODY_NAMES,
    "drainage" => EC.DRAINAGE_NAMES, "barriers" => String[], "permeable" => EC.PERMEABLE_NAMES,
    "permeable_rates" => Dict(n => 20.0 for n in EC.PERMEABLE_NAMES),
    "infiltration_groups" => String[], "effective_culverts" => EC.EFFECTIVE_CULVERT_NAMES)
cache = EC._compute_layer_cache(ws, selection)
ts    = EC._compute_tstruct(ws, selection, cache)
infil = get(cache, "infiltration_rate_grid", zeros(ws["resolution"]))

ids = cache["building_ids"]
LI  = LinearIndices(size(ids))
cells_of = Dict{Int,Vector{Int}}()
for c in eachindex(ids)
    ids[c] > 0 && push!(get!(cells_of, ids[c], Int[]), LI[c])
end
largest = sort(collect(keys(cells_of)); by = b -> -length(cells_of[b]))
footprints = [cells_of[b] for b in largest[1:K]]

JLD2.jldsave(OUT; ts, infil, footprints)
println("saved $OUT")
