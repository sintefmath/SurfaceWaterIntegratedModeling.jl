# Save the spillanalysis inputs of a webapp model as the webapp assembles them (_compute_tstruct:
# buildings raised, waterbodies, courtyard and manhole sinks, culverts), so the trap structure
# can be built from them in SWIM's own environment, e.g. by two code versions.
# Run from SUrbArea/WebInterface, pinned:
#   taskset -c 0,2,4,6,8,10,12,14 julia --project=. <path>/save_vakaas_inputs.jl [MODEL] [OUT]
# MODEL is the webapp's model name (default "vakas"); layers a model lacks are skipped.
ENV["GENIE_ENV"] = "test"
push!(LOAD_PATH, abspath("src"))
using Genie
Genie.loadapp()
const EC   = Main.WebInterface.ExplorersController
const SWIM = EC.SWIM
const JLD2 = EC.JLD2
MODEL = length(ARGS) >= 1 ? ARGS[1] : "vakas"
OUT   = length(ARGS) >= 2 ? ARGS[2] : "/tmp/$(MODEL)_inputs.jld2"

ws = Dict{String,Any}("model_name" => MODEL)
ws["bbox"], = EC.boundingbox_and_downsampling_suggestion(MODEL)
dtm = JLD2.load(joinpath("model_data", "geometry", "$(EC.model_file_prefix(MODEL))_dtm_flattened.jld2"), "dtm")
ws["resolution"] = size(dtm)
ws["jsondata"] = EC.load_model_jsondata(MODEL)
selection = Dict{String,Any}(
    "buildings" => EC.BUILDING_LAYER_NAMES, "landcover" => EC.WATER_BODY_NAMES,
    "drainage" => EC.DRAINAGE_NAMES, "barriers" => String[], "permeable" => EC.PERMEABLE_NAMES,
    "permeable_rates" => Dict(n => 20.0 for n in EC.PERMEABLE_NAMES),
    "infiltration_groups" => String[], "effective_culverts" => EC.EFFECTIVE_CULVERT_NAMES)
cache = EC._compute_layer_cache(ws, selection)

# as _compute_tstruct
grid = copy(dtm)
get(cache, "building_grid", nothing) === nothing ||
    SWIM.raise_buildings!(grid, cache["building_grid"], 10.0)
waterbodies = get(cache, "water_grid", nothing)
sinks = fill(false, size(dtm))
for cells in get(cache, "courtyards", Vector{CartesianIndex{2}}[])
    sinks[cells[argmax(grid[cells])]] = true
end
haskey(cache, "drains_grid") && (sinks .|= (cache["drains_grid"] .== EC.MANHOLE_DRAIN_CODE))
culverts = EC._unique_culverts(vcat(get(cache, "culverts", Tuple{CartesianIndex{2},CartesianIndex{2}}[]),
                                    get(cache, "effective_culverts", Tuple{CartesianIndex{2},CartesianIndex{2}}[])))
JLD2.jldsave(OUT; grid, waterbodies, sinks, culverts)
println("saved $OUT: $(count(sinks)) sinks, $(length(culverts)) culverts")
