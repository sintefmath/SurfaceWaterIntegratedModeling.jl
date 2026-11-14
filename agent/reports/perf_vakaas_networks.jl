# How many dynamic networks do K green roofs on vakaas give, and how big are they?  Builds the
# webapp case (as perf_vakaas_webapp.jl), runs fill_sequence without NBS for the filled states,
# and reports the networks setup_network builds for the roofs at the start, at the end of the
# rain and at the end.  Run from SUrbArea/WebInterface, pinned to the P-cores:
#   taskset -c 0,2,4,6,8,10,12,14 julia --project=. <path>/perf_vakaas_networks.jl [K]
ENV["GENIE_ENV"] = "test"
push!(LOAD_PATH, abspath("src"))
using Genie
Genie.loadapp()
using Printf
const EC   = Main.WebInterface.ExplorersController
const SWIM = EC.SWIM
const JLD2 = EC.JLD2

K = isempty(ARGS) ? 10 : parse(Int, ARGS[1])

const MODEL = "vakas"
ws = Dict{String,Any}("model_name" => MODEL)
ws["bbox"], = EC.boundingbox_and_downsampling_suggestion(MODEL)
ws["resolution"] = size(JLD2.load(joinpath("model_data", "geometry", "vakaas_dtm_flattened.jld2"), "dtm"))
ws["jsondata"] = EC.load_model_jsondata(MODEL)
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
roof() = SWIM.elhadiGreenRoof(10.0, 5.0, 0.5, 0.5, 1.0, 1.0, 0.0, 1.0)
nbs = [SWIM.DynNBSPlacement(roof(), cells_of[b], 2, CartesianIndex{2}[]) for b in largest[1:K]]
SWIM._prepare_nbs!(nbs, ts)

rain = 120.0 / (1000 * 60)
w = [SWIM.WeatherEvent(0.0, rain), SWIM.WeatherEvent(30.0, 0.0), SWIM.WeatherEvent(120.0, 0.0)]
seq = SWIM.fill_sequence(ts, w; infiltration = infil)
at(t) = SWIM.filled_at(seq, something(findlast(ev -> ev.timestamp <= t, seq), 1))

println("\nK = $K green roofs; roof sizes (cells): ", [length(n.footprint) for n in nbs])
for (label, t) in (("start", 0.0), ("end of rain (t=30)", 30.0), ("end (t=120)", 120.0))
    filled = at(t)
    comps  = SWIM.setup_network(ts, findall(filled); nbs = nbs)
    @printf("\n%s: %d full traps in the terrain, %d network(s)\n", label, count(filled), length(comps))
    println("  net | NBS | traps | covered | paths | path cells | footprint cells of its traps")
    for (i, net) in enumerate(comps)
        covered = sum(1 + length(SWIM._descendants(ts, t.trap_ix)) for t in net.traps; init = 0)
        pcells  = sum(length(p.cells) for p in net.flow_paths; init = 0)
        fcells  = sum(length(ts.footprints[t.trap_ix]) for t in net.traps; init = 0)
        @printf("  %3d | %3d | %5d | %7d | %5d | %10d | %d\n", i, length(net.nbs), length(net.traps),
                covered, length(net.flow_paths), pcells, fcells)
    end
end
