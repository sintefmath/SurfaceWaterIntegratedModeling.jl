# Where does fill_sequence spend its time on vakaas with K green roofs?  Builds the webapp case as
# perf_vakaas_webapp.jl does, then samples fill_sequence with the Profile stdlib for a fixed wall
# time and writes a flat profile, sorted by sample count, to OUT.
#
# Needs a second thread (`--threads=1,1`) for the task that stops the sampling and writes the
# result.  Run from SUrbArea/WebInterface:
#   taskset -c 0,2,4,6,8,10,12,14 julia --threads=1,1 --project=. <path>/perf_vakaas_profile.jl [K] [seconds] [OUT]
ENV["GENIE_ENV"] = "test"
push!(LOAD_PATH, abspath("src"))
using Genie
Genie.loadapp()
using Profile
const EC   = Main.WebInterface.ExplorersController
const SWIM = EC.SWIM
const JLD2 = EC.JLD2

K    = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
SECS = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 900.0
OUT  = length(ARGS) >= 3 ? ARGS[3] : "/tmp/perf_vakaas_profile.txt"

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

rain = 120.0 / (1000 * 60)
w = [SWIM.WeatherEvent(0.0, rain), SWIM.WeatherEvent(30.0, 0.0), SWIM.WeatherEvent(120.0, 0.0)]

# Stop after SECS of sampling: write the flat profile and end the process.  The main thread is
# busy in fill_sequence and never yields, so this runs on the default-pool thread (in Julia 1.12
# the main thread is the interactive one) and waits with an OS sleep: Julia's `sleep` needs the
# event loop, which is not serviced while the main thread never yields.
Threads.nthreads(:default) >= 1 || error("start julia with --threads=1,1")
Threads.@spawn :default begin
    t_end = time() + SECS
    while time() < t_end
        Libc.systemsleep(2.0)
    end
    Profile.stop_timer()
    open(OUT, "w") do io
        println(io, "K = $K green roofs, $(SECS) s of fill_sequence sampled")
        Profile.print(IOContext(io, :displaysize => (10000, 220)); format = :flat,
                      sortedby = :count, mincount = 200, noisefloor = 0)
    end
    ccall(:_exit, Cvoid, (Cint,), 0)
end

Profile.init(n = 10^8, delay = 0.01)
println("sampling fill_sequence with $K roofs for $SECS s"); flush(stdout)
@profile SWIM.fill_sequence(ts, w; infiltration = infil, nbs = nbs)
Profile.stop_timer()
println("fill_sequence finished before the time limit")
open(OUT, "w") do io
    Profile.print(IOContext(io, :displaysize => (10000, 220)); format = :flat,
                  sortedby = :count, mincount = 200, noisefloor = 0)
end
