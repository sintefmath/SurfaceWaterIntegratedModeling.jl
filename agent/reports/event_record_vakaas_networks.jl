# Record fill_sequence's event history on the vakaas webapp case with networks: the 10-roof
# networks' traps as dyn_traps (no NBS), and the 10 roofs; the network counterpart of
# event_record.jl, read by event_record_compare.jl.  Run from SUrbArea/WebInterface:
#   SWIM_INFIL_MMH=20 taskset -c 0,2,4,6,8,10,12,14 julia --project=. <path>/event_record_vakaas_networks.jl OUT
include("vakaas_webapp_setup.jl")
record(seq) = [(ev.timestamp, [(u.index, u.value) for u in filled_updates_of(ev)]) for ev in seq]
filled_updates_of(ev) = eltype(ev.filled) <: SWIM.IncrementalUpdate ? ev.filled :
                        [SWIM.IncrementalUpdate(i, v) for (i, v) in enumerate(ev.filled)]

nbs = placements(10)
dyn = sort!(unique([t.trap_ix for net in SWIM.setup_network(ts, Int[]; nbs) for t in net.traps]))
results = Dict{String,Any}()
for (name, kw) in [("vakaas dyn_traps", (; dyn_traps = dyn)), ("vakaas 10 roofs", (; nbs))]
    stats = @timed SWIM.fill_sequence(ts, w; infiltration = infil, kw...)
    @printf("%-28s %7.2f s  %6d events\n", name, stats.time, length(stats.value))
    results[name] = (stats.time, record(stats.value))
end
JLD2.save(ARGS[1], "results", results)
