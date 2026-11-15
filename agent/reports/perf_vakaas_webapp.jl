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
include("vakaas_webapp_setup.jl")

seq = phase("fill_sequence (no NBS)", () -> SWIM.fill_sequence(ts, w; infiltration = infil))
@printf("events %d\n", length(seq))
flush(stdout)

for K in (isempty(ARGS) ? [1, 10] : parse.(Int, ARGS))
    nbs = placements(K)
    @printf("K = %d green roofs, %d footprint cells\n", K, sum(length(p.footprint) for p in nbs))
    # SWIM_VERBOSE=1 prints fill_sequence's progress lines, for the progress monitor
    s = phase("fill_sequence (K=$K NBS)", () -> SWIM.fill_sequence(ts, w; infiltration = infil, nbs = nbs,
                                                                  verbose = get(ENV, "SWIM_VERBOSE", "0") == "1"))
    @printf("events %d\n", length(s))
    flush(stdout)
end
