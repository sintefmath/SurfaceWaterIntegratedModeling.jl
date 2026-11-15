# Base SWIM against the dynamic networks on the same physics, on vakaas: fill_sequence without NBS,
# once on the constant-rate path and once with the traps of the 10-roof networks solved as networks
# (`dyn_traps`).  Both solve the same equations, so the same traps should fill, at the same times,
# and the trap volumes at sampled times should agree to solver tolerance.  The large-terrain
# version of the test "fill_sequence dyn_traps parity".  Run from SUrbArea/WebInterface:
#   SWIM_INFIL_MMH=20 taskset -c 0,2,4,6,8,10,12,14 julia --project=. <path>/parity_vakaas.jl
include("vakaas_webapp_setup.jl")

# the traps of the networks the 10 roofs start with (nothing full yet)
dyn = sort!(unique([t.trap_ix for net in SWIM.setup_network(ts, Int[]; nbs = placements(10))
                              for t in net.traps]))
@printf("dyn_traps: %d traps\n", length(dyn))

base = phase("fill_sequence, base path", () -> SWIM.fill_sequence(ts, w; infiltration = infil))
nets = phase("fill_sequence, networks", () -> SWIM.fill_sequence(ts, w; infiltration = infil,
                                                                 dyn_traps = dyn))
@printf("events: base %d, networks %d\n", length(base), length(nets))

# trap -> first time it is filled, from the event log (a period's first event is a full snapshot)
function filltimes(seq)
    d = Dict{Int,Float64}()
    for e in seq
        filled = e.filled isa AbstractVector{Bool} ? findall(e.filled) :
                                                     [u.index for u in e.filled if u.value]
        foreach(t -> get!(d, t, e.timestamp), filled)
    end
    return d
end
fb, fn = filltimes(base), filltimes(nets)
only_b, only_n = setdiff(keys(fb), keys(fn)), setdiff(keys(fn), keys(fb))
@printf("traps filled: base %d, networks %d; only in base %d, only in networks %d\n",
        length(fb), length(fn), length(only_b), length(only_n))
both = collect(intersect(keys(fb), keys(fn)))
dt   = [abs(fb[t] - fn[t]) for t in both]
worst = both[argmax(dt)]
@printf("fill times: max |difference| %.3g (trap %d, %.4f vs %.4f), traps differing > 1e-3: %d\n",
        maximum(dt), worst, fb[worst], fn[worst], count(>(1e-3), dt))

# where the differing traps are: in the network set, or elsewhere, and when they fill
dynset = Set(dyn)
covered(t) = t in dynset || any(a -> a in dynset, ts.supertraps_of[ts.regions[ts.footprints[t][1]]])
diff = sort([(abs(fb[t] - fn[t]), t) for t in both if abs(fb[t] - fn[t]) > 1e-3], rev = true)
@printf("differing traps: %d in or under a dyn_trap, %d elsewhere; earliest base fill time among them %.4f\n",
        count(x -> covered(x[2]), diff), count(x -> !covered(x[2]), diff), minimum(fb[t] for (_, t) in diff))
println("  trap  | base fill | nets fill | diff   | dyn? covered? parent?")
for (d, t) in diff[1:min(15, end)]
    @printf("  %5d | %9.4f | %9.4f | %6.3f | %s %s %s\n", t, fb[t], fn[t], d, t in dynset, covered(t),
            t > SWIM.numregions(ts))
end
early = sort(diff, by = x -> fb[x[2]])[1:min(10, end)]
println("earliest differing:")
for (d, t) in early
    @printf("  %5d | %9.4f | %9.4f | %6.3f | %s %s %s\n", t, fb[t], fn[t], d, t in dynset, covered(t),
            t > SWIM.numregions(ts))
end
# volumes at sampled times
tps = collect(0.0:10.0:300.0)
sb  = SWIM.all_states_at_timepoints(ts, base, tps; infiltration = infil, verbose = false)
sn  = SWIM.all_states_at_timepoints(ts, nets, tps; infiltration = infil, dyn_traps = dyn,
                                    verbose = false)
own(t) = ts.trapvolumes[t] - ts.subvolumes[t]
println("   t | stored base  stored nets  rel diff | max |dV|/capacity (trap) | filled differ | V < 0")
for (k, t) in enumerate(tps)
    a, b = sb[k].amounts, sn[k].amounts
    rel  = [own(i) > 0 ? abs(a[i] - b[i]) / own(i) : abs(a[i] - b[i]) for i in eachindex(a)]
    i    = argmax(rel)
    @printf("%4.0f | %11.2f  %11.2f  %8.1e | %8.1e (%5d)            | %5d         | %d\n",
            t, sum(a), sum(b), abs(sum(a) - sum(b)) / max(sum(a), eps()), rel[i], i,
            count(sb[k].filled .!= sn[k].filled), count(<(0), b))
end
