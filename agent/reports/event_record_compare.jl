# Compare two event records (event_record.jl): bit-identical event histories, and timings.
#   julia --project=. agent/reports/event_record_compare.jl REF NEW
using JLD2, Printf
ref, new = JLD2.load(ARGS[1], "results"), JLD2.load(ARGS[2], "results")
for name in sort(collect(keys(ref)))
    (t0, r0), (t1, r1) = ref[name], new[name]
    @printf("%-28s %-10s %7.2f s -> %7.2f s  (%+5.1f %%)\n", name,
            r0 == r1 ? "identical" : "DIFFERENT", t0, t1, 100 * (t1 / t0 - 1))
end
