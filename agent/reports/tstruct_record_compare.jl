# Compare two TrapStructure records (tstruct_record.jl) field by field, and their timings.
#   julia --project=. agent/reports/tstruct_record_compare.jl REF NEW
using SurfaceWaterIntegratedModeling, JLD2, Printf
ref, new = JLD2.load(ARGS[1], "results"), JLD2.load(ARGS[2], "results")
for name in sort(collect(keys(ref)))
    (t0, f0), (t1, f1) = ref[name], new[name]
    diff = sort([k for k in keys(f0) if !isequal(f0[k], f1[k])])
    @printf("%-22s %-40s %7.2f s -> %7.2f s  (%+5.1f %%)\n", name,
            isempty(diff) ? "identical" : "DIFFERENT: " * join(diff, ", "), t0, t1, 100 * (t1 / t0 - 1))
end
