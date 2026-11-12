using SurfaceWaterIntegratedModeling, LazyArtifacts, JLD2
const SWIM = SurfaceWaterIntegratedModeling
d = joinpath(artifact"swim_testdata", "data")
for (name, f) in (("mini", "small/mini.txt"), ("bay", "small/bay.txt"), ("synsurf", "synthetic/synsurf.txt"))
    g = SWIM.loadgrid(joinpath(d, f))
    ts = spillanalysis(g)
    z = findall(ts.trapvolumes .== 0.0)
    println(name, ": traps=", numtraps(ts), " zero-volume=", length(z), " (leaf: ", count(<=(numregions(ts)), z), ")")
end
surb = get(ENV, "SURBAREA_DIR", normpath(joinpath(@__DIR__, "..", "..", "..", "..", "SUrbArea")))
dtm = JLD2.load(joinpath(surb, "debug", "synthetic_dtm_flattened.jld2"), "dtm")
ts = spillanalysis(dtm)
println("SUrbArea synthetic: zero-volume=", count(ts.trapvolumes .== 0.0))
