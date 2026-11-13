# Runs every test file / testset touched by the NBS restructure, each isolated so a known
# failure (MaskingTests, bay-grid drift, raise_buildings!) cannot hide the others.
using Test, SurfaceWaterIntegratedModeling, LazyArtifacts, Random
const TD = joinpath(@__DIR__, "..", "..", "..", "test")
include(joinpath(TD, "basicTestFuns.jl"))
include(joinpath(TD, "watercourses_test.jl"))
files = ["culvert_regions_test.jl", "dynamics_test.jl", "dynamic_membership_test.jl",
         "network_driver_test.jl", "nbs_routing_test.jl", "network_states_test.jl",
         "network_watercourses_test.jl", "fill_sequence_rates_test.jl"]
fails = String[]
for f in files
    try
        include(joinpath(TD, f))
    catch e
        push!(fails, f)
    end
end
# testsets from runtests.jl that exercise fill_sequence (skip the known-broken ones)
lines = readlines(joinpath(TD, "runtests.jl"))
grids = join(lines[15:24], "\n")
include_string(Main, grids, joinpath(TD, "runtests.jl"))
for name in ("Upstream area", "Sequencing")
    a = findfirst(l -> occursin("@testset \"$name\"", l), lines)
    b = a + findfirst(l -> l == "end", lines[a+1:end])
    try
        include_string(Main, join(lines[a:b], "\n"), joinpath(TD, "runtests.jl"))
    catch e
        push!(fails, name)
    end
end
println("\nRUN_TESTS FAILED: ", isempty(fails) ? "none" : join(fails, ", "))
