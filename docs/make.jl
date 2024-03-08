using Documenter, Literate, SurfaceWaterIntegratedModeling
using LazyArtifacts

push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "../examples/")

# Prepare example scripts
Literate.markdown("../examples/urban.jl", "src/"; execute=false)
Literate.markdown("../examples/flat_areas.jl", "src/"; execute=false)
Literate.markdown("../examples/synthetic.jl", "src/"; execute=false)

# Build documentation
makedocs(
    modules = [SurfaceWaterIntegratedModeling], 
    sitename = "SWIM",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    draft=false,
    pages = [
        "Introduction" => "index.md",
        "Static analysis" => "static.md",
        "Dynamic analysis" => "dynamic.md",
        "Utilities and visualization" => "utils.md",
        "Examples" => ["urban.md", "synthetic.md", "flat_areas.md"],
        "Index" => "indexlist.md"
    ]
)