using Documenter, Literate, SurfaceWaterIntegratedModeling
using LazyArtifacts

push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "../examples/")

# Prepare example scripts
function set_to_cairo(content)
    content = replace(content, "GLMakie" => "CairoMakie")
    return content
end

Literate.markdown("examples/urban.jl", "docs/src/"; execute=false, preprocess = set_to_cairo)
Literate.markdown("examples/flat_areas.jl", "docs/src/"; execute=false, preprocess = set_to_cairo)
Literate.markdown("examples/synthetic.jl", "docs/src/"; execute=false, preprocess = set_to_cairo)

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
        "Examples" => ["urban.md", "synthetic.md", "flat_areas.md"],
        "Static analysis" => "static.md",
        "Dynamic analysis" => "dynamic.md",
        "Utilities and visualization" => "utils.md",
        "Index" => "indexlist.md"
    ]
)

deploydocs(;
    repo="github.com/sintefmath/SurfaceWaterIntegratedModeling.jl",
    devbranch="main",
)
