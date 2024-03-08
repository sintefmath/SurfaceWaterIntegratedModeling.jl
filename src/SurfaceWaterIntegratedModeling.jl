module SurfaceWaterIntegratedModeling

using LazyArtifacts
# ----------------------------------------------------------------------------
include("domain.jl")
include("spillfield.jl")
include("spillregions.jl")
include("spillpoints.jl")
include("trapvolumes.jl")
# ----------------------------------------------------------------------------
include("spillanalysis.jl")
include("TrapStructure.jl")
include("sshierarchy.jl")
# ----------------------------------------------------------------------------
include("fill_sequence/weatherevent.jl")
include("fill_sequence/spillevent.jl")
include("fill_sequence/spillgraph.jl")
include("fill_sequence/rateinfo.jl")
include("fill_sequence/flow.jl")
include("fill_sequence/fill_sequence.jl")
# ----------------------------------------------------------------------------
include("watercourses.jl")
include("utils.jl")
# ----------------------------------------------------------------------------
include("IOandplot.jl") 
# ----------------------------------------------------------------------------
include("artifacts.jl")

end # module SurfaceWaterIntegratedModeling
