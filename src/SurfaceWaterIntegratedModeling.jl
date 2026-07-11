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
# NBS storage model + Dyn* element types: defined before the fill_sequence /
# watercourses layers because those reference DynNBSPlacement in their signatures.
# nbs_elements (NBSSystem/NBSLayer) must precede elements (DynNBSPlacement.system).
include("dynamics/nbs_elements.jl")
include("dynamics/elements.jl")
# ----------------------------------------------------------------------------
include("fill_sequence/weatherevent.jl")
include("fill_sequence/spillevent.jl")
include("fill_sequence/spillgraph.jl")
include("fill_sequence/rateinfo.jl")
include("fill_sequence/flow.jl")
include("fill_sequence/fill_sequence.jl")
# ----------------------------------------------------------------------------
include("watercourses.jl")
include("dynamics/culvert_rate.jl")
include("utils.jl")
include("dynamics/networksolver.jl")
include("dynamics/network_context.jl")
include("dynamics/network_utils.jl")
include("dynamics/network_reachability.jl")
# ----------------------------------------------------------------------------
include("IOandplot.jl")
# ----------------------------------------------------------------------------
include("artifacts.jl")

end # module SurfaceWaterIntegratedModeling
