module SurfaceWaterIntegratedModeling

using LazyArtifacts
# ----------------------------------------------------------------------------
include("domain.jl")
include("spillfield.jl")
include("spillregions.jl")
include("spillpoints.jl")
include("trapvolumes.jl")
# ----------------------------------------------------------------------------
# NBS element type definitions are pulled in ahead of TrapStructure so the latter
# can name `NBSPlacement` in a field (see `nbs` field).  These files define only
# types/functions (no include-time reference to TrapStructure), so loading them
# here is safe.
include("dynamics/elements.jl")
include("dynamics/nbs_elements.jl")
include("nbs_placement.jl")
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
include("dynamics/culvert_rate.jl")
include("utils.jl")
include("dynamics/networksolver.jl")
include("dynamics/network_context.jl")
# ----------------------------------------------------------------------------
include("IOandplot.jl")
# ----------------------------------------------------------------------------
include("artifacts.jl")

end # module SurfaceWaterIntegratedModeling
