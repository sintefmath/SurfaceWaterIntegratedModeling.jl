# Static topographic analysis

This page contains the documentation of the main functions related to the static analysis.

## TrapStructure

The `TrapStructure` contains all the fundamental results from running
[`spillanalysis`](@ref).  This includes the identified traps, their volumes,
footprints, how they are connected, and the associated watersheds and spill point
locations. `TrapStructure` also contains the *spillfield* (how flow is directed
between gridcells), as well as a copy of the topography, buildings and sinks
that went into the analysis.

```@autodocs
Modules = [SurfaceWaterIntegratedModeling]
Pages = ["TrapStructure.jl"]
```
## Spill analysis

The function `spillanalysis` runs a complete analysis on a provided terrain grid
and associated buildings and sinks, and returns a [`TrapStructure`](@ref) as
output.  This is the central function for static topographic analysis in SWIM.

```@autodocs
Modules = [SurfaceWaterIntegratedModeling]
Pages = ["spillanalysis.jl"]
```
## Lower-level functions

The functions listed below are used by the [`spillanalysis`](@ref)
function.  They can however also be run independently.

### Spillfield

The [`spillfield`](@ref) function computes the basic gravity-driven flow pattern
on a terrain.  For each cell, it determines the direction of flow and the
correpsonding neighbor downstream cell, using either a 4-element or 8-element stencil.

```@autodocs
Modules=[SurfaceWaterIntegratedModeling]
Pages = ["spillfield.jl"]
```

### Spillregions

For a given spillfield (as computed by the [`spillfield`](@ref) function), the
[`spillregions`](@ref) function separates the terain into individual watersheds,
each associated with a local minimum in the terrain surface. 

```@autodocs
Modules=[SurfaceWaterIntegratedModeling]
Pages = ["spillregions.jl"]
```

### Spillpoints

The [`Spillpoint`](@ref) struct contains the information for a particular spill
point, i.e. the point where a trap/lake spills over.  `Spillpoint` contains
the index of the downstream region, the location of the spillpoint itself, as
well as its vertical elevation.

The [`spillpoints`](@ref) function is used to compute all the spillpoints
associated with a set of *spill regions* (which has been previously computed
using [`spillregions`](@ref)).

```@autodocs
Modules=[SurfaceWaterIntegratedModeling]
Pages = ["spillpoints.jl"]
```
### Subtrap-supertrap hierarchy

The [`sshierarchy!`](@ref) identifies higher-level traps from a set of
subtraps. In other words, it identifies those traps that emerges when smaller
traps coalesce as they fill up.  It requires the output of
[`spillregions`](@ref) and [`spillpoints`](@ref) as part of it input.

```@autodocs
Modules=[SurfaceWaterIntegratedModeling]
Pages = ["sshierarchy.jl"]
```
### Trap volumes

The function [`trapvolumes`][@ref] is used to compute the volumes of all traps
and subtraps identified on an analysed surface. It requires the result of 
[`spillregions`](@ref) and [`spillpoints`](@ref) as input parameters.

```@autodocs
Modules=[SurfaceWaterIntegratedModeling]
Pages = ["trapvolumes.jl"]
```

### Domain

The `Domain2D` struct represent a rectangular subpart of a terrain grid.  It is
mainly used to specify subsets of a terrain that should be updated or processed
in parallel.

```@autodocs
Modules=[SurfaceWaterIntegratedModeling]
Pages = ["domain.jl"]
```

