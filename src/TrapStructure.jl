import Graphs

export TrapStructure, numregions, numtraps, subtrapsof, parentof


"""
    struct TrapStructure{T<:Real}

A struct representing a watershed drainage trap for topographical analysis.

# Fieldnames

- `topography::Matrix{T}`: raster grid of terrain height values.
- `flowgraph::Graphs.SimpleDiGraph`: graph describing flow routing between cells
- `regions::Matrix{Int}`: raster grid with region numbers (see [`spillregions`](@ref))
- `spillpoints::Vector{Spillpoint{T}}`: vector with spill point information per trap
                                     (see [`spillpoints`](@ref))
- `trapvolumes::Vector{T}`: computed trap volumes (see [`trapvolumes`](@ref))
- `subvolumes::Vector{T}`: the part of each trap's volume that is fully contained within
                           its subtraps
- `footprints::Vector{Vector{Int}}`: one entry per trap, listing the (linear) indices
                                     of all its cells 
- `lowest_subtraps_for::Vector{Vector{Int}}`: one entry per trap, listing all its
                                              lowest-level subtraps
- `supertraps_of::Vector{Vector{Int}}`: one entry per lowest-level trap, listing all
                                        the supertraps it belongs to (including itself).
- `agglomerations::Graphs.SimpleDiGraph`: hierarchy of sub/supertraps, presented as a
                                          graph structure
- `clip_mask::Union{Matrix{Bool}, Nothing}`: clip mask (0: active terrain, 1: clipped away)
- `sinks::Vector{CartesianIndex{2}}`: list of sinks, in term of grid
                                      cell coordinates
- `cut_edges::Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}}`:
                                              Dict of cut edges (barriers), i.e. edges
                                              that block flow between adjacent cells
"""
mutable struct TrapStructure{T<:Real}

    topography::Matrix{T}           # raster grid of terrain height values
    flowgraph::Graphs.SimpleDiGraph  # graph describing flow routing between cells
    trap_bottoms::Vector{CartesianIndex{2}} # list of trap bottom cells in grid
    regions::Matrix{Int}            # raster grid with region numbers
    spillpoints::Vector{Spillpoint{T}} # vector with spill point information per trap
    trapvolumes::Vector{T}          # computed trap volumes
    subvolumes::Vector{T}           # the part of a trap's volume that is fully contained
                                    # within its subtraps
    footprints::Vector{Vector{Int}} # one entry per trap, giving all its cells
    lowest_subtraps_for::Vector{Vector{Int}} # one entry per trap, giving all
                                             # its lowest-level subtraps
    supertraps_of::Vector{Vector{Int}} # one entry per _lowest-level_ trap,
                                       # listing all the supertraps it belongs
                                       # to (including itself)
    agglomerations::Graphs.SimpleDiGraph # hierarchy of sub/supertraps
    clip_mask::Union{Matrix{Bool}, Nothing} # clip mask (0: active terrain, 1: clipped away)
    sinks::Vector{CartesianIndex{2}}      # list cells being sinks
    waterbodies::Vector{CartesianIndex{2}} # list of cells belonging to waterbodies
    cut_edges::Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}} # Dict of cut edges (barriers)
end

# ----------------------------------------------------------------------------
"""
    numregions(tstruct)

Returns the number of regions identified in the [`TrapStructure`](@ref) `tstruct`.
The number of regions is identical to the number of lowest-level traps.
"""
function numregions(tstruct::TrapStructure)
    # there should be exactly one entry in 'supertraps_of' for each lowest-level
    # trap, i.e. each region
    return length(tstruct.supertraps_of)
end

# ----------------------------------------------------------------------------
"""
    numtraps(tstruct)

Returns the number of traps identified in the [`TrapStructure`](@ref) `tstruct`.
This includes the lowest-level traps as well as all higher-level traps.
"""
function numtraps(tstruct::TrapStructure)
    # There should be exactly one spillpoint per trap
    return length(tstruct.spillpoints)
end

# ----------------------------------------------------------------------------
"""
    subtrapsof(tstruct, trap_ix)

Given a [`TrapStructure`](@ref) tstruct, returns the indices to the immediate
subtraps of the trap indexed `trap_ix`. This is usually zero or two, but can be
higher than two in degenerate cases.
"""
function subtrapsof(tstruct::TrapStructure, trap_ix::Int)
    return Graphs.inneighbors(tstruct.agglomerations, trap_ix)
end

# ----------------------------------------------------------------------------
"""
    parentof(tstruct, trap_ix)

Given a [`TrapStructure`](@ref) tstruct, returns the index of the parent trap
(immediate higher-level trap) of the trap indexed `trap_ix`.  If this trap does 
not have a parent trap, nothing is returned.
"""
function parentof(tstruct::TrapStructure, trap_ix::Int)
    out = Graphs.outneighbors(tstruct.agglomerations, trap_ix)
    @assert length(out) <= 1
    return isempty(out) ? nothing : out[1]
end

# ----------------------------------------------------------------------------
function regioncells(regions::Matrix{Int})
    # returns a vector of vectors, where the i'th entry is a vector of the
    # linear indices of all cells in region i.  Does not consider regions with
    # zero or negative numbers (out-of-domain regions).
    numregions = maximum(regions)
    cells = [Int[] for i in 1:numregions]

    # iterate through cells in 'regions' and collect the indices of cells for each region
    for idx in eachindex(regions)
        region_ix = regions[idx]
        if region_ix > 0
            push!(cells[region_ix], idx)
        end
    end
    
    return cells
end
