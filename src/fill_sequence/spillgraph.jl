import Graphs
export SpillGraph, compute_complete_spillgraph, update_spillgraph!


"""
    mutable struct SpillGraph

A struct representing a tree-like spillgraph.  It specifies how water flow from
one (filled) trap to the next.  When a trap is filled, it may spill into its
downstream trap, or into its parent trap (depending on whether it shares a parent
with its downstream trap, and the downstream trap is filled).

The spill graph is recorded as a set of edges, stored in the Dict `edges`, so that
if trap 'i' spills into trap 'j', then `edges[i] = k`.
"""
mutable struct SpillGraph
    edges::Dict{Int, Int} # Trap indices.  The first 'spills into' the second.
end

"""
Default constructor for SpillGraph creates an empty spillgraph.
"""
function SpillGraph()
    return SpillGraph(Dict{Int, Int}())
end

"""
    get_graph(spill_graph)

Converts the [`SpillGraph`](@ref) `spill_graph` into a `Graphs.SimpleDiGraph`
from the `Graphs` package (https://juliagraphs.org).
"""
function get_graph(sg::SpillGraph)
    edges = Vector{Graphs.SimpleEdge{Int}}(
        [Graphs.SimpleEdge{Int}(k, v) for (k, v) ∈ sg.edges])
    return Graphs.SimpleDiGraph(edges)
end

function _set_outedge!(sg::SpillGraph, from::Int, target::Int)
    old_target = get(sg.edges, from, 0)
    sg.edges[from] = target
    # return values: (start node, old endnode, new endnode)
    return (from, old_target, target)
end

function _remove_outedge!(sg::SpillGraph, from::Int)

    @assert from ∈ keys(sg.edges) "Tried to remove inexistant edge"
    if from ∈ keys(sg.edges)
        old_target = pop!(sg.edges, from, 0)
        # return values: (start node, old endnode, new endnode)
        return (from, old_target, 0) 
    end
end

function _target_of(sg::SpillGraph, ix::Int)
    return get(sg.edges, ix, nothing)
end


# ----------------------------------------------------------------------------
# the following function is slow but should be failsafe
"""
    compute_complete_spillgraph(tstruct, full_traps)

Compute the complete spillgraph, given a trapping structure and a (sub)set of
filled traps.

Returns an object of type [`SpillGraph`](@ref).  

The `full_traps` vector must be *consistent*, i.e. no filled trap can have
non-filled subtraps.

# Arguments
- `tstruct::TrapStructure{T<:Real}`: trapping structure
- `full_traps::Vector{Bool}`: a vector with one entry per trap, specifying if that
                              trap has been completely filled or not.

See also [`TrapStructure`](@ref), [`update_spillgraph!`](@ref).
"""
function compute_complete_spillgraph(tstruct::TrapStructure{T},
                                     full_traps::Vector{Bool}) where T <: Real

    num_traps = numtraps(tstruct)
    @assert num_traps == length(full_traps)

    # defining two simple shorthands
    children_of = x -> Graphs.inneighbors(tstruct.agglomerations, x)
    parent_of = x -> Graphs.outneighbors(tstruct.agglomerations, x)
    
    # identify the traps that are filled 
    @assert _valid_trap_status(full_traps, tstruct)

    # identify traps that are active 
    active_traps = [full_traps[x] | all(full_traps[children_of(x)]) for x in 1:num_traps]

    # determine where the filled traps spill
    full_trap_ixs = findall(full_traps)

    function target(x)
        parent = parent_of(x)
        ds_reg = tstruct.spillpoints[x].downstream_region
        if !isempty(parent) && active_traps[parent[1]]
            return parent[1]
        end
        return ds_reg > 0 ? ds_reg : num_traps + 1 # latter value flags out-of-domain
    end

    spillgraph = SpillGraph()
    for x in full_trap_ixs
        _set_outedge!(spillgraph, x, target(x))
    end

    return spillgraph
end

# ----------------------------------------------------------------------------
function _valid_trap_status(filled_traps, tstruct)
    # check that no unfilled trap has a filled supertrap
    filled_ixs = findall(filled_traps)
    for trap in filled_ixs
        if any(.!filled_traps[Graphs.inneighbors(tstruct.agglomerations, trap)])
            # found filled trap with non-filled child
            return false
        end
    end
    return true
end

# ----------------------------------------------------------------------------
"""
    update_spillgraph!(spillgraph, fill_changes, tstruct)

Update an existing spillgraph after there have been changes in which traps are
filled or not.

The argument `spillgraph` will be modified.

# Arguments
- `spillgraph::SpillGraph`: the spillgraph to be updated
- `fill_changes::Vector{IncrementalUpdate{Bool}}`: changes in the fill status of 
       one or more traps, presented as a vector of incremental updates
- `tstruct::TrapStructure{T<:Real}`: the [`TrapStructure`](@ref) representing the
       trapping structure of the terrain the spill graph is associated with

See also [`IncrementalUpdate`](@ref), [`SpillGraph`](@ref) and 
[`compute_complete_spillgraph`](@ref).
"""
function update_spillgraph!(spillgraph::SpillGraph,
                             fill_changes::Vector{IncrementalUpdate{Bool}},
                             tstruct::TrapStructure{T}) where T <: Real

    num_traps = numtraps(tstruct)
    changes = Vector{Tuple{Int, Int, Int}}() # source, old target, new target

    for fc in fill_changes
        if fc.value # filled trap
            dsreg = tstruct.spillpoints[fc.index].downstream_region
            # num_traps+1 represents 'out of domain'
            target = dsreg > 0 ? dsreg : num_traps+1
            push!(changes, _set_outedge!(spillgraph, fc.index, target))
        else # un-filled trap
            push!(changes, _remove_outedge!(spillgraph, fc.index))
        end
    end
    
    # Now that all updates have been initially applied, check for sibling cycles
    # that must be updated.  (We do this in a second pass, since it is not
    # guaranteed that the `fill_changes` have been presented in chronological order).
    for fc in fill_changes
        
        # identify siblings (if any)
        siblings, parent = _siblings_of(tstruct, fc.index)
        
        isempty(siblings) && continue # if the trap has no sibling, nothing to do

        if fc.value
            # This trap just filled.  Check if it completed a 'cycle' so that
            # flow should be redirected to its parent
            if all([_isfull(x, spillgraph) for x in siblings])
                # all siblings are full.  Ensure they are all spilling to parent
                for s in siblings
                    push!(changes, _set_outedge!(spillgraph, s, parent))
                end
            end
        else
            # This trap just stopped being filled.  Redirect any flow from other
            # siblings that still flow into parent.
            for s in siblings
                if s ∈ keys(spillgraph.edges)
                    @assert spillgraph.edges[s] == parent
                    dsreg = tstruct.spillpoints[s].downstream_region
                    target = dsreg > 0 ? dsreg : num_traps + 1
                    push!(changes, _set_outedge!(spillgraph, s, target))
                end
            end
        end
    end

    # Make the list of incremental updates.  If there are multiple changes
    # registered for a trap, merge them to ensure the old and new targets are
    # consistent with the old and new graph
    tmp = Dict{Int, Tuple{Int, Int}}()
    for ch in changes
        existing = get(tmp, ch[1], nothing)
        if isnothing(existing)
            tmp[ch[1]] = (ch[2], ch[3]) # from/to
        else
            tmp[ch[1]] = (existing[1], ch[3]) # re-use the old from value
        end
    end
    return [IncrementalUpdate(k, v) for (k, v) ∈ tmp]
end

# ----------------------------------------------------------------------------
function _siblings_of(tstruct, trap)
    parent = Graphs.outneighbors(tstruct.agglomerations, trap)
    if isempty(parent)
        return [], nothing
    end
    @assert(length(parent) == 1)
    parent = parent[1]
    siblings = Graphs.inneighbors(tstruct.agglomerations, parent)
    return siblings, parent
end

# ----------------------------------------------------------------------------
function _isfull(trap, spillgraph)
    # a trap is full if it spills to another trap (or to outside of domain)
    return !isnothing(_target_of(spillgraph, trap))
end

