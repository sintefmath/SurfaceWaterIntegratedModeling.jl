import Graphs
export compute_flow

"""
    compute_flow(spillgraph, precipitation, infiltration, tstruct, verbose=false)

For a given spillgraph, precipitation and infiltration rates, compute the runoff
and trap inflows, in the form of a [`RateInfo`](@ref) object.

# Arguments
- `spillgraph::SpillGraph`: the current spillgraph, i.e. a tree graph representing 
                            which traps spill into which
- `precipitation::Union{Real, Matrix{<:Real}}`: specifies the precipitation rate. 
        Can be given either as a single scalar, or cell-wise in the form of a Matrix
- `infiltration::Union{Real, Matrix{<:Real}}`: specifies the infiltration rate.
        Can be given either as a single scalar, or cell-wise in the form of a Matrix
- `tstruct::TrapStructure{<:Real}`: object representing the trap sturcture
- `verbose::Bool`: if `true`, dump progress information during computation

See also [`SpillGraph`](@ref), [`TrapStructure`](@ref).
"""
function compute_flow(spillgraph::SpillGraph,
                      precipitation::Union{Real, Matrix{<:Real}},
                      infiltration::Union{Real, Matrix{<:Real}},
                      tstruct::TrapStructure{<:Real},
                      verbose::Bool=false)
    num_traps = length(tstruct.spillpoints)
    num_regions = length(tstruct.supertraps_of)

    # compute initial spillfield with all traps empty
    println("compute initial rateinfo with all traps empty")
    rateinfo = _compute_initial_rateinfo(precipitation, infiltration, tstruct)

    # --- Add influence of traps spilling over ---

    # make a sorted spill graph, so that traversing it in order will not
    # invalidate the result
    g = get_graph(spillgraph)
    sortedg = Graphs.topological_sort_by_dfs(g)

    icount = 0 # @@@
    for sourcenode in sortedg
        icount = icount + 1
        verbose && (mod(icount, 100) == 0) && println("compute_flow iteration: ", icount)

        # target may either be downstream trap, parent trap or out of domain
        targetnode = Graphs.outneighbors(g, sourcenode)
        if isempty(targetnode)
            continue
        end
        targetnode = targetnode[1]

        if _is_parent(targetnode, sourcenode, tstruct)
            # target is parent, no flow tracking over terrain necessary
            setinflow!(rateinfo,
                       targetnode,
                       getinflow(rateinfo, targetnode) + getinflow(rateinfo, sourcenode))
        else
            # track flow downstream until trap, sink or domain boundary is encountered
            signed_flow = getinflow(rateinfo, sourcenode) - getsmax(rateinfo, sourcenode)
            outflow = max(signed_flow, 0.0) # outflow is always positive
            # changes the runoff values of the affected cells, and the inflow and Smin/Smax of affected traps
            _track_flow!(rateinfo, sourcenode, outflow, tstruct) # update 'rateinfo'
        end
    end
    return rateinfo
end

# ----------------------------------------------------------------------------
"""
    _is_parent(target, source, tstruct) -> Bool

True when `target` is `source`'s immediate parent in the agglomeration hierarchy, i.e. `source`
spills by merging into the supertrap that subsumes it rather than by running over terrain.  That
distinction decides whether flow is handed straight over or tracked downstream.  Nothing is
mutated.
"""
function _is_parent(target, source, tstruct)
    parent = Graphs.outneighbors(tstruct.agglomerations, source)
    @assert isempty(parent) || length(parent) == 1
    return !isempty(parent) && parent[1] == target
end

# ----------------------------------------------------------------------------
"""
    _is_trap_bottom(cell, tstruct) -> Bool

True when `cell` is the bottom of a trap, i.e. where arriving flow accumulates.

A cell with no downstream neighbour that is not a sink and belongs to a non-negative region must
be a trap bottom: a cell spilling out of the domain would belong to a negative region instead.
Nothing is mutated.
"""
function _is_trap_bottom(cell, tstruct)
    @assert cell > 0
    return isempty(Graphs.outneighbors(tstruct.flowgraph, cell)) && 
        (tstruct.regions[cell] > 0) &&
        !_is_sink(cell, tstruct)
end

# # ----------------------------------------------------------------------------
# function _is_trap_bottom(cell, tstruct)
#     #return cell > 0 && tstruct.spillfield[cell] == -1
#     return cell > 0 && isempty(outneighbors(tstruct.flowgraph, cell))
# end

# ----------------------------------------------------------------------------
"""
    _is_sink(cell, tstruct) -> Bool

True when `cell` is one of the terrain's designated sinks, where water leaves the model rather
than accumulating.  Always `false` when the structure defines no sinks.  Nothing is mutated.
"""
function _is_sink(cell, tstruct)
    # @@@ TODO: should this be optimized, to avoid creation of intermediary
    # CartesianIndices objects every time?
    return cell > 0 && (tstruct.sinks !== nothing) &&
        CartesianIndices(size(tstruct.topography))[cell] in tstruct.sinks
end
# ----------------------------------------------------------------------------
"""
    _track_flow!(rateinfo, node, amount, tstruct) -> Float64

Propagate a signed `amount` of flow from trap `node`'s spillpoint downstream over the terrain,
adding it to each cell's runoff until it is exhausted or the flow reaches a trap bottom, a sink,
or the domain boundary.  Returns the amount still carried on arrival (`0.0` if it ran out).

Each cell's runoff changes by the amount arriving, and what passes on is
`max(R + amount, 0) − max(R, 0)` (`R` the cell's runoff before): added flow first fills a cell's
remaining infiltration capacity (negative runoff), and removed flow takes away at most what the
cell was passing on.  The same rule in both directions is what makes a removal along the same
cells exactly undo an earlier addition, as `_update_flow!` relies on when it retracts a spill.
The network solver's `_attenuate_range` applies the same rule to the grid this produces.

**Mutates `rateinfo`**: the runoff of every cell traversed, the `Smin`/`Smax` of the supertraps
over the arrival region, and — when the flow lands on a trap bottom — that region's inflow.
"""
function _track_flow!(rateinfo, node, amount, tstruct)

    spoint = tstruct.spillpoints[node]
    
    # first, check that the spillpoint is not directly out of the domain, in which
    # case we return immediately
    if spoint.current_region_cell == spoint.downstream_region_cell
        return 0.0
    end

    cell = spoint.downstream_region_cell

    while cell > 0
        prev_value = getrunoff(rateinfo, cell)
        setrunoff!(rateinfo, cell, prev_value + amount)
        amount = max(prev_value + amount, 0.0) - max(prev_value, 0.0)
        amount == 0.0 && break # nothing left to propagate/remove
        ds = Graphs.outneighbors(tstruct.flowgraph, cell)
        isempty(ds) && break # no downstream cell, stop here (either trap bottom, sink, or 
                             # domain boundary)
        @assert length(ds) == 1
        cell = ds[1]
    end

    # recompute inflow and footprint infiltration of affected supertraps
    @assert cell > 0
    reg = tstruct.regions[cell]
    if reg > 0
        _update_Smin_Smax!(rateinfo, tstruct, tstruct.supertraps_of[reg])
    end

    if reg > 0 && _is_trap_bottom(cell, tstruct)
        setinflow!(rateinfo, reg, getinflow(rateinfo, reg) + amount)
    end
    
    return amount
end

# ----------------------------------------------------------------------------
"""
    _compute_initial_rateinfo(precipitation, infiltration, tstruct) -> RateInfo

The [`RateInfo`](@ref) for a terrain with every trap empty: the baseline runoff field and each
lowest-level trap's inflow, before any trap spills.  Scalar `precipitation` / `infiltration` are
broadcast to the grid.

Higher-level traps start at zero inflow — they only receive once their children fill.  The
runoff field is NBS-oblivious; an NBS's effect is applied later as a signed correction inside
the network solver.  Nothing is mutated; the `RateInfo` is freshly built.
"""
function _compute_initial_rateinfo(precipitation, infiltration, tstruct)
    if typeof(precipitation) <: Real
        precipitation = precipitation .* ones(size(tstruct.regions))
    end
    if typeof(infiltration) <: Real
        infiltration = infiltration .* ones(size(tstruct.regions))
    end

    # compute the basic runoff field with all traps empty (NBS-oblivious; the NBS effect
    # is applied later as a signed correction inside the network solver)
    runoff, reg_accum, _, _ =
        watercourses(tstruct, [false],
                     precipitation=precipitation,
                     infiltration=infiltration)

    Smin = zeros(length(tstruct.footprints))
    Smax = zeros(length(tstruct.footprints))

    # all lowest-level traps are assigned the inflow from their region,
    # higher-level traps left at zero for now.
    trap_inflow = vcat(reg_accum, zeros(numtraps(tstruct) - numregions(tstruct)))
    rateinfo = RateInfo(runoff, Smax, Smin, trap_inflow)
    _update_Smin_Smax!(rateinfo, tstruct, 1:length(tstruct.footprints))

    return rateinfo
end

# ----------------------------------------------------------------------------
"""
    _ponding_infiltration(rateinfo, tstruct, trap) -> Float64

Maximum remaining infiltration capacity of `trap`: its remaining per-cell capacity summed over
the footprint, **excluding** cells whose terrain lies at or above the trap's spillpoint.  This
is the trap's `Smax`.  Nothing is mutated.

An excluded cell never holds standing water while part of the trap — water reaching the spill
level flows out rather than pooling — so it carries no trap infiltration.  Excluding it keeps
the full-trap loss continuous at capacity and consistent with the dynamic network solver,
removing the fill/unfill chatter the discontinuity used to cause.  `fill_trap_until` and
`_build_trap_geometry` apply the same `topography >= spillpoint` rule to their own footprints.

The test is on the cell's actual terrain height; there is no need to raise the bottom to the
subtrap spillpoint as the volume/level code does, because a child's spillpoint is always `<=`
the parent's, so `max(topo, child_sp) < sp` reduces to `topo < sp` — and where they would
differ (a degenerate zero-own-volume parent) the plain test is the correct one.  The loop is
explicit and allocation-free: it runs per trap and is refreshed per event, so it must stay cheap
on large terrains.
"""
function _ponding_infiltration(rateinfo, tstruct, trap)
    sp = Float64(tstruct.spillpoints[trap].elevation)   # concrete: Spillpoint.elevation is ::Real
    s  = 0.0
    @inbounds for c in tstruct.footprints[trap]
        Float64(tstruct.topography[c]) < sp && (s -= min(getrunoff(rateinfo, c), 0.0))
    end
    return s
end

"""
    _update_Smin_Smax!(rateinfo, tstruct, traps) -> nothing

Recompute the infiltration bounds of each trap in `traps`, **mutating `rateinfo`**.

`Smax` is the trap's full remaining capacity ([`_ponding_infiltration`](@ref)), incurred once it
is full and its whole footprint is wetted.  `Smin` is the least it can lose — for a parent, the
`Smax` of its immediate children, whose footprints are already submerged; zero for a leaf.  A
trap's net inflow therefore lies between `inflow - Smax` and `inflow - Smin`, which is what makes
the fill-sequence changetime estimate a bracket rather than an exact time.

The two passes cannot be merged: every `Smax` must be final before any `Smin` reads its
children's.
"""
function _update_Smin_Smax!(rateinfo, tstruct, traps)

    for i in traps
        setsmax!(rateinfo, i, _ponding_infiltration(rateinfo, tstruct, i))
    end

    for i in traps
        setsmin!(rateinfo, i, sum(getsmax(rateinfo, subtrapsof(tstruct, i))))
    end
end

# ----------------------------------------------------------------------------
"""
    _update_flow!(rateinfo, graph_updates, tstruct, sgraph) -> nothing

Bring `rateinfo` up to date after the spillgraph changed, given `graph_updates` as one
`(trap, (from, to))` diff per redirected edge.  **Mutates `rateinfo`**; `sgraph` is read as the
already-updated graph.

Runs in two passes, and the order is load-bearing: every old flow is retracted before any new
flow is added, so a trap whose inflow feeds a redirected edge is not counted under both routings.
The retraction re-derives the old routing from the updates' `from` fields, since `sgraph` already
holds the new one.  As each edge's outflow is fully deducted its entry is zeroed, stopping
further upstream deduction across it.
"""
function _update_flow!(rateinfo, graph_updates, tstruct, sgraph)

    # When redirecting flow, we need to know the spill graph before updates happened.
    # We deduce this from the updates, and use it with _propagate_amount! below
    # when correcting flows
    old_edges = Dict([(update.index, update.value[1]) for update in graph_updates])

    # subtract flow from previous target (this must be completed before we start
    # adding new flows)
    for gup in graph_updates
        trap, from, to = gup.index, gup.value[1], gup.value[2]
        signed_outflow = getinflow(rateinfo, trap) - getsmax(rateinfo, trap)

        if from > 0
            outflow = _is_parent(from, trap, tstruct) ? getinflow(rateinfo, trap) :
                                                        max(signed_outflow, 0.0)
            _propagate_amount!(rateinfo, trap, -1 * outflow, tstruct,
                               sgraph, old_edges=old_edges)
            if haskey(old_edges, trap)
                # total outflow across this edge has now been deducted, any further
                # deduction from upstream traps stop here
                old_edges[trap] = 0
            end
        end
    end

    # add flow to new target
    for gup in graph_updates
        trap, from, to = gup.index, gup.value[1], gup.value[2]
        signed_outflow = getinflow(rateinfo, trap) - getsmax(rateinfo, trap)
        
        if to > 0
            outflow = _is_parent(to, trap, tstruct) ? getinflow(rateinfo, trap) :
                                                      max(signed_outflow, 0.0)
            _propagate_amount!(rateinfo, trap, outflow, tstruct, sgraph)
        end
    end
end

# ----------------------------------------------------------------------------
"""
    _propagate_amount!(rateinfo, trap, outflow, tstruct, sgraph; old_edges=nothing) -> nothing

Push a signed `outflow` from `trap` along the spill chain, cascading downstream for as long as
each target is itself full, and stopping at the first non-filled trap or the domain boundary.

Each hop either hands the amount straight to a parent trap (a merge, no terrain involved) or
tracks it over terrain via [`_track_flow!`](@ref), which may exhaust it en route — once it
reaches zero the walk ends.  A negative `outflow` retracts a previous spill along the same
chain.

**Mutates `rateinfo`** (runoff, trap inflows, and `Smin`/`Smax` via `_track_flow!`).  Pass
`old_edges` to walk the *pre-update* spillgraph, which is what makes retraction follow the route
the water actually took; it falls back to `sgraph` for traps it does not name.
"""
function _propagate_amount!(rateinfo, trap, outflow, tstruct, sgraph; old_edges=nothing)

    while outflow != 0.0

        target = isnothing(old_edges) ? 
            get(sgraph.edges, trap, 0) :
            get(old_edges, trap, get(sgraph.edges, trap, 0))

        (target == 0) && break # we reached a non-filled trap.  Stop propagating here

        if _is_parent(target, trap, tstruct)
            setinflow!(rateinfo, target, getinflow(rateinfo, target) + outflow)
            trap = target 
        else
            outflow = _track_flow!(rateinfo, trap, outflow, tstruct)
            (target > numtraps(tstruct)) && break # exit if we reach domain boundary

            trap = target
        end
    end
end
