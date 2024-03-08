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
    rateinfo = _compute_initial_rateinfo(precipitation, infiltration, tstruct)

    # --- Add influence of traps spilling over ---

    # make a sorted spill graph, so that traversing it in order will not
    # invalidate the result
    g = get_graph(spillgraph)
    sortedg = Graphs.topological_sort_by_dfs(g)

    icount = 0 # @@
    for sourcenode in sortedg
        icount = icount + 1
        verbose && (mod(icount, 100) == 0) && println("compute_flow iteration: ", icount)

        # target may either be downstream trap, parent trap or out of domain
        targetnode = Graphs.outneighbors(g, sourcenode)
        if isempty(targetnode)
            continue
        end
        targetnode = targetnode[1]

        signed_flow = getinflow(rateinfo, sourcenode) - getsmax(rateinfo, sourcenode)
        if _is_parent(targetnode, sourcenode, tstruct)
            # target is parent, no flow tracking over terrain necessary
            setinflow!(rateinfo,
                       targetnode,
                       getinflow(rateinfo, targetnode) + signed_flow)
        else
            # track flow downstream until trap, sink or domain boundary is encountered
            outflow = max(signed_flow, 0.0) # outflow is always positive
            _track_flow!(rateinfo, sourcenode, outflow, tstruct) # update 'rateinfo'
        end
    end
    return rateinfo
end

# ----------------------------------------------------------------------------
function _is_parent(target, source, tstruct)
    parent = Graphs.outneighbors(tstruct.agglomerations, source)
    @assert isempty(parent) || length(parent) == 1
    return !isempty(parent) && parent[1] == target
end

# ----------------------------------------------------------------------------
function _is_trap_bottom(cell, tstruct)
    return cell > 0 && tstruct.spillfield[cell] == -1
end
# ----------------------------------------------------------------------------
function _is_sink(cell, tstruct)
    return cell > 0 && tstruct.spillfield[cell] == -3
end
# ----------------------------------------------------------------------------
function _track_flow!(rateinfo, node, amount, tstruct)

    initial_sign = sign(amount)
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

        if initial_sign == 1 # (we are adding flow)
            # In case there was remaining infiltration, deduct this from the amount
            amount += min(prev_value, 0.0) # gradually reduce towards zero
        else # (we are subtracting flow
            amount -= min(prev_value, 0.0) # gradually increase towards zero
        end

        (sign(amount) != initial_sign) && break # stop when there is no more to
                                                # propagate/remove
        _is_trap_bottom(cell, tstruct) && break # return if the cell is a trap bottom

        cell, = SurfaceWaterIntegratedModeling._downstream_cell(tstruct.spillfield, cell)

        _is_sink(cell, tstruct) && break # return if the new cell is a sink
    end

    # recompute inflow and footprint infiltration of affected supertraps
    reg = cell > 0 ? tstruct.regions[cell] : -1
    if reg > 0
        _update_Smin_Smax!(rateinfo, tstruct, tstruct.supertraps_of[reg])
    end

    amount = (sign(amount) != initial_sign) ? 0.0 : amount

    if reg > 0 && _is_trap_bottom(cell, tstruct)
        setinflow!(rateinfo, reg, getinflow(rateinfo, reg) + amount)
    end
    return amount
end

# ----------------------------------------------------------------------------
function _compute_initial_rateinfo(precipitation, infiltration, tstruct)
    if typeof(precipitation) <: Real
        precipitation = precipitation .* ones(size(tstruct.regions))
    end
    if typeof(infiltration) <: Real
        infiltration = infiltration .* ones(size(tstruct.regions))
    end

    num_traps = length(tstruct.spillpoints)
    
    # compute the basic runoff field with all traps empty
    runoff, reg_accum = watercourses(tstruct, [false],
                                     precipitation=precipitation,
                                     infiltration=infiltration)

    Smin = zeros(length(tstruct.footprints))
    Smax = zeros(length(tstruct.footprints))
    trap_inflow = vcat(reg_accum, zeros(numtraps(tstruct) - numregions(tstruct)))
    rateinfo = RateInfo(runoff, Smax, Smin, trap_inflow)
    _update_Smin_Smax!(rateinfo, tstruct, 1:length(tstruct.footprints))

    return rateinfo
end

# ----------------------------------------------------------------------------
function _update_Smin_Smax!(rateinfo, tstruct, traps)
        
    for i in traps
        # maximum remaining infiltration capacity rate within trap footprints
        setsmax!(rateinfo,
                 i,
                 -sum(min.(getrunoff(rateinfo, tstruct.footprints[i]), 0.0)))
    end

    for i in traps
        # minimum extra infiltration in submerged areas. For parent traps, this
        # amounts to the sum of the Smax of the immediate children traps
        setsmin!(rateinfo, i, sum(getsmax(rateinfo, subtrapsof(tstruct, i))))
    end
end

# ----------------------------------------------------------------------------
function _spills_to_parent(trap, tstruct, sgraph)
    target = get(sgraph.edges, trap, 0)
    parent = Graphs.outneighbors(tstruct.agglomerations, trap)

    return !isempty(parent) && parent[1] == target
end

# ----------------------------------------------------------------------------
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
            outflow = _is_parent(from, trap, tstruct) ? signed_outflow :
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
            outflow = _is_parent(to, trap, tstruct) ? signed_outflow :
                                                      max(signed_outflow, 0.0)
            _propagate_amount!(rateinfo, trap, outflow, tstruct, sgraph)
        end
    end
end

# ----------------------------------------------------------------------------
"""
will modify rateinfo 
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
