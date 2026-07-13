import Interpolations
using DifferentialEquations: solve, ODEProblem, VectorContinuousCallback, terminate!
export fill_sequence

"""
    fill_sequence(tstruct, weather_events, time_slack=0.0, infiltration=nothing, verbose=false)

Compute the sequence of events that describes how water on the terrain evolves over time.

For a given set of weather events, and a given terrain with associated trap structure,
determine the sequence of events that describes how the flow and accumulation of water
on the terrain changes over time.

Returns a `Vector{SpillEvent}` that expresses the discrete points in time
when different traps fills/empties, and the resulting changes on the surface flow
patterns.

# Arguments
- `tstruct::TrapStructure{Real}`: trap structure object describing the terrain traps
- `weather_events::Vector{WeatherEvent}`: Vector of weather events describing
                                          changes in weather over time
- `time_slack::Real`: tolerance for when to merge events that are close to each other
                      in time.  Should be set to zero or a small number.
                      @@@ NB: Support for this currently unimplemented.
- `infiltration::Union{Matrix{Real}, Nothing}`: 
                      grid of same shape as the terrain, giving the infiltration rate
                      at each gridcell.
- `verbose::Bool`: if `true`, dump progress information during computation

See also [`TrapStructure`](@ref), [`WeatherEvent`](@ref), [`SpillEvent`](@ref)
"""
function fill_sequence(tstruct::TrapStructure{<:Real},
            weather_events::Vector{WeatherEvent};
            time_slack::Real = 0.0,
            infiltration::Union{Matrix{<:Real}, Nothing} = nothing,
            dyn_traps::Vector{Int} = Int[],
            culverts::Vector{DynCulvert} = DynCulvert[],
            nbs::Vector{DynNBSPlacement} = DynNBSPlacement[],
            verbose::Bool=false)::Vector{SpillEvent}
    @assert !isempty(weather_events)

    num_traps = numtraps(tstruct)
    (num_traps == 0) && return SpillEvent[] # no traps -> empty sequence (return
                                            # type is Vector{SpillEvent})
    
    # initialize infiltration map from user input
    infiltration =
        (typeof(infiltration) == Nothing) ? zeros(size(tstruct.topography)) :
        (typeof(infiltration) <: Real)  ? ones(size(tstruct.topography)) * infiltration :
                                          infiltration
    # compute tables to support computation of trap water volume as function of
    # water level
    z_vol_tables = _compute_z_vol_tables(tstruct)

    # set initial filled_traps, cur_amounts and spillgraph
    filled_traps = Vector{Bool}(tstruct.trapvolumes .== 0.0)
    cur_amounts = fill(FilledAmount(0.0, weather_events[1].timestamp), num_traps)    
    sgraph = compute_complete_spillgraph(tstruct, filled_traps) 
    
    # start with empty sequence
    seq = Vector{SpillEvent}()

    # Persistent NBS layer storage, keyed by placement index — the single source of
    # truth for each NBS's water content across events and weather periods (see
    # network_context.jl's _nbs_layer_block / _store_nbs_state!).  Empty when no NBS.
    nbs_state = Dict{Int,Vector{Float64}}()

    # compute development within the duration of each weather event
    for (wix, we) in enumerate(weather_events)
        cur_time = we.timestamp
        end_time =
            (wix == length(weather_events)) ? Inf : weather_events[wix+1].timestamp
    
        @assert(all([ca.time == cur_time for ca ∈ cur_amounts]))

        # compute inflow/runoff/infiltration rates corresponding to the fill
        # graph and new rain rate
        rateinfo = compute_flow(sgraph, we.rain_rate, infiltration, tstruct, verbose; nbs=nbs)

        # build the dynamic networks for this weather period (§3) with the incremental
        # membership driver; empty when no dyn_traps/culverts/nbs are supplied, in which
        # case behaviour is unchanged.
        driver = build_network_driver(tstruct,
                                      _dyn_seeds(tstruct, dyn_traps, DynCulvert[]),
                                      culverts, findall(filled_traps), cur_amounts,
                                      rateinfo, infiltration, z_vol_tables,
                                      cur_time, end_time;
                                      nbs_placements = nbs, nbs_state = nbs_state)
        net_contexts    = driver.contexts
        net_trap_set    = _net_trap_set(net_contexts)
        net_covered_set = _net_covered_set(net_contexts, tstruct)

        # compute initial time estimates for when a trap become filled, or split
        # into subtraps (network traps are overridden from their network prediction)
        changetimeest = _set_initial_changetime_estimates(rateinfo, cur_amounts,
                                                          cur_time, filled_traps,
                                                          tstruct, net_contexts,
                                                          net_covered_set)

        # register the start of this weather event as a new, fully computed, spill event
        push!(seq, SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
                              copy(rateinfo.trap_inflow), copy(we.rain_rate),
                              copy(rateinfo.runoff)))

        # Will add new events to `seq`.  `sgraph`, `rateinfo`, `changetimeest`,
        # `filled_traps` and `cur_amounts` are also modified in the process
        _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
                                          filled_traps, cur_amounts, z_vol_tables,
                                          tstruct, infiltration, end_time, time_slack,
                                          net_contexts, net_trap_set, net_covered_set,
                                          verbose, driver, nbs, nbs_state)
    end

    return seq
end

# ----------------------------------------------------------------------------
function _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
                                           filled_traps, cur_amounts, z_vol_tables,
                                           tstruct, infiltration, endtime, time_slack,
                                           net_contexts, net_trap_set, net_covered_set,
                                           verbose, driver,
                                           nbs = DynNBSPlacement[],
                                           nbs_state = Dict{Int,Vector{Float64}}())
    cur_time = cur_amounts[1].time

    fill_updates = Vector{IncrementalUpdate{Bool}}()
    graph_updates = Vector{IncrementalUpdate{Int}}()

    count = 0
    while cur_time < endtime
        verbose && (mod(count+=1, 10) == 0) && println("Fill sequence iteration: ", count)

        cur_time, fill_updates =
            _identify_next_status_change!(changetimeest, cur_amounts, rateinfo,
                                          filled_traps, tstruct, z_vol_tables,
                                          cur_time, endtime,
                                          net_trap_set, net_covered_set)

        (cur_time > endtime || isempty(fill_updates)) && break # do not register
                                                               # more events
        # A network parent that emptied exposes its immediate children as transitory
        # traps; add them to the fill updates so they leave `filled_traps` too (§:empty).
        fill_updates = _expand_empty_fill_updates(fill_updates, net_contexts, tstruct)
        for u in fill_updates
            filled_traps[u.index] = u.value
        end
        # given changes in fill state, update spill graph.  Network traps (nodes AND
        # subsumed descendants) are kept out of the spillgraph — their routing is handled
        # by the ODE — so only non-covered fill updates drive it (§4/§8).
        non_net_fill_updates = filter(u -> u.index ∉ net_covered_set, fill_updates)
        graph_updates = update_spillgraph!(sgraph, non_net_fill_updates, tstruct)

        # given the updates ot the spill graph, update flow information in `rateinfo`
        setsavepoint!(rateinfo)
        _update_flow!(rateinfo, graph_updates, tstruct, sgraph)

        # Touch the affected networks (commit → rebuild → predict), refreshing net_trap_set /
        # net_covered_set and their changetimeest entries.  Runs AFTER _update_flow!.
        #
        # Touch gate (plan D4/§8): a network changes ONLY when a member fired (topology) or its
        # external inflow changed (dynamics) — and its growth boundary (the terminal unfilled
        # trap) is itself a node, so these two signals catch everything.  If no network meets
        # either, the whole pass is skipped: cached state/prediction stay exact (extern_inflow
        # unchanged).  Quiet events then cost no ODE solves — the dominant runtime win.
        old_covered   = net_covered_set
        net_committed = Dict{Int,Float64}()
        network_touched = false
        if !isempty(net_contexts)
            inflow_updated = Set(u.index for u in getinflowupdates(rateinfo))
            network_touched = any(net_contexts) do ctx
                any(u.index ∈ ctx.global_ix for u in fill_updates) ||
                    !isdisjoint(ctx.inflow_sources, inflow_updated)
            end
        end
        if network_touched
            net_contexts, net_trap_set, net_covered_set, net_committed =
                _touch_networks_driver!(driver, changetimeest, sgraph, tstruct,
                                        filled_traps, cur_amounts, rateinfo, z_vol_tables,
                                        infiltration, fill_updates, old_covered,
                                        cur_time, endtime)
            # traps that LEFT the networks need a fresh constant-rate changetime estimate
            for t in setdiff(old_covered, net_covered_set)
                changetimeest[t] = _compute_changetime_estimate(t, cur_amounts, cur_time,
                                                                rateinfo, filled_traps, tstruct)
            end
        end

        # update water amount in traps whose inflow rate is about to change, or that
        # just filled.  Network-covered traps are excluded — their amounts come from the
        # network state instead (§9).
        amount_updates = _update_affected_amounts(rateinfo, cur_amounts, filled_traps,
                                                  tstruct, z_vol_tables, cur_time,
                                                  net_covered_set)
        append!(amount_updates,
                [IncrementalUpdate(tix, FilledAmount(tstruct.trapvolumes[tix] -
                    tstruct.subvolumes[tix], cur_time))
                 for tix in [u.index for u in fill_updates]
                 if tix ∉ net_covered_set && tix ∉ old_covered])
        # network-trap amounts (nodes from ODE state, subsumed full at capacity, just
        # exited from their committed boundary value).  Only emitted when the networks
        # were touched: an untouched network's `state` is not committed to `cur_time`,
        # so its amounts are left to interpolate between touch events (plan §9).
        if network_touched && !(isempty(net_covered_set) && isempty(old_covered))
            append!(amount_updates,
                    _network_amount_updates(net_contexts, union(old_covered, net_covered_set),
                                            net_committed, tstruct, cur_amounts, cur_time))
        end

        # integrate the changes into the continously updated `cur_amounts` vector
        _apply_updates!(cur_amounts, amount_updates)

        # add current state to result
        push!(seq, SpillEvent(cur_time, amount_updates, fill_updates,
                              getinflowupdates(rateinfo), nothing,
                              getrunoffupdates(rateinfo)))
    end

    # make sure all amounts are exactly computed at end.  Network-covered traps are
    # excluded here — their end-of-period volumes follow the multi-trap ODE, not the
    # constant-rate projection, and are set by _finalize_networks! below (§10).
    for (trap, cur_fill) ∈ enumerate(cur_amounts)
        (trap ∈ net_covered_set) && continue
        if cur_fill.time < endtime
            cur_amounts[trap] =
                FilledAmount(_compute_exact_fill(rateinfo, cur_amounts, trap,
                                                 filled_traps, tstruct, endtime,
                                                 z_vol_tables, false),
                             min(cur_time, endtime))
        end
    end

    # advance every network to `endtime` and read its traps' boundary volumes from the
    # settled ODE state — the exact amounts the next weather period rebuilds from (§10)
    _finalize_networks!(cur_amounts, net_contexts, tstruct, infiltration,
                        z_vol_tables, cur_time, endtime, nbs_state)
end

# ----------------------------------------------------------------------------
function _update_affected_amounts(rateinfo, cur_amounts, filled_traps, tstruct,
                                  z_vol_tables, cur_time, net_covered_set = nothing)

    results = Vector{IncrementalUpdate{FilledAmount}}()

    for iu ∈ getinflowupdates(rateinfo)
        # network-covered traps get their amount from the network state, not the
        # constant-rate fill projection (§9)
        net_covered_set !== nothing && iu.index ∈ net_covered_set && continue
        amount = _compute_exact_fill(rateinfo, cur_amounts, iu.index,
                                     filled_traps, tstruct, cur_time, z_vol_tables, true)
        push!(results, IncrementalUpdate(iu.index, FilledAmount(amount, cur_time)))
    end
    return results
end

# ----------------------------------------------------------------------------
function _apply_updates!(amounts, updates)
    for up in updates
        amounts[up.index] = up.value
    end
end

# ----------------------------------------------------------------------------
function _compute_changetime_estimate(trap, cur_amounts, cur_time, rateinfo, filled_traps, tstruct)

    min_net_inflow = tr -> getinflow(rateinfo, tr) - getsmax(rateinfo, tr)
    max_net_inflow = tr -> getinflow(rateinfo, tr) - getsmin(rateinfo, tr)

    if filled_traps[trap]
        # Trap is currently full.  Return time when trap starts emptying
        if min_net_inflow(trap) >= 0
            # trap is currently full, but will stay that way as it has a
            # positive net inflow
            return ChangeTimeEstimate(false, Inf, Inf)
        end
        # The trap has a negative inflow, and will eventually empty.  If trap
        # does not have a parent, it will start emtpying right away.  Otherwise
        # it will start emptying as soon as it is not submerged by parent.
        parent = parentof(tstruct, trap)

        if isnothing(parent)
            return ChangeTimeEstimate(false, cur_time, cur_time) # start emptying now
        elseif filled_traps[parent]
            # parent must lose its 'filled' status before we can start
            # estimating how much time is required to empty it
            return ChangeTimeEstimate(false, Inf, Inf)
        else
            # The trap will become unfilled as soon as no longer submerged by
            # parent.  Compute when that will happen.
            parent_min_net_inflow = min_net_inflow(parent)
            parent_max_net_inflow = max_net_inflow(parent)
            if parent_max_net_inflow > 0
                # parent will not empty all the way to expose its subtraps
                return ChangeTimeEstimate(false, Inf, Inf)
            else
                parent_volume = cur_amounts[parent].amount
                min_time =
                    (parent_volume > 0.0) ? -parent_volume / parent_min_net_inflow : 0.0
                max_time =
                    (parent_volume > 0.0) ? -parent_volume / parent_max_net_inflow : 0.0

                return ChangeTimeEstimate(false,
                                          cur_time + min_time,
                                          cur_time + max_time)
            end
        end
    else
        # Trap is not yet full.  Return time when it switches to full.
        if min_net_inflow(trap) < 0.0 # inflow will become negative before trap
                                      # has been filled
            return ChangeTimeEstimate(true, Inf, Inf)
        else
            ownvolume = tstruct.trapvolumes[trap] - tstruct.subvolumes[trap]
            remaining = ownvolume - cur_amounts[trap].amount # remains to be filled
            starttime = cur_amounts[trap].time

            min_time = (remaining > 0.0) ? remaining / max_net_inflow(trap) : 0.0
            max_time = (remaining > 0.0) ? remaining / min_net_inflow(trap) : 0.0
            return ChangeTimeEstimate(true, min_time + starttime, max_time + starttime)
        end
    end
end


# ----------------------------------------------------------------------------
function _update_changetime_estimates!(changetimeest, cur_amounts, cur_time,
                                       rateinfo, filled_traps, tstruct,
                                       net_covered_set = nothing)

    inflow_updates = getinflowupdates(rateinfo)
    for update ∈ inflow_updates
        trap = update.index
        # covered traps (network nodes + subsumed descendants) are governed by the
        # network prediction, not the constant-rate estimate — skip them (§7)
        net_covered_set !== nothing && trap ∈ net_covered_set && continue
        changetimeest[trap] =
            _compute_changetime_estimate(trap, cur_amounts, cur_time,
                                         rateinfo, filled_traps, tstruct)
    end

end

# ----------------------------------------------------------------------------
function _set_initial_changetime_estimates(rateinfo, cur_amounts, cur_time,
                                           filled_traps, tstruct,
                                           net_contexts = nothing,
                                           net_covered_set = nothing)

    changetimeest =
        [_compute_changetime_estimate(trap, cur_amounts, cur_time, rateinfo,
                                      filled_traps, tstruct)
         for trap ∈ 1:numtraps(tstruct)]

    # override network traps from their network prediction (exact min == max), and
    # force every covered (node + subsumed) trap that is not a predicted trigger to Inf
    if net_contexts !== nothing && !isempty(net_contexts)
        _apply_network_changetimeest!(changetimeest, net_contexts, net_covered_set)
    end
    return changetimeest
end

# ----------------------------------------------------------------------------
# True iff every subtrap of `ix` is filled (allocation-free; short-circuits).
function all_subtraps_filled(tstruct, ix, filled_traps)
    @inbounds for c in subtrapsof(tstruct, ix)
        filled_traps[c] || return false
    end
    return true
end

function _identify_next_status_change!(changetimeest, cur_amounts, rateinfo,
                                       filled_traps, tstruct, z_vol_tables,
                                       cur_time, tmax,
                                       net_trap_set = nothing,
                                       net_covered_set = nothing)

    # update changetimeest for traps that have had their inflow changed
    _update_changetime_estimates!(changetimeest, cur_amounts, cur_time,
                                  rateinfo, filled_traps, tstruct, net_covered_set)
    # initialize return variables
    earliest_changetime = tmax

    # whether `ix` is a network node (keeps its exact prediction as a candidate)
    is_node(ix) = net_trap_set !== nothing && ix ∈ net_trap_set
    # whether `ix` is a subsumed descendant (covered but not a node): never a candidate
    is_subsumed(ix) = net_covered_set !== nothing && ix ∈ net_covered_set && !is_node(ix)

    # identify traps that may change their status before the earliest identified
    # changetime.  Explicit loop (no per-trap `filled_traps[subtrapsof(ix)]` temporary,
    # which allocated O(num_traps) small arrays every event and dominated GC on large
    # terrains); `all_subtraps_filled` short-circuits without allocating.
    num_traps = numtraps(tstruct)
    candidates = Int[]
    for ix in 1:num_traps
        (!is_subsumed(ix) && changetimeest[ix].min < earliest_changetime &&
         all_subtraps_filled(tstruct, ix, filled_traps)) && push!(candidates, ix)
    end

    candidate_mintimes = [changetimeest[ix].min for ix in candidates]

    cur_best_candidates = Vector{Int}()
    
    while !isempty(candidates)
        best_candidate_ix = argmin(candidate_mintimes)
        cur_candidate = candidates[best_candidate_ix]
        ctest = _compute_exact_changetime(cur_candidate, changetimeest, cur_amounts,
                                          rateinfo, tstruct, filled_traps, z_vol_tables)
        # update 'changetimeest'
        changetimeest[cur_candidate] = ctest
        # check if we found an improvement
        #(NB: ctest.min = ctest.max = exact changetime)
        if ctest.min < earliest_changetime
            cur_best_candidates = [cur_candidate] # discard previous, better found
            earliest_changetime = ctest.min
        elseif ctest.min == earliest_changetime
            push!(cur_best_candidates, cur_candidate)
        end
        # remove current candidate from list of candidates (since we have already
        # examined it), eliminate any other candidate with no possibility of
        # improving on currently found earliest changetime
        eliminate_ix = findall(candidate_mintimes .> earliest_changetime)
        unique!(sort!(push!(eliminate_ix, best_candidate_ix)))
        deleteat!(candidates, eliminate_ix)
        deleteat!(candidate_mintimes, eliminate_ix)
    end

    # determine fill updates. For candidates that filled, this refers to the
    # candidate itself.  For candidates that emptied below subtrap level, this
    # refers to any subtrap with negative inflow
    fill_updates = Vector{IncrementalUpdate{Bool}}()

    for cand ∈ cur_best_candidates
        push!(fill_updates, IncrementalUpdate{Bool}(cand, changetimeest[cand].filling))

        # this trap will not change again unless there is a weather change (in which
        # case all changetime estimates will be recomputed), so set it to infinity
        changetimeest[cand] = ChangeTimeEstimate(false, Inf, Inf)

        # Network nodes are re-predicted by the touch step (which also resets their
        # subsumed children), so the constant-rate subtrap recompute below is both
        # unnecessary and wrong for them — skip it.
        is_node(cand) && continue

        # Recompute changetimes for subtraps (which may change when parent
        # change its filled status)
        filled_traps[cand] = !filled_traps[cand] # temporary flip it while
                                                 # recomputing estimates
        children = subtrapsof(tstruct, cand)
        for c in children
            changetimeest[c] =
                _compute_changetime_estimate(c, cur_amounts, earliest_changetime,
                                             rateinfo, filled_traps, tstruct)
        end
        filled_traps[cand] = !filled_traps[cand] # flip it back so as not to have
                                                 # this function argument changed.
    end
    return earliest_changetime, fill_updates
end

# ----------------------------------------------------------------------------
function _compute_exact_fill(rateinfo, cur_amounts, trap, filled_traps, tstruct,
                             time, z_vol_tables, use_saved::Bool)
    if filled_traps[trap]
        return tstruct.trapvolumes[trap] - tstruct.subvolumes[trap]
    elseif !all(filled_traps[subtrapsof(tstruct, trap)])
        # this parent trap has unfilled childen, and should thus be empty
        return 0.0
    end
    
    vol, tstop = fill_trap_until(trap, rateinfo, cur_amounts[trap], time,
                                 tstruct, z_vol_tables, use_saved=use_saved)
    return vol
end
    
# ----------------------------------------------------------------------------
function _compute_exact_changetime(trap, changetimes, cur_amounts, rateinfo,
                                   tstruct, filled_traps, z_vol_tables)

    if changetimes[trap].min == changetimes[trap].max
        # we know the exact changetime already
        return changetimes[trap]
    end

    vol, starttime, tstop = 0.0, 0.0, nothing
    if filled_traps[trap]
        # We need to compute the exact time when this trap ends being completely filled.
        # If we got here, the exact changetime is not yet known, which should only happen
        # if the trap is submerged by a parent.  Determine how long it takes for that
        # parent to drain.
        @assert !changetimes[trap].filling
        par = parentof(tstruct, trap)
        starttime = cur_amounts[par].time
        @assert !isnothing(par) # if we got here, the trap should have a parent
        vol, tstop = fill_trap_until(par, rateinfo, cur_amounts[par],
                                     Inf, tstruct, z_vol_tables, use_saved=false)
    else
        # This trap is not yet filled.  The trap is however expected to be filling up.
        @assert changetimes[trap].filling
        starttime = cur_amounts[trap].time
        vol, tstop = fill_trap_until(trap, rateinfo, cur_amounts[trap],
                                     Inf, tstruct, z_vol_tables, use_saved=false)
    end
    tchange = isnothing(tstop) ? Inf : tstop + starttime
    return ChangeTimeEstimate(changetimes[trap].filling, tchange, tchange)
end

# ----------------------------------------------------------------------------
# for each trap, compute a table expressing increasing waterlevel (z) as
# function of trap water volume.  This is useful when we need to rapidly compute
# water z level as function of volume, as is done in _update_waterlevel().
function _compute_z_vol_tables(tstruct)

    N = numtraps(tstruct)
    
    zvt = [(Vector{Float64}(), Vector{Float64}()) for i in 1:N]

    for tix in 1:N
        trap_bottom = tstruct.topography[tstruct.footprints[tix]]
        children = subtrapsof(tstruct, tix)
        if !isempty(children)
            # This trap is a parent.  Its bottom is defined as being above its subtraps
            trap_bottom = max.(trap_bottom, tstruct.spillpoints[children[1]].elevation)
        end
        zsorted = sort(trap_bottom)
        append!(zsorted, tstruct.spillpoints[tix].elevation)
        vsorted = (zsorted .* (1:length(zsorted))) .- cumsum(zsorted) # trapvolume
        # get rid of duplicates by keeping only the last of each duplicate
        # (which also ensures that the corresponding volume in 'vsorted' will be
        # the correct one)
        keep = fill(true, length(zsorted))
        for i = 1:length(zsorted)-1
            keep[i] = zsorted[i+1] > zsorted[i]
        end
        zvt[tix] = (zsorted[keep], vsorted[keep])
    end
    return zvt
end

# ----------------------------------------------------------------------------
# Spillpoint-exclusion rule for a trap's infiltration: a footprint cell whose (raised)
# terrain bottom lies at or above the trap's spillpoint never ponds as part of the trap
# (water reaching that level spills out rather than pooling), so it carries no trap
# infiltration.  The canonical, allocation-free implementation is `_ponding_infiltration`
# (flow.jl), used where no local bottom is at hand (Smax updates, the
# `_compute_Smin_Smax_for_specific_trap!` utility).  The per-trap/per-event hot paths that
# already compute the raised bottom — `fill_trap_until` below and the network solver's
# `_build_trap_geometry` — apply the same `bottom >= spillpoint` test to that bottom
# directly.  (Excluding these cells is what keeps the wetted-area loss continuous at
# V = capacity and removes the fill/unfill chatter of a trap pinned at its spillpoint.)

# ----------------------------------------------------------------------------
# first return value: new amount of water in trap: any value between 0.0 (empty) and
#                     the trap volume (or trap self-volume, if subtraps are excluded)
# second return value: timepoint where the trap filled (if it became full) or emptied
#                      (if it emptied).  If it neither reached filled or emptied status
#                      during  before `endtim`, then the second return value will be
#                      `nothing`.
function fill_trap_until(trap, rateinfo, cur_amount, endtime, tstruct, z_vol_tables;
                         use_saved=false)

    footprint = tstruct.footprints[trap]
    trap_bottom = tstruct.topography[footprint]
    tvolume = tstruct.trapvolumes[trap] - tstruct.subvolumes[trap]
    
    Smin   = use_saved ? getsavedsmin(rateinfo, trap)   : getsmin(rateinfo, trap)
    Smax   = use_saved ? getsavedsmax(rateinfo, trap)   : getsmax(rateinfo, trap)
    inflow = use_saved ? getsavedinflow(rateinfo, trap) : getinflow(rateinfo, trap)

    if (trap > numregions(tstruct))
        children = subtrapsof(tstruct, trap)
        trap_bottom = max.(trap_bottom, tstruct.spillpoints[children[1]].elevation)
    else
        tvolume = tstruct.trapvolumes[trap]
    end

    if Smax == Smin
        # net rate will not depend on the degree of fill, and we do not need to
        # solve an ODE to get the time to trap filled (or emptied)
        accum_rate = inflow - Smax
        (accum_rate == 0.0) && return (cur_amount.amount, nothing) # no change in fill
                                                                   # amount
        dt = (accum_rate > 0) ?
                (tvolume - cur_amount.amount) / accum_rate : # time to full
                cur_amount.amount / abs(accum_rate)          # time to empty
        t = cur_amount.time + dt
        reached = (t <= endtime)

        return (reached) ?
            (accum_rate > 0.0 ? tvolume : 0.0, t) :
            (cur_amount.amount + (endtime - cur_amount.time) * accum_rate, nothing)
    end
    # if we got here, the amount of infiltration at any time depends on how much
    # the trap has been filled (since the footprint of its water content will vary).
    # We must solve an ODE to determine how much it fills or empties over the time period.
    fprint_infil =
        use_saved ? [-min(getsavedrunoff(rateinfo, i), 0.0) for i in footprint] :
                    -min.(getrunoff(rateinfo, footprint), 0.0)
    # Cells whose terrain lies at or above the spillpoint never pond as part of the trap:
    # drop their infiltration so the wetted-area loss is continuous at V = capacity
    # (removing the fill/unfill chatter the discontinuity used to cause), consistent with
    # `_ponding_infiltration` and the network solver.  Test the actual terrain height, not
    # the raised `trap_bottom` (which would over-null a degenerate zero-own-volume parent).
    _sp = Float64(tstruct.spillpoints[trap].elevation)   # concrete: Spillpoint.elevation is ::Real
    @inbounds for k in eachindex(fprint_infil)
        tstruct.topography[footprint[k]] >= _sp && (fprint_infil[k] = 0.0)
    end

    infilfun = ixs -> sum(fprint_infil[ixs])
    v0 = [cur_amount.amount]
    dvdt = _setup_dvdt(trap_bottom, tvolume, infilfun, inflow,
                       tstruct.spillpoints[trap], z_vol_tables[trap])
    dv0 = [0.0]
    dvdt(dv0, v0, 0, 0) # compute derivative at start of integration (used to
                        # detect sign changes in derivative later
    # Empty traps w/negative rate or full traps w/positive rate are already finished
    (v0[1] == 0.0 && dv0[1] <= 0) && return (0.0, cur_amount.time)
    (v0[1] == tvolume && dv0[1] >= 0) && return (tvolume, cur_amount.time)
    
    function condition(out, v, t, integrator)
        out[1] = tvolume - v[1] # trap full
        out[2] = v[1]           # trap empty
        deriv = [0.0]
        dvdt(deriv, v, 0, t)
        out[3] = dv0[1] * deriv[1] # stagnation: sign of time derivative changed, meaning 
    end                            #             it reached zero along the way

    condition_reached = [0]
    
    function affect!(integrator, ix)
        # DiffEq v8+ passes ix as Vector{Int8} flags; older versions passed a scalar index
        condition_reached[] = isa(ix, AbstractVector) ? findfirst(!iszero, ix) : ix
        terminate!(integrator)
    end
    cb = VectorContinuousCallback(condition, affect!, 3)
    # Cap dt: Inf tspan causes DiffEq's adaptive stepper to reach t ~ 1e307, at which
    # point the internal event-check range(t, Inf, n) overflows.  1e12 is far beyond
    # any physical horizon; if no callback fires within it, the caller sees tstop=nothing.
    dt = min(endtime - cur_amount.time, 1e12)
    sol = solve(ODEProblem(dvdt, v0, [0, dt]), callback = cb, abstol=1e-6, reltol=1e-4)

    return (sol.u[end][1], # amount of water at end of integration
            (condition_reached[1] ∈ [1,2]) ? sol.t[end] : nothing) # integration duration
end

# ----------------------------------------------------------------------------
function _setup_dvdt(trap_bottom, trapvolume, infilfun, inflow, spillpoint, zvtable)
    # note: for parent traps, trap_bottom will represent the bottom topography
    # except within the footprint of its children's traps, where trap_bottom
    # will have the value of the upper child's spillpoint

    zmin = minimum(trap_bottom)
    zspan = spillpoint.elevation - zmin

    v2z = length(zvtable[2]) == 1 ?
        x -> zmin : # degenerate case
        Interpolations.linear_interpolation(zvtable[2], zvtable[1],
                                             extrapolation_bc=Interpolations.Line());
        # Interpolations.LinearInterpolation(zvtable[2], zvtable[1],
        #                                      extrapolation_bc=Interpolations.Line());
    function _dvdt(dv, v, p, t)
        z = (v[1] <= 0)          ? zmin :
            (v[1] >= trapvolume) ? Inf :
                                  v2z(v[1])
        return dv[1] = inflow - infilfun(trap_bottom .<= z)
    end

    return _dvdt # return a closure
end
