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
                      @@ NB: Support for this currently unimplemented.
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
            verbose::Bool=false)::Vector{SpillEvent} 
    @assert !isempty(weather_events)

    num_traps = numtraps(tstruct)
    (num_traps == 0) && return # if the terrain has no traps, there is nothing to do
    
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

    # compute development within the duration of each weather event
    for (wix, we) in enumerate(weather_events)
        cur_time = we.timestamp
        end_time =
            (wix == length(weather_events)) ? Inf : weather_events[wix+1].timestamp
    
        @assert(all([ca.time == cur_time for ca ∈ cur_amounts]))

        # compute inflow/runoff/infiltration rates corresponding to the fill
        # graph and new rain rate
        rateinfo = compute_flow(sgraph, we.rain_rate, infiltration, tstruct, verbose)

        # compute initial time estimates for when a trap become filled, or split
        # into subtraps
        changetimeest = _set_initial_changetime_estimates(rateinfo, cur_amounts,
                                                          cur_time, filled_traps,
                                                          tstruct)

        # register the start of this weather event as a new, fully computed, spill event
        push!(seq, SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
                              copy(rateinfo.trap_inflow), copy(we.rain_rate),
                              copy(rateinfo.runoff)))

        # Will add new events to `seq`.  `sgraph`, `rateinfo`, `changetimeest`,
        # `filled_traps` and `cur_amounts` are also modified in the process
        _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
                                          filled_traps, cur_amounts, z_vol_tables,
                                          tstruct, infiltration, end_time, time_slack,
                                          verbose)
    end

    return seq
end

# ----------------------------------------------------------------------------
function _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
                                           filled_traps, cur_amounts, z_vol_tables,
                                           tstruct, infiltration, endtime, time_slack,
                                           verbose)
    cur_time = cur_amounts[1].time

    fill_updates = Vector{IncrementalUpdate{Bool}}()
    graph_updates = Vector{IncrementalUpdate{Int}}()

    count = 0
    while cur_time < endtime
        verbose && (mod(count+=1, 10) == 0) && println("Fill sequence iteration: ", count)
        
        cur_time, fill_updates =
            _identify_next_status_change!(changetimeest, cur_amounts, rateinfo,
                                          filled_traps, tstruct, z_vol_tables,
                                          cur_time, endtime)

        (cur_time > endtime || isempty(fill_updates)) && break # do not register
                                                               # more events
        for u in fill_updates
            filled_traps[u.index] = u.value
        end
        # given changes in fill state, update spill graph
        graph_updates = update_spillgraph!(sgraph, fill_updates, tstruct)

        # given the updates ot the spill graph, update flow information in `rateinfo`
        setsavepoint!(rateinfo)
        _update_flow!(rateinfo, graph_updates, tstruct, sgraph)
        
        # update water amount in traps whose inflow rate is about to change, or
        # that just filled
        amount_updates = _update_affected_amounts(rateinfo, cur_amounts, filled_traps,
                                                  tstruct, z_vol_tables, cur_time)
        append!(amount_updates,
                [IncrementalUpdate(tix, FilledAmount(tstruct.trapvolumes[tix] -
                    tstruct.subvolumes[tix], cur_time))
                 for tix in [u.index for u in fill_updates]])

        # integrate the changes into the continously updated `cur_amounts` vector
        _apply_updates!(cur_amounts, amount_updates)

        # add current state to result
        push!(seq, SpillEvent(cur_time, amount_updates, fill_updates,
                              getinflowupdates(rateinfo), nothing,
                              getrunoffupdates(rateinfo)))
    end

    # make sure all amounts are exactly computed at end
    for (trap, cur_fill) ∈ enumerate(cur_amounts)
        if cur_fill.time < endtime
            cur_amounts[trap] =
                FilledAmount(_compute_exact_fill(rateinfo, cur_amounts, trap,
                                                 filled_traps, tstruct, endtime,
                                                 z_vol_tables, false),
                             min(cur_time, endtime))
        end
    end
end

# ----------------------------------------------------------------------------
function _update_affected_amounts(rateinfo, cur_amounts, filled_traps, tstruct,
                                  z_vol_tables, cur_time)

    results = Vector{IncrementalUpdate{FilledAmount}}()

    for iu ∈ getinflowupdates(rateinfo)
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

    new_inflow = getinflow(rateinfo, trap)
    # NB: when computing min_net_inflow and max_net_inflow, we must keep in mind
    # that the flow into a parent trap has already had the infiltration of the
    # subtrap deducted.  The maximum inflow is therefore equal the inflow, and
    # the minimum inflow will have the additional amount deducted that is the
    # max infiltration of the parent trap (Smax) less the quantity that has
    # already been deducted from the inflow (Smin).
    min_net_inflow = tr -> getinflow(rateinfo, tr) -
                           (getsmax(rateinfo, tr) - getsmin(rateinfo, tr))
    max_net_inflow = tr -> getinflow(rateinfo, tr)

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
                starttime = cur_amounts[parent].time
                min_time =
                    (parent_volume > 0.0) ? -parent_volume / parent_min_net_inflow : 0.0
                max_time =
                    (parent_volume > 0.0) ? -parent_volume / parent_max_net_inflow : 0.0

                return ChangeTimeEstimate(false,
                                          min_time + starttime,
                                          max_time + starttime)
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
                                       rateinfo, filled_traps, tstruct)

    inflow_updates = getinflowupdates(rateinfo)
    for update ∈ inflow_updates
        trap = update.index
        changetimeest[trap] =
            _compute_changetime_estimate(trap, cur_amounts, cur_time,
                                         rateinfo, filled_traps, tstruct)
    end

end

# ----------------------------------------------------------------------------
function _set_initial_changetime_estimates(rateinfo, cur_amounts, cur_time,
                                           filled_traps, tstruct)

    return [_compute_changetime_estimate(trap, cur_amounts, cur_time, rateinfo,
                                         filled_traps, tstruct)
            for trap ∈ 1:numtraps(tstruct)]
end

# ----------------------------------------------------------------------------
function _identify_next_status_change!(changetimeest, cur_amounts, rateinfo, 
                                       filled_traps, tstruct, z_vol_tables,
                                       cur_time, tmax)

    # update changetimeest for traps that have had their inflow changed
    _update_changetime_estimates!(changetimeest, cur_amounts, cur_time,
                                  rateinfo, filled_traps, tstruct)
    # initialize return variables
    earliest_changetime = tmax

    # identify traps that may change their status before the earliest identified
    # changetime
    num_traps = numtraps(tstruct)
    candidates = findall([all(filled_traps[subtrapsof(tstruct, ix)]) && 
                          changetimeest[ix].min < earliest_changetime
                          for ix ∈ 1:num_traps])
    
    candidate_mintimes = [x.min for x in changetimeest[candidates]]

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

        # Recompute changetimes for subtraps (which may change when parent
        # change its filled status)
        filled_traps[cand] = !filled_traps[cand] # temporary flip it while
                                                 # recomputing estimates
        children = subtrapsof(tstruct, cand)
        for c in children
            changetimeest[c] =
                _compute_changetime_estimate(c, cur_amounts, cur_time,
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
        accum_rate = inflow # NB: Smin (=Smax) has already been deducted from inflow
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

    infilfun = ixs -> sum(fprint_infil[ixs]) - Smin # discount Smin since it has been
                                                    # deducted when computing inflow
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
        condition_reached[] = ix
        terminate!(integrator)
    end
    cb = VectorContinuousCallback(condition, affect!, 3)
    dt = endtime - cur_amount.time
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
