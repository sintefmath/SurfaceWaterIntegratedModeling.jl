import Interpolations
using DifferentialEquations: solve, ODEProblem, VectorContinuousCallback, terminate!
export fill_sequence

"""
    fill_sequence(tstruct, weather_events; time_slack=0.0, infiltration=nothing,
                  dyn_traps=Int[], culverts=DynCulvert[], nbs=DynNBSPlacement[], verbose=false)

Compute the sequence of events that describes how water on the terrain evolves over time.

For a given set of weather events, and a given terrain with associated trap structure,
determine the sequence of events that describes how the flow and accumulation of water
on the terrain changes over time.

When any of `dyn_traps`, `culverts`, or `nbs` is supplied, the affected traps are solved
as coupled multi-trap ODE networks (see `dynamics/`) rather than the constant-rate
analytic path; the rest of the terrain is unaffected.

Returns a `Vector{SpillEvent}` that expresses the discrete points in time
when different traps fills/empties, and the resulting changes on the surface flow
patterns.

# Arguments
- `tstruct::TrapStructure{Real}`: trap structure object describing the terrain traps
- `weather_events::Vector{WeatherEvent}`: Vector of weather events describing
                                          changes in weather over time
- `time_slack::Real`: tolerance for when to merge events that are close to each other
                      in time.  Should be set to zero or a small number.
                      @@@ NB: Support for this currently verified.
- `infiltration::Union{Matrix{Real}, Nothing}`:
                      grid of same shape as the terrain, giving the infiltration rate
                      at each gridcell.  NBS footprints are forced to zero infiltration
                      (the layer model accounts for it instead).
- `dyn_traps::Vector{Int}`: trap indices to solve as dynamic networks (e.g. traps whose
                      inflow is coupled through culverts or NBS).
- `culverts::Vector{DynCulvert}`: culverts linking cells across the terrain; both
                      endpoints seed a dynamic network.
- `nbs::Vector{DynNBSPlacement}`: nature-based solution placements.  Each applies a signed
                      correction to the surface flow inside the network solver; the rest of
                      the code is NBS-oblivious.
- `verbose::Bool`: if `true`, dump progress information during computation

See also [`TrapStructure`](@ref), [`WeatherEvent`](@ref), [`SpillEvent`](@ref),
[`DynCulvert`](@ref), [`DynNBSPlacement`](@ref)
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
    (numtraps(tstruct) == 0) && return SpillEvent[] # no traps -> empty sequence
    
    # initialize infiltration map from user input
    infiltration = (typeof(infiltration) == Nothing) ? zeros(size(tstruct.topography)) :
                   (typeof(infiltration) <: Real)    ? ones(size(tstruct.topography)) * infiltration :
                                                      copy(infiltration)
    # ensuring no infiltration overlaps with NBS footprints
    infiltration[vcat([n.footprint for n in nbs]...)] .= 0.0

    # tables mapping water level -> trap volume
    z_vol_tables = _compute_z_vol_tables(tstruct)

    # set initial filled_traps, cur_amounts and spillgraph
    filled_traps = Vector{Bool}(tstruct.trapvolumes .== 0.0)
    cur_amounts = fill(FilledAmount(0.0, weather_events[1].timestamp), numtraps(tstruct))    
    sgraph = compute_complete_spillgraph(tstruct, filled_traps) # traps are nodes here
    
    # start with empty sequence.  `nbs_state` represents the persistent NBS layer storage
    seq = Vector{SpillEvent}()
    nbs_state = Dict{Int,Vector{Float64}}()

    # compute development within the duration of each weather event
    for (wix, we) in enumerate(weather_events)
        cur_time = we.timestamp
        end_time =
            (wix == length(weather_events)) ? Inf : weather_events[wix+1].timestamp
        @assert(all([ca.time == cur_time for ca ∈ cur_amounts]))

        # flow rates for the current spillgraph and rain rate
        rateinfo = compute_flow(sgraph, we.rain_rate, infiltration, tstruct, verbose)

        # supporting structures for the dynamic networks (empty if there are no seeds)
        driver = build_network_driver(tstruct,
                                      _dyn_seeds(tstruct, dyn_traps, DynCulvert[]),
                                      culverts, findall(filled_traps), cur_amounts,
                                      rateinfo, infiltration, z_vol_tables,
                                      cur_time, end_time;
                                      nbs_placements = nbs, precipitation = we.rain_rate,
                                      nbs_state = nbs_state)
        # initial changetime estimates (network traps overridden from their prediction)
        changetimeest = _set_initial_changetime_estimates(rateinfo, cur_amounts,
                                                          cur_time, filled_traps, tstruct,
                                                          driver.contexts,
                                                          _net_covered_set(driver.contexts, tstruct))

        # the weather event's start is itself a full spill event
        push!(seq, SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
                              copy(rateinfo.trap_inflow), copy(we.rain_rate), copy(rateinfo.runoff)))

        # extends seq; also mutates sgraph, rateinfo, changetimeest, filled_traps, cur_amounts
        _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest, filled_traps,
                                          cur_amounts, z_vol_tables, tstruct, infiltration,
                                          end_time, time_slack, verbose, driver)
    end

    return seq
end

# ----------------------------------------------------------------------------
"""
    _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest, filled_traps,
                                      cur_amounts, z_vol_tables, tstruct, infiltration,
                                      endtime, time_slack, verbose, driver) -> nothing

Run the event loop for one weather period, from the time held in `cur_amounts` up to `endtime`.

Each iteration resolves the next status change, flips the affected traps, republishes the
spillgraph and flow, touches any affected dynamic network, and appends one `SpillEvent`.  At
the period boundary the non-covered traps are settled analytically and the networks are
advanced to `endtime`, so their volumes seed the next period.

Mutates nearly everything it is handed: `seq` (extended by one event per iteration),
`filled_traps`, `cur_amounts`, `changetimeest`, `sgraph`, `rateinfo`, and `driver`.
`tstruct`, `infiltration` and `z_vol_tables` are read only.  `time_slack` is accepted but not
yet used (see [`fill_sequence`](@ref)).
"""
function _fill_sequence_for_weather_event!(seq, sgraph, rateinfo, changetimeest,
                                           filled_traps, cur_amounts, z_vol_tables,
                                           tstruct, infiltration, endtime, time_slack,
                                           verbose, driver)
    cur_time = cur_amounts[1].time

    # network membership caches, kept in sync with driver.contexts (the single source of
    # truth); _touch_affected_networks! updates them in place as the topology changes
    net_trap_set    = _net_trap_set(driver.contexts)
    net_covered_set = _net_covered_set(driver.contexts, tstruct)

    count = 0
    while cur_time < endtime
        verbose && (mod(count+=1, 10) == 0) && println("Fill sequence iteration: ", count)

        cur_time, fill_updates =
            _identify_next_status_change!(changetimeest, cur_amounts, rateinfo,
                                          filled_traps, tstruct, z_vol_tables,
                                          cur_time, endtime,
                                          net_trap_set, net_covered_set)

        (cur_time > endtime || isempty(fill_updates)) && break # no more events to register

        # an emptied dynamic network parent exposes its children as transitory traps; flush
        fill_updates = _expand_empty_network_fill_updates(fill_updates, driver.contexts, tstruct)
        for u in fill_updates
            filled_traps[u.index] = u.value
        end

        # network-covered traps route via the ODE, so keep them out of the spillgraph
        non_net_fill_updates = filter(u -> u.index ∉ net_covered_set, fill_updates)
        graph_updates = update_spillgraph!(sgraph, non_net_fill_updates, tstruct)
        setsavepoint!(rateinfo)
        _update_flow!(rateinfo, graph_updates, tstruct, sgraph) # agnostic of dynamic networks

        # snapshot the covered set before the touch mutates it (needed below for the amounts)
        old_covered = copy(net_covered_set)
        network_touched =
            _touch_affected_networks!(net_trap_set, net_covered_set, driver, changetimeest,
                                      sgraph, rateinfo, filled_traps, cur_amounts, z_vol_tables,
                                      infiltration, fill_updates, tstruct, old_covered,
                                      cur_time, endtime)

        amount_updates = _collect_amount_updates(rateinfo, cur_amounts, filled_traps, tstruct,
                                                 z_vol_tables, cur_time, fill_updates, driver,
                                                 net_covered_set, old_covered, network_touched)
        _apply_updates!(cur_amounts, amount_updates)

        push!(seq, SpillEvent(cur_time, amount_updates, fill_updates,
                              getinflowupdates(rateinfo), nothing, getrunoffupdates(rateinfo)))
    end

    _settle_noncovered_at_endtime!(cur_amounts, rateinfo, filled_traps, tstruct,
                                   z_vol_tables, net_covered_set, cur_time, endtime)
    # advance every network to endtime; its settled ODE volumes seed the next weather period
    _finalize_networks!(cur_amounts, driver.contexts, tstruct, infiltration,
                        z_vol_tables, cur_time, endtime, driver.nbs_state)
end

# ----------------------------------------------------------------------------
"""
    _touch_affected_networks!(net_trap_set, net_covered_set, driver, changetimeest, sgraph,
                              rateinfo, filled_traps, cur_amounts, z_vol_tables, infiltration,
                              fill_updates, tstruct, old_covered, cur_time, endtime) -> Bool

Touch the dynamic networks affected by this event.  A network is touched when one of its member
traps fired or one of its inflow regions changed; an untouched network keeps its exact cached
state, so quiet events cost no ODE solves.  On a touch, commit / rebuild / re-predict through
the driver and give every trap that just left a network a fresh constant-rate changetime
estimate.

Mutates `net_trap_set` and `net_covered_set` in place, refreshing them to the new membership,
plus `driver`, `changetimeest`, `sgraph` and `rateinfo` via `_touch_networks_driver!`.
`filled_traps`, `cur_amounts`, `fill_updates` and `old_covered` (the covered set before this
event) are read only.

Returns `true` if a touch happened; `false` leaves every argument untouched.
"""
function _touch_affected_networks!(net_trap_set, net_covered_set, driver, changetimeest,
                                   sgraph, rateinfo, filled_traps, cur_amounts, z_vol_tables,
                                   infiltration, fill_updates, tstruct, old_covered,
                                   cur_time, endtime)
    isempty(driver.contexts) && return false
    inflow_updated = Set(u.index for u in getinflowupdates(rateinfo))
    # a network is touched if any of its member traps fired or any of its inflow regions changed
    touched = any(driver.contexts) do ctx
        any(u.index ∈ ctx.global_ix for u in fill_updates) ||
            !isdisjoint(ctx.inflow_regions, inflow_updated)
    end
    touched || return false

    nts, ncs = _touch_networks_driver!(driver, changetimeest, sgraph, tstruct,
                                       filled_traps, cur_amounts, rateinfo, z_vol_tables,
                                       infiltration, fill_updates, old_covered,
                                       cur_time, endtime)

    # for those traps that just left the dynamic networks, recompute their changetime estimates
    # in the regular constant-rate way
    for t in setdiff(old_covered, ncs)
        changetimeest[t] = _compute_changetime_estimate(t, cur_amounts, cur_time,
                                                        rateinfo, filled_traps, tstruct)
    end
    empty!(net_trap_set);    union!(net_trap_set, nts)
    empty!(net_covered_set); union!(net_covered_set, ncs)
    return true
end

# ----------------------------------------------------------------------------
"""
    _collect_amount_updates(rateinfo, cur_amounts, filled_traps, tstruct, z_vol_tables,
                            cur_time, fill_updates, driver, net_covered_set, old_covered,
                            network_touched) -> Vector{IncrementalUpdate{FilledAmount}}

The amount updates for one event, in three parts: constant-rate traps whose inflow changed,
traps that just flipped fill status (which sit at exactly their own capacity at that instant),
and — when `network_touched` — the network-owned amounts read from the driver.

Network-covered traps are excluded from the first two parts: their volumes come from the ODE
state, not from a constant-rate projection.  `old_covered` is the covered set before the touch.
Nothing is mutated; the updates are returned for the caller to apply.
"""
function _collect_amount_updates(rateinfo, cur_amounts, filled_traps, tstruct, z_vol_tables,
                                 cur_time, fill_updates, driver, net_covered_set, old_covered,
                                 network_touched)

    # amounts of non-network traps whose inflow changed
    amount_updates = _update_affected_amounts(rateinfo, cur_amounts, filled_traps, tstruct,
                                              z_vol_tables, cur_time, net_covered_set)
    # a trap sits at its own capacity at the instant its fill-status flips; record it
    # (network traps excluded — their amounts come from the network state below)
    append!(amount_updates,
            [IncrementalUpdate(tix, FilledAmount(tstruct.trapvolumes[tix] -
                tstruct.subvolumes[tix], cur_time))
             for tix in [u.index for u in fill_updates]
             if tix ∉ net_covered_set && tix ∉ old_covered])

    # update amounts of network traps if a network was touched (a touched network always has
    # at least one covered trap)
    if network_touched
        append!(amount_updates,
                _network_amount_updates(driver.contexts, union(old_covered, net_covered_set),
                                        driver.vol_by_trapix, tstruct, cur_amounts, cur_time))
    end
    return amount_updates
end

# ----------------------------------------------------------------------------
"""
    _settle_noncovered_at_endtime!(cur_amounts, rateinfo, filled_traps, tstruct, z_vol_tables,
                                   net_covered_set, cur_time, endtime) -> cur_amounts

Project every non-covered trap's amount to exactly `endtime`, closing out the weather period.
Traps in `net_covered_set` are skipped — they follow the multi-trap ODE and are settled by
`_finalize_networks!` instead.

Mutates and returns `cur_amounts`; every other argument is read only.
"""
function _settle_noncovered_at_endtime!(cur_amounts, rateinfo, filled_traps, tstruct,
                                        z_vol_tables, net_covered_set, cur_time, endtime)
    for (trap, cur_fill) in enumerate(cur_amounts)
        (trap ∈ net_covered_set) && continue
        cur_fill.time < endtime || continue
        cur_amounts[trap] = FilledAmount(
            _compute_exact_fill(rateinfo, cur_amounts, trap, filled_traps, tstruct,
                                endtime, z_vol_tables, false),
            min(cur_time, endtime))
    end
    return cur_amounts
end

# ----------------------------------------------------------------------------
"""
    _update_affected_amounts(rateinfo, cur_amounts, filled_traps, tstruct, z_vol_tables,
                             cur_time, net_covered_set = nothing)
        -> Vector{IncrementalUpdate{FilledAmount}}

One amount update per trap whose inflow changed, each projected to `cur_time` under the trap's
*saved* rates — the rates it was evolving under since the last event, not the newly-updated
ones.

Traps in `net_covered_set` are skipped, taking their amount from the network state instead;
`nothing` (the default) means no networks are active.  Nothing is mutated.
"""
function _update_affected_amounts(rateinfo, cur_amounts, filled_traps, tstruct,
                                  z_vol_tables, cur_time, net_covered_set = nothing)

    results = Vector{IncrementalUpdate{FilledAmount}}()

    for iu ∈ getinflowupdates(rateinfo)
        # covered traps get their amount from the network state, not this projection
        net_covered_set !== nothing && iu.index ∈ net_covered_set && continue
        amount = _compute_exact_fill(rateinfo, cur_amounts, iu.index,
                                     filled_traps, tstruct, cur_time, z_vol_tables, true)
        push!(results, IncrementalUpdate(iu.index, FilledAmount(amount, cur_time)))
    end
    return results
end

# ----------------------------------------------------------------------------
"""
    _apply_updates!(amounts, updates) -> nothing

Write each `IncrementalUpdate` in `updates` into `amounts` at its own index.  Mutates
`amounts`.
"""
function _apply_updates!(amounts, updates)
    for up in updates
        amounts[up.index] = up.value
    end
end

# ----------------------------------------------------------------------------
"""
    _compute_changetime_estimate(trap, cur_amounts, cur_time, rateinfo, filled_traps, tstruct)
        -> ChangeTimeEstimate

Bracket when `trap` next flips fill status, as a `[min, max]` window of absolute times.

The bounds come from the trap's net inflow range: infiltration varies with the wetted
footprint, so the net rate lies between `inflow - Smax` and `inflow - Smin`.  When those
coincide the rate is fill-independent and the estimate is *exact* (`min == max`), which lets
`_compute_exact_changetime` short-circuit.  A full trap is bracketed on when it starts
emptying, an unfilled one on when it fills; a submerged trap waits on its parent draining.

Nothing is mutated.  This function is unaffected by the presence of dynamic networks.
"""
function _compute_changetime_estimate(trap, cur_amounts, cur_time, rateinfo, filled_traps, tstruct)

    min_net_inflow = tr -> getinflow(rateinfo, tr) - getsmax(rateinfo, tr)
    max_net_inflow = tr -> getinflow(rateinfo, tr) - getsmin(rateinfo, tr)

    if filled_traps[trap]
        # full trap: return when it starts emptying
        if min_net_inflow(trap) >= 0
            return ChangeTimeEstimate(false, Inf, Inf) # positive net inflow: stays full
        end
        # negative inflow, so it will empty — right away if it has no parent, otherwise
        # once the parent no longer submerges it
        parent = parentof(tstruct, trap)

        if isnothing(parent)
            return ChangeTimeEstimate(false, cur_time, cur_time) # start emptying now
        elseif filled_traps[parent]
            return ChangeTimeEstimate(false, Inf, Inf) # wait until the parent loses 'filled'
        else
            # unfills once the parent has drained enough to expose it; find when
            parent_min_net_inflow = min_net_inflow(parent)
            parent_max_net_inflow = max_net_inflow(parent)
            if parent_max_net_inflow > 0
                return ChangeTimeEstimate(false, Inf, Inf) # parent won't empty to expose subtraps
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
        # not yet full: return when it switches to full
        if min_net_inflow(trap) < 0.0 # inflow turns negative before it fills
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
"""
    _update_changetime_estimates!(changetimeest, cur_amounts, cur_time, rateinfo, filled_traps,
                                  tstruct, net_covered_set = nothing) -> nothing

Refresh the changetime estimate of every trap whose inflow changed in `rateinfo`.

Traps in `net_covered_set` are skipped: they are governed by the network prediction, not by
this constant-rate estimate.  `nothing` (the default) means no networks are active.  Mutates
`changetimeest` at the affected indices.
"""
function _update_changetime_estimates!(changetimeest, cur_amounts, cur_time,
                                       rateinfo, filled_traps, tstruct,
                                       net_covered_set = nothing)

    inflow_updates = getinflowupdates(rateinfo)
    for update ∈ inflow_updates
        trap = update.index
        # covered traps are governed by the network prediction, not this estimate — skip them
        net_covered_set !== nothing && trap ∈ net_covered_set && continue
        changetimeest[trap] =
            _compute_changetime_estimate(trap, cur_amounts, cur_time,
                                         rateinfo, filled_traps, tstruct)
    end

end

# ----------------------------------------------------------------------------
"""
    _set_initial_changetime_estimates(rateinfo, cur_amounts, cur_time, filled_traps, tstruct,
                                      net_contexts = nothing, net_covered_set = nothing)
        -> Vector{ChangeTimeEstimate}

One changetime estimate per trap of `tstruct`, for the start of a weather period.

When networks are active, the covered traps are then overridden from their network prediction —
an exact estimate for each predicted trigger, `Inf` for the rest — so the branch-and-bound in
`_identify_next_status_change!` resolves them from the ODE rather than the constant-rate logic.
Nothing is mutated; the vector is freshly built.
"""
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
"""
    all_subtraps_filled(tstruct, ix, filled_traps) -> Bool

True iff every subtrap of trap `ix` is filled.  Allocation-free and short-circuiting: a trap
cannot start filling before its children have, so this gates the candidate scan in
`_identify_next_status_change!` on every event.
"""
function all_subtraps_filled(tstruct, ix, filled_traps)
    @inbounds for c in subtrapsof(tstruct, ix)
        filled_traps[c] || return false
    end
    return true
end

# ----------------------------------------------------------------------------
"""
    _identify_next_status_change!(changetimeest, cur_amounts, rateinfo, filled_traps, tstruct,
                                  z_vol_tables, cur_time, tmax, net_trap_set = nothing,
                                  net_covered_set = nothing) -> (earliest_changetime, fill_updates)

Find the next moment any trap flips fill status, and which traps flip at it.

Estimates are cheap and exact changetimes are not, so this is a branch and bound: take the
candidate with the smallest lower bound, resolve it exactly, and drop every candidate whose
lower bound now exceeds the best exact time, until none remain.  Ties are kept — several traps
can flip at the same instant.

Network nodes stay candidates and carry their exact ODE prediction; subsumed descendants are
excluded, and nodes skip the constant-rate subtrap recompute since the touch step re-predicts
them.

# Returns
`(earliest_changetime, fill_updates)` — the absolute time of the next change (`tmax` if none
falls before it), and one `IncrementalUpdate{Bool}` per trap flipping at it.  Empty updates mean
no further event this period.

Mutates `changetimeest`: resolved candidates are pinned to their exact time, traps that fired
are set to `Inf`, and the children of a fired trap are re-estimated.  `filled_traps` is flipped
temporarily while re-estimating those children (a child's changetime depends on its parent's
status) and restored before returning — so it is consistent on return, but not exception-safe.
"""
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
        changetimeest[cur_candidate] = ctest # ctest.min == ctest.max == exact changetime
        if ctest.min < earliest_changetime
            cur_best_candidates = [cur_candidate] # discard previous, better found
            earliest_changetime = ctest.min
        elseif ctest.min == earliest_changetime
            push!(cur_best_candidates, cur_candidate)
        end
        # drop the examined candidate and any that can no longer beat earliest_changetime
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

        # won't change again until the next weather event recomputes everything
        changetimeest[cand] = ChangeTimeEstimate(false, Inf, Inf)

        # network traps are re-predicted by the touch step, so the constant-rate subtrap
        # recompute below is both unnecessary and wrong for them
        is_node(cand) && continue

        # a subtrap's changetime depends on the parent's filled status, so recompute the
        # children with `cand` temporarily flipped (restored afterwards)
        filled_traps[cand] = !filled_traps[cand]
        children = subtrapsof(tstruct, cand)
        for c in children
            changetimeest[c] =
                _compute_changetime_estimate(c, cur_amounts, earliest_changetime,
                                             rateinfo, filled_traps, tstruct)
        end
        filled_traps[cand] = !filled_traps[cand]
    end
    return earliest_changetime, fill_updates
end

# ----------------------------------------------------------------------------
"""
    _compute_exact_fill(rateinfo, cur_amounts, trap, filled_traps, tstruct, time,
                        z_vol_tables, use_saved) -> Float64

The exact water volume in `trap` at absolute `time`, projected from `cur_amounts[trap]` under a
constant rate.

Short-circuits the projection for the two states whose volume is already known: a filled trap
sits at its own capacity, and a parent with any unfilled child is empty.  Otherwise the volume
comes from [`fill_trap_until`](@ref), with `use_saved` selecting the saved rather than the
current rates.  Nothing is mutated.
"""
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
"""
    _compute_exact_changetime(trap, changetimes, cur_amounts, rateinfo, tstruct, filled_traps,
                              z_vol_tables) -> ChangeTimeEstimate

Resolve `trap`'s bracketed changetime to an exact one, returned as a degenerate estimate
(`min == max`) at the absolute time it flips, or `Inf` if it never does.

A filled trap is resolved by projecting its *parent* — it can only stop being full once the
parent drains enough to expose it — and an unfilled one by projecting itself.  Either way the
projection is [`fill_trap_until`](@ref), whose elapsed result is placed on the absolute clock.

Nothing is mutated.  This code is dynamic network-agnostic: a network-covered trap arrives with
`min == max` already set from its network prediction, so the guard returns that straight back
and none of the constant-rate logic runs.
"""
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
"""
    _compute_z_vol_tables(tstruct) -> Vector{Tuple{Vector{Float64},Vector{Float64}}}

Per trap, a table mapping trap water volume to water level (z), for fast level lookups.

Returns one `(z, v)` pair per trap of `tstruct`, both ascending and of equal length, so a level
is read by interpolating `v -> z`.  `z` runs from the trap's bottom to its spillpoint
elevation, with duplicate levels dropped; a parent trap's bottom is raised to its subtraps'
spillpoint, so `v` counts only its own volume.  `tstruct` is not mutated.
"""
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
        # drop duplicate z levels, keeping the last (its volume is the correct one)
        keep = fill(true, length(zsorted))
        for i = 1:length(zsorted)-1
            keep[i] = zsorted[i+1] > zsorted[i]
        end
        zvt[tix] = (zsorted[keep], vsorted[keep])
    end
    return zvt
end

# ----------------------------------------------------------------------------
"""
    fill_trap_until(trap, rateinfo, cur_amount, endtime, tstruct, z_vol_tables;
                    use_saved=false) -> (amount, tstop)

Project `trap` forward under a constant inflow, filling or emptying it over
`[cur_amount.time, endtime]`.  With a fill-independent net rate the time to full/empty is
solved analytically; otherwise the wetted footprint (and so the infiltration loss) varies with
the fill level and an ODE is integrated instead.

# Arguments
- `trap`: trap index into `tstruct`.
- `rateinfo`: the rate state supplying the trap's inflow and infiltration bounds.
- `cur_amount`: the `FilledAmount` to start from — its `amount` and its `time`.
- `endtime`: absolute time to project to (`Inf` to run to the next full/empty).
- `tstruct`, `z_vol_tables`: terrain geometry and the cached volume->level tables.
- `use_saved`: read the *saved* rates rather than the current ones, projecting the trap under
  the inflow it has been evolving under rather than a newly-updated one.

# Returns
`(amount, tstop)`, mutating nothing.  `amount` is the water volume, between 0 and the trap's
own capacity.  `tstop` is the time it took to become full or empty — an **elapsed duration**,
relative to `cur_amount.time`, not an absolute time like `endtime`; add `cur_amount.time` to
place it on the absolute clock.  It is `nothing` if the trap became neither full nor empty
before `endtime`, and `0.0` if it was already full or empty on entry.
"""
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
        # net rate is fill-independent, so no ODE needed for the time to full/empty
        accum_rate = inflow - Smax
        (accum_rate == 0.0) && return (cur_amount.amount, nothing) # no change
        dt = (accum_rate > 0) ?
                (tvolume - cur_amount.amount) / accum_rate : # time to full
                cur_amount.amount / abs(accum_rate)          # time to empty
        reached = (cur_amount.time + dt <= endtime)

        return (reached) ?
            (accum_rate > 0.0 ? tvolume : 0.0, dt) :
            (cur_amount.amount + (endtime - cur_amount.time) * accum_rate, nothing)
    end
    # infiltration depends on fill level (wetted footprint varies), so solve an ODE
    fprint_infil =
        use_saved ? [-min(getsavedrunoff(rateinfo, i), 0.0) for i in footprint] :
                    -min.(getrunoff(rateinfo, footprint), 0.0)
    # Footprint cells at or above the spillpoint never pond, so drop their infiltration: this
    # keeps the wetted-area loss continuous at V = capacity and removes fill/unfill chatter
    # (same rule as `_ponding_infiltration` and the network solver).  Test the actual terrain
    # height, not the raised `trap_bottom`, which would over-null a zero-own-volume parent.
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
    (v0[1] == 0.0 && dv0[1] <= 0) && return (0.0, 0.0)
    (v0[1] == tvolume && dv0[1] >= 0) && return (tvolume, 0.0)
    
    function condition(out, v, t, integrator)
        out[1] = tvolume - v[1] # trap full
        out[2] = v[1]           # trap empty
        deriv = [0.0]
        dvdt(deriv, v, 0, t)
        out[3] = dv0[1] * deriv[1] # stagnation: net rate changed sign (passed through zero)
    end

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
"""
    _setup_dvdt(trap_bottom, trapvolume, infilfun, inflow, spillpoint, zvtable) -> Function

Build the ODE rate closure `_dvdt(dv, v, p, t)` for a single trap's water volume, for
[`fill_trap_until`](@ref) to integrate.

The rate is `inflow` minus infiltration through the *wetted* footprint only: the closure reads
the water level off `zvtable`, selects the cells at or below it, and sums their infiltration via
`infilfun`.  That makes it a step function of `v` — hence the ODE rather than a closed form —
but autonomous, ignoring `t`.  The returned closure writes `dv[1]` and mutates nothing else.
"""
function _setup_dvdt(trap_bottom, trapvolume, infilfun, inflow, spillpoint, zvtable)
    # note: for parent traps, trap_bottom will represent the bottom topography
    # except within the footprint of its children's traps, where trap_bottom
    # will have the value of the upper child's spillpoint

    zmin = minimum(trap_bottom)
    zspan = spillpoint.elevation - zmin

    v2z = length(zvtable[2]) == 1 ?
        x -> zmin : # degenerate case
        Interpolations.linear_interpolation(zvtable[2], zvtable[1],
                                             extrapolation_bc=Interpolations.Line())
    function _dvdt(dv, v, p, t)
        z = (v[1] <= 0)          ? zmin :
            (v[1] >= trapvolume) ? Inf :
                                  v2z(v[1])
        return dv[1] = inflow - infilfun(trap_bottom .<= z)
    end

    return _dvdt
end
