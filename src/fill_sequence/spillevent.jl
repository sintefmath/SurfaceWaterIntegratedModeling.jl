export IncrementalUpdate, SpillEvent, ChangeTimeEstimate, FilledAmount, amount_at,
       filled_at, inflow_at, rainrate_at, runoff_at

"""
Structure used when estimating when a [`SpillEvent`](@ref) will occur, i.e.
when a given non-filled trap will completely fill up, or a given filled trap 
will start emptying.  

The structure provides an upper and lower bound for when this event should happen.
"""
struct ChangeTimeEstimate
    filling::Bool # indicate if trap is filling or emptying
    min::Float64 # earliest time the trap will fill / empty below subtrap level
    max::Float64 # lastest time the trap will fill / empty below subtrap level
end

"""
Structure specifying that a trap will contain an `amount` volume of water at 
timepoint `time`.
"""
struct FilledAmount
    amount::Float64
    time::Float64
end

"""
Structure specifying an incremental update in the quantity of something.  It is 
used in several different settings.  It contain an index (to refer to a specific 
element to be updated) and a corresponding value (the value to update).
"""
struct IncrementalUpdate{ValType}
    index::Int
    value::ValType
end

"""
    struct SpillEvent

A struct representing a "spill event", i.e. a moment in time when a trap changes its
fill status.  A trap can either transitory or filled, and the spill event refers to the
moment in time when the trap passes from one to the other.

# Fields
- `timestamp::Float64`: the time when the event occurs

- `amount::Union{Vector{FilledAmount}, Vector{IncrementalUpdate{FilledAmount}}}`: 
   amount of water in the traps.  Can be presented as a vector with one entry per
   trap, or as a vector of incremental updates to certain traps (as compared to the
   last registered SpillEvent).

- `filled::Union{Vector{Bool}, Vector{IncrementalUpdate{Bool}}}`: 
      the traps that are filled when the event occurs.  Can be presented as a vector 
      with one (boolean) value per trap, or as an incremental update representing
      the change since last registered SpillEvent.

- `inflow::Union{Vector{Float64}, Vector{IncrementalUpdate{Float64}}}`: 
      water inflow to each trap, presented either as a vector with the individual
      values of each and every trap, or as a set of incremental updates representing
      the change since the last registered SpillEvent.  Inflow is a combination of
      direct presipitation, runoff from the trap's watershed or spillover from 
      upstream overflowing traps.

- `rain_rate::Union{Matrix{Float64}, Float64, Nothing}`: 
      if this event consisted of a weather change, it will be represented here, 
      either as a single floating point value (considered constant over the grid),
      or as a grid with individual values for each cell.

- `runoff::Union{Matrix{Float64}, Vector{IncrementalUpdate{Float64}}}`: 
      matrix keeping track of current overland flow and infiltration rate for each
      gridcell.  Positive values represent overland flow, negative values infiltration.
      Can be presented either as a matrix with values for each and every grid cell,
      or as a set of incremental changes compared with last SpillEvent.
"""
struct SpillEvent
    timestamp::Float64

    # keep track of amount of water in traps
    amount::Union{Vector{FilledAmount},
                  Vector{IncrementalUpdate{FilledAmount}}}

    # keep track of filled traps
    filled::Union{Vector{Bool},
                  Vector{IncrementalUpdate{Bool}}}

    # keep track of water inflow to each trap
    inflow::Union{Vector{Float64},
                  Vector{IncrementalUpdate{Float64}}}

    # keep track of changes in weather and runoff/infiltration
    rain_rate::Union{Matrix{Float64}, Float64, Nothing}

    runoff::Union{Matrix{Float64}, Vector{IncrementalUpdate{Float64}}}
end

# ----------------------------------------------------------------------------
"""
   amount\\_at(seq, end_ix=-1)

Given a sequence of chronological SpillEvents `seq`, determine the `FilledAmount`
for each trap at the end of the sequence (end\\_ix=-1), or at an earlier 
event (end\\_ix ∈ [1:length(seq)]).  

This function computes its result by running through the sequence and accumulating
differences until it has identified the correct entry for all traps.

Note that the returned `FilledAmount` for each trap does not necessarily
correspond to the timestamp of the corresponding `FillEvent`, but may represent
a value computed for an earlier point in time.  The `RateInfo` for that trap 
should however have remained constant since the last time the value was computed,
so it should be relatively straightforward to compute the present value if needed.
"""
function amount_at(seq::Vector{SpillEvent}, end_ix::Int=-1)
    return _reconstruct(seq, :amount, end_ix)
end

# ----------------------------------------------------------------------------
"""
   filled\\_at(seq, end_ix=-1)

Given a sequence of chronological SpillEvents `seq`, determine which traps are 
filled at the end of the sequence (end\\_ix=-1), or at an earlier 
event (end\\_ix ∈ [1:length(seq)]).  

This function computes its result by running through the sequence and accumulating
differences until it has identified the correct entry for all traps.
"""
function filled_at(seq::Vector{SpillEvent}, end_ix::Int=-1)
    return _reconstruct(seq, :filled, end_ix)
end

# ----------------------------------------------------------------------------
"""
   inflow\\_at(seq, end_ix=-1)

Given a sequence of chronological SpillEvents `seq`, determine the rate of water
flowing into a trap (from its watershed, from upstream, or directly from 
precipitation) at the end of the sequence (end\\_ix=-1), or at an earlier 
event (end\\_ix ∈ [1:length(seq)]).  

This function computes its result by running through the sequence and accumulating
differences until it has identified the correct entry for all traps.
"""
function inflow_at(seq::Vector{SpillEvent}, end_ix::Int=-1)
    return _reconstruct(seq, :inflow, end_ix)
end

# ----------------------------------------------------------------------------
"""
   runoff\\_at(seq, end_ix=-1)

Given a sequence of chronological SpillEvents `seq`, determine the runoff 
across each point of the terrain (i.e. each gridcell) at the end of the sequence 
(end\\_ix=-1), or at an earlier event (end\\_ix ∈ [1:length(seq)]).  

Positive values represent infiltration excess flow (runoff), whereas negative
values represent remaining infiltration capacity.

This function computes its result by running through the sequence and accumulating
differences until it has identified the correct entry for all traps.
"""
function runoff_at(seq::Vector{SpillEvent}, end_ix::Int=-1)
    return _reconstruct(seq, :runoff, end_ix)
end

# ----------------------------------------------------------------------------
"""
    rainrate\\_at(seq, end_ix=-1)

Given a sequence of chronological SpillEvents `seq`, determine the rain rate
in effect at the end of the sequence (end\\_ix=-1), or at an earlier event
 (end\\_ix ∈ [1:length(seq)]).  

This function computes its result by running through the sequence and identifying
the latest rainrate that went into effect before the sequence entry specified by
`end_ix`.
"""
function rainrate_at(seq::Vector{SpillEvent}, end_ix::Int=-1)

    ix = (end_ix == -1) ? length(seq) : end_ix
    while isnothing(getfield(seq[ix], :rain_rate))
        ix -= 1
        @assert ix > 0 "Invalid sequence. Missing requested field info."
    end
    return getfield(seq[ix], :rain_rate)
end

# ----------------------------------------------------------------------------
function _reconstruct(seq::Vector{SpillEvent}, fname, end_ix)

    end_ix = (end_ix == -1) ? length(seq) : end_ix
    
    ix = end_ix
    while eltype(getfield(seq[ix], fname)) <: IncrementalUpdate
        ix -= 1
        @assert ix > 0 "Invalid sequence.  Requested field not fully computed."
    end

    result = copy(getfield(seq[ix], fname))

    for jx = ix+1:end_ix
        updates = getfield(seq[jx], fname) # should be a Vector{Tuple{Int, type}}
        for u in updates
            result[u.index] = u.value
        end
    end
    return result
end
