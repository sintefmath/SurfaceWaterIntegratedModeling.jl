export RateInfo, getrunoff, getsmin, getsmax, getinflow, setrunoff!, setsmin!, setsmax!,
setinflow!, setsavepoint!, disablesavepoint!, getsavedrunoff, getsavedsmin, getsavedsmax,
getsavedinflow, getrunoffupdates, getinflowupdates, copy

"""
RateInfo

A mutable struct representing rates of water (inflow, infiltration) associated
with the terrain and its traps.

The fields of this struct should in general not be accessed directly, but through
the set of associated methods, listed below. 

Fields
- `runoff::Matrix{Float64}`: net overland flow or remaining infiltration capacity, 
                             one value per gridcell in the terrain grid.
- `Smax::Vector{Float64}`: maximum remaining infiltration within trap footprint, one
                           value per trap.
- `Smin::Vector{Float64}`: minimum infiltration within trap footprint, one value
                           per trap. (Only nonzero for traps with subtraps, and then
                            equals the sum of `Smax` for its subtraps).
- `trap_inflow::Vector{Float64}`: total inflow to each trap.

- `stored_runoff_values::Dict{Int,Float64}`: if values have changed and `save_active`
                                             is set to `true`, this field stores the
                                             original runoff values for those that 
                                             were later modified.
- `stored_Smin_values::Dict{Int,Float64}`: if values have changed and `save_active`
                                             is set to `true`, this field stores the
                                             original Sminvalues for those that 
                                             were later modified.
                                           
- `stored_Smax_values::Dict{Int,Float64}`: if values have changed and `save_active`
                                             is set to `true`, this field stores the
                                             original Smax values for those that 
                                             were later modified.
                                           
- `stored_inflow_values::Dict{Int,Float64}`: if values have changed and `save_active`
                                             is set to `true`, this field stores the
                                             original Smax values for those that 
                                             were later modified.
- `save_active::Bool`: flag specifying if the initial state should be preseved.

Associated methods for access
-  `getrunoff`
-  `getsmin`
-  `getsmax`
-  `getinflow`
-  `setrunoff!`
-  `setsmin!`
-  `setsmax!`
-  `setinflow!`
-  `setsavepoint!`
-  `disablesavepoint!`
-  `getsavedrunoff`
-  `getsavedsmin`
-  `getsavedsmax`
-  `getsavedinflow`
-  `getrunoffupdates`
-  `getinflowupdates` 

## Constructor
`RateInfo(runoff, Smax, Smin, trap_inflow)`

Creates a new RateInfo object with the provided parameters.
"""
mutable struct RateInfo
    runoff::Matrix{Float64} # net overland flow (one per region, may be modified
                            # as traps spill over), or remaining infiltration
                            # capacity (negative values)
    Smax::Vector{Float64} # maximum remaining infiltration within trap footprint
                          # (one per trap)
    Smin::Vector{Float64} # minimum infiltration within trap footprint (one per
                          # trap; only nonzero for parent traps)
    trap_inflow::Vector{Float64} # total inflow to each trap

    # the RateInfo object can be set to store its state at a specific time.
    # Since the terrain can be big, we only store the incremental differences
    # between current and saved state
    stored_runoff_values::Dict{Int,Float64}
    stored_Smin_values::Dict{Int, Float64}
    stored_Smax_values::Dict{Int, Float64}
    stored_inflow_values::Dict{Int, Float64}
    save_active::Bool

    function RateInfo(runoff, Smax, Smin, trap_inflow)
        new(runoff, Smax, Smin, trap_inflow, Dict{Int, Float64}(),
            Dict{Int, Float64}(), Dict{Int, Float64}(), Dict{Int, Float64}(), false)
    end
end

"""
Get the current runoff value for the gridcell with linear index 'ix'.
"""
function getrunoff(ri::RateInfo, ix)
    return ri.runoff[ix]
end

"""
Get the current Smin value for the trap with index 'ix'.
"""
function getsmin(ri::RateInfo, ix)
    return ri.Smin[ix]
end

"""
Get the current Smax value for the trap with index 'ix'.
"""
function getsmax(ri::RateInfo, ix)
    return ri.Smax[ix]
end

"""
Get the current inflow value for the trap with index 'ix'.
"""
function getinflow(ri::RateInfo, ix)
    return ri.trap_inflow[ix]
end

"""
Set the current inflow value for the trap with index 'ix'.
"""
function setinflow!(ri::RateInfo, ix, val)
    if ri.save_active && !haskey(ri.stored_inflow_values, ix)
        ri.stored_inflow_values[ix] = ri.trap_inflow[ix]
    end
    ri.trap_inflow[ix] = val
end

"""
Set the current runoff value for the gridcell with linear index 'ix'.
"""
function setrunoff!(ri::RateInfo, ix, val)
    if ri.save_active && !haskey(ri.stored_runoff_values, ix)
        ri.stored_runoff_values[ix] = ri.runoff[ix]
    end
    ri.runoff[ix] = val
end

"""
Set the current Smin value for the trap with index 'ix'.
"""
function setsmin!(ri::RateInfo, ix, val)
    if ri.save_active && !haskey(ri.stored_Smin_values, ix)
        ri.stored_Smin_values[ix] = ri.Smin[ix]
    end
    ri.Smin[ix] = val
end

"""
Set the current Smax value for the trap with index 'ix'.
"""
function setsmax!(ri::RateInfo, ix, val)
    if ri.save_active && !haskey(ri.stored_Smax_values, ix)
        ri.stored_Smax_values[ix] = ri.Smax[ix]
    end
    ri.Smax[ix] = val
end

"""
Set the savepoint here. (State at this point will be remembered, even
in the face of later updates, until the savepoint is reset or disabled).
"""
function setsavepoint!(ri::RateInfo)
    _clear_saved_dicts!(ri)
    ri.save_active = true
end

function _clear_saved_dicts!(ri::RateInfo)
    ri.stored_runoff_values = Dict{Int, Float64}()
    ri.stored_Smin_values = Dict{Int, Float64}()
    ri.stored_Smax_values = Dict{Int, Float64}()
    ri.stored_inflow_values = Dict{Int, Float64}()
end
    
"""
Disable the current savepoint.  Saved state will no longer be maintained.
"""
function disablesavepoint!(ri::RateInfo)
    _clear_saved_dicts!(ri)
    ri.save_active = false
end

"""
Get the saved runoff value for the gridcell with linear index 'ix'.
"""
function getsavedrunoff(ri::RateInfo, ix)
    return get(ri.stored_runoff_values, ix, ri.runoff[ix])
end

"""
Get the saved Smin value for the trap with index 'ix'.
"""
function getsavedsmin(ri::RateInfo, ix)
    return get(ri.stored_Smin_values, ix, ri.Smin[ix])
end

"""
Get the saved Smax value for the trap with index 'ix'.
"""
function getsavedsmax(ri::RateInfo, ix)
    return get(ri.stored_Smax_values, ix, ri.Smax[ix])
end

"""
Get the saved inflow value for the trap with index 'ix'.
"""
function getsavedinflow(ri::RateInfo, ix)
    return get(ri.stored_inflow_values, ix, ri.trap_inflow[ix])
end

"""
Get a vector with the runoff updates that have occured since last savepoint was set.
"""
function getrunoffupdates(ri::RateInfo)
    return [IncrementalUpdate(k, getrunoff(ri, k))
            for k in keys(ri.stored_runoff_values)]
end

"""
Get a vector with the inflow updates that have occured since the last savepoint was set.
"""
function getinflowupdates(ri::RateInfo)
    return [IncrementalUpdate(k, getinflow(ri, k))
            for k in keys(ri.stored_inflow_values)]
end

"""
A copy of a RateInfo object is always a deep copy.
"""
Base.copy(ri::RateInfo) = RateInfo(deepcopy(ri.runoff),
                                   deepcopy(ri.Smax),
                                   deepcopy(ri.Smin),
                                   deepcopy(ri.stored_runoff_values),
                                   deepcopy(ri.stored_Smin_values),
                                   deepcopy(ri.stored_Smax_values),
                                   deepcopy(ri.stored_inflow_values),
                                   ri.save_active)
