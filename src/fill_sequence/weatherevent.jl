export WeatherEvent

"""
    WeatherEvent

A struct representing a weather event, i.e. a change of weather.

# Fields
- `timestamp::Float64`: Start time point of the weather event

- `rain_rate::Union{Matrix{Float64}, Float64}`: 
      Rain rate set by the event.  It can be given either as a single
      floating-point number that represents a uniform rate across the 
      terrain, or a matrix of floating-point number for individual rain rates
      per grid cell.

The unit of rain rate is given in _topographical grid height unit_ per _time unit_.
E.g., if the topography is presented in a height unit of meters, and time is 
counted in hours, then the specified rate value will be interpreted as meters 
per hour.

See also [`fill_sequence`](@ref), which takes a `Vector{WeatherEvent}` as one
of its inputs.
"""
struct WeatherEvent
    # Start time point of this weather event
    timestamp::Float64

    # Rain rate set by this weather event.  Either a single 'Float64' for a
    # uniform rate across terrain, or a matrix of 'Float64', giving an
    # individual rate for each grid cell in the terrain.  Rate is given in terms
    # of millimeters per time unit.
    rain_rate::Union{Matrix{Float64}, Float64} 
end
