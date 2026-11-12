# Rebuild the SUrbArea webapp "synthetic" NBS case outside the webapp, using the
# current SWIM checkout.  Mirrors ExplorersController: _compute_tstruct (no layers),
# _build_nbs_placements / _nbs_feature_gridgeom, _parse_precip_events, _simulate_sequence.
#
# Run with the WebInterface environment:
#   julia --project=<SUrbArea>/WebInterface run1.jl      (then run2.jl)
#
# The SUrbArea checkout is assumed to sit next to this repository; override with the
# environment variable SURBAREA_DIR.  Inputs come from <SUrbArea>/debug/ (not in git).

using SurfaceWaterIntegratedModeling
const SWIM = SurfaceWaterIntegratedModeling
using SUrbArea, JLD2, JSON

const SURB  = get(ENV, "SURBAREA_DIR", normpath(joinpath(@__DIR__, "..", "..", "..", "..", "SUrbArea")))
const BBOX  = [0.0, 0.0, 200.0, 200.0]

load_dtm_flat() = JLD2.load(joinpath(SURB, "debug", "synthetic_dtm_flattened.jld2"), "dtm")

function load_features(file)
    d = JSON.parsefile(file)
    return d["features"]
end

# GreenRoof only (the synthetic setups use only that model)
function build_placements(features, gridsize)
    pl = DynNBSPlacement[]
    id = 0
    for f in features
        p = f["properties"]
        get(p, "objekttypenavn", "") == "GreenRoof" || continue
        ring = f["geometry"]["coordinates"][1]
        mask = zeros(Int, gridsize)
        ring_grid, _ = SUrbArea.geocoord2gridix(ring, BBOX, gridsize)
        SUrbArea._fill_polygon!(mask, ring_grid, 1)
        footprint = [LinearIndices(gridsize)[c] for c in findall(!iszero, mask)]
        gix, inside = SUrbArea.geocoord2gridix([p["outlet"]], BBOX, gridsize)
        @assert inside[1]
        outlet = CartesianIndex(gix[1][1], gix[1][2])
        num(k, d) = Float64(get(p, k, d))
        sys = elhadiGreenRoof(num("Smin_soil", 10.0), num("Smin_drain", 5.0),
                              num("kINF", 0.5), num("kDRAIN", 0.5),
                              num("nINF", 1.0), num("nDRAIN", 1.0),
                              num("EVC", 0.0), num("EVS11", 1.0))
        id += 1
        np = DynNBSPlacement(sys, footprint, [outlet])
        np.id = id
        push!(pl, np)
    end
    return pl
end

function gauges(features, gridsize)
    g = CartesianIndex{2}[]
    for f in features
        get(f["properties"], "objekttypenavn", "") == "Gaugepoint" || continue
        gix, inside = SUrbArea.geocoord2gridix([f["geometry"]["coordinates"]], BBOX, gridsize)
        inside[1] && push!(g, CartesianIndex(gix[1][1], gix[1][2]))
    end
    return g
end

# mm/hr -> m/min, as _parse_precip_events
weather(events) = begin
    t = 0.0; ev = WeatherEvent[]
    for (mmh, dur) in events
        push!(ev, WeatherEvent(t, mmh / (1000.0 * 60.0))); t += dur
    end
    ev
end
