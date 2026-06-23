# Visual check of `culvert_rate`: sweep the boundary conditions of a couple of
# culverts and plot their flow capacity, so the weir->orifice switch, the
# inlet/outlet-control crossover, and reverse flow can be eyeballed.
#
# Run from the REPL (needs a display; shares the examples/ project):
#
#   using Pkg
#   Pkg.activate("examples")
#   include("examples/verification/culvert_rate.jl")
#   plot_culvert_rates()

using SurfaceWaterIntegratedModeling
import GLMakie

# `culvert_rate` only reads `tstruct.topography`, so a NamedTuple is a fine stand-in.
# A 1x2 grid with a small drop puts the inlet invert just above the outlet invert.
# The drop is kept modest (1 m) on purpose: with a steep drop inlet control governs
# regardless of tailwater and the two outlet conditions would plot identically.
const _DROP_TS = (; topography = [1.0 0.0])   # z(inlet)=1 m, z(outlet)=0 m
const _INLET   = CartesianIndex(1, 1)
const _OUTLET  = CartesianIndex(1, 2)
const _TW_HIGH = 3.0                           # high tailwater for the submerged case

# Sweep inlet head and return the rate for one outlet condition.  The inlet is
# treated as submerged once the head reaches the barrel diameter (mimicking the
# event the solver would raise), which is where weir flow gives way to orifice.
function _rate_sweep(cv; outlet_submerged, outlet_head, heads, allow_reverse = false)
    D = 2 * cv.r
    return [culvert_rate(cv, _DROP_TS;
                         inlet_submerged  = h >= D, inlet_head = h,
                         outlet_submerged = outlet_submerged, outlet_head = outlet_head,
                         allow_reverse = allow_reverse)
            for h in heads]
end

"""
    plot_culvert_rates() -> GLMakie.Figure

Plot, for two contrasting culverts (columns):
  * top row    -- capacity vs inlet head under a free outfall and a submerged
                  outlet (downhill-only, the default `allow_reverse = false`);
  * bottom row -- the submerged-outlet case with `allow_reverse` off vs on, so
                  the negative (outlet->inlet) reverse branch in the drowned
                  region is visible.
Returns the figure.
"""
function plot_culvert_rates(; heads = range(0.0, 4.0; length = 400))
    # Two contrasting culverts (built with the convenience constructor):
    #   A: small, smooth concrete pipe.   B: larger, rough corrugated-metal pipe.
    cvA = DynCulvert(_DROP_TS, _INLET, _OUTLET; r = 0.5,  n = 0.012)
    cvB = DynCulvert(_DROP_TS, _INLET, _OUTLET; r = 0.75, n = 0.024)

    fig = GLMakie.Figure(size = (960, 800))
    for (col, (cv, name)) in enumerate([(cvA, "A: r=0.5 m, smooth"),
                                        (cvB, "B: r=0.75 m, rough")])
        D = 2 * cv.r

        # --- top: outlet boundary conditions, downhill-only -------------------
        ax1 = GLMakie.Axis(fig[1, col];
                           title  = "Culvert $name",
                           xlabel = "inlet head (m)", ylabel = "capacity (m³/s)")
        GLMakie.lines!(ax1, heads, _rate_sweep(cv; outlet_submerged = false,
                       outlet_head = 0.0, heads = heads);
                       label = "free outfall", linewidth = 2)
        GLMakie.lines!(ax1, heads, _rate_sweep(cv; outlet_submerged = true,
                       outlet_head = _TW_HIGH, heads = heads);
                       label = "submerged outlet (TW = $(_TW_HIGH) m)", linewidth = 2)
        GLMakie.vlines!(ax1, [D]; color = (:black, 0.3), linestyle = :dash)  # weir->orifice
        GLMakie.axislegend(ax1; position = :lt)

        # --- bottom: reverse flow in the drowned region -----------------------
        ax2 = GLMakie.Axis(fig[2, col];
                           title  = "submerged outlet (TW = $(_TW_HIGH) m): reverse flow",
                           xlabel = "inlet head (m)", ylabel = "rate (m³/s)")
        GLMakie.lines!(ax2, heads, _rate_sweep(cv; outlet_submerged = true,
                       outlet_head = _TW_HIGH, heads = heads, allow_reverse = false);
                       label = "allow_reverse = false", linewidth = 2)
        GLMakie.lines!(ax2, heads, _rate_sweep(cv; outlet_submerged = true,
                       outlet_head = _TW_HIGH, heads = heads, allow_reverse = true);
                       label = "allow_reverse = true", linewidth = 2)
        GLMakie.hlines!(ax2, [0.0]; color = (:black, 0.4), linestyle = :dash)
        GLMakie.axislegend(ax2; position = :lt)
    end

    GLMakie.display(fig)
    return fig
end
