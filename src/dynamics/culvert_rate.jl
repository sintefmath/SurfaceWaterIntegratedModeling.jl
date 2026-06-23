# Culvert flow-rate hydraulics.
#
# Computes the volumetric capacity (m^3/s, SI) of a `DynCulvert` given the water
# level at each end.  The model is a deliberately simplified version of the
# FHWA HDS-5 procedure (see agent/reports/culvert_hydraulics_reference.md and the
# "Clarifications" heading in agent/prompts/culvert_rate_implementation.org):
#
#   * Inlet control: a hard switch on inlet submergence -- weir below, orifice
#     once submerged.  No weir<->orifice transition blend.
#   * Outlet control: full-barrel pressure flow.  When the downstream end is not
#     submerged a free-outfall effective tailwater is assumed.
#   * The governing (more restrictive, i.e. smaller) of the two rates is returned.
#
# A culvert is hydraulically symmetric, so reverse flow (downstream pool higher
# than upstream) is the same model with the two ends swapped and the sign flipped
# -- see `_directional_capacity` and the `allow_reverse` option below.
#
# All heads are measured as the water level above the inlet/outlet cell elevation
# (the invert) -- no diameter/centroid correction.

export culvert_rate

# Non-negative capacity for a single flow direction, from the `up`(stream) end to
# the `down`(stream) end.  Each end is given as (submerged?, head above invert,
# invert elevation).  The entrance-loss coefficient `cv.Ke` is used symmetrically
# (same value whichever end is the entrance).
function _directional_capacity(cv::DynCulvert, g, D, A,
                               up_submerged::Bool,   up_head,   up_z,
                               down_submerged::Bool, down_head, down_z)
    # --- inlet control at the upstream end ------------------------------------
    # Hard switch on submergence (no transition blend); weir width = D.
    Q_inlet = up_submerged ?
        cv.Cd * A * sqrt(2g * up_head) :    # orifice (submerged entrance)
        cv.Cw * D * up_head^1.5             # weir    (free entrance)

    # --- outlet control (full-barrel pressure flow) ---------------------------
    up_elev = up_z + up_head                # upstream water-surface elevation
    if down_submerged
        down_elev = down_z + down_head      # tailwater surface elevation
    else
        # Free outfall: no real tailwater, so use an effective one from the
        # critical depth at the downstream end.
        # @@@ yc evaluated one-shot at Q_inlet (no Q<->yc iteration); the
        #     free-outfall outlet-control branch almost never governs anyway.
        # @@@ rough critical-depth approximation (HDS-5 §5); exact iterative
        #     circular solution is available if more accuracy is ever needed.
        yc = D * 0.325 * (Q_inlet / D^2.5)^(2/3) + 0.083 * D
        down_elev = down_z + max(yc, (yc + D) / 2)
    end
    dH = max(up_elev - down_elev, 0.0)      # net driving head (0 if not driving)
    Q_outlet = A * sqrt(2g * dH) / sqrt(1 + cv.Ke + cv.Kf)

    # @@@ slope/Froude applicability gate (HDS-5 §3) intentionally omitted: the
    #     min() below already lets inlet control win when the barrel cannot fill.
    return min(Q_inlet, Q_outlet)
end

"""
    culvert_rate(cv, tstruct;
                 inlet_submerged,  inlet_head,
                 outlet_submerged, outlet_head,
                 allow_reverse = false) -> Float64

Volumetric flow through culvert `cv` in m^3/s (SI), given the submergence state
and water head (metres above the cell invert) at each end.  `tstruct` supplies
the inlet/outlet invert elevations used by the outlet-control balance.

A positive result is flow from the inlet to the outlet.  By default
(`allow_reverse = false`) the culvert is treated as downhill-only: if the outlet
pool stands higher than the inlet pool the culvert is "drowned" and `0` is
returned.  With `allow_reverse = true` the reverse direction is computed by
swapping the two ends (the model is symmetric, with `Ke` reused symmetrically)
and returned as a **negative** rate.

Within each direction the returned rate is `min(Q_inlet_control, Q_outlet_control)`
-- the bottleneck that physically governs.  See the module header for the
modelling assumptions.

!!! note
    @@@ The network solver does not yet consume a negative (reverse) rate -- the
    routing layer still assumes downhill flow (see `CLAUDE.md`/`AGENTS.md`).  Use
    `allow_reverse = true` only where the caller is prepared to handle it.
"""
function culvert_rate(cv::DynCulvert, tstruct;
                      inlet_submerged::Bool,  inlet_head::Real,
                      outlet_submerged::Bool, outlet_head::Real,
                      allow_reverse::Bool = false)
    g = 9.81
    D = 2 * cv.r              # barrel diameter
    A = pi * cv.r^2           # full cross-sectional area

    hi = max(float(inlet_head),  0.0)   # head above inlet invert
    ho = max(float(outlet_head), 0.0)   # head above outlet invert
    z_in  = tstruct.topography[cv.inlet]
    z_out = tstruct.topography[cv.outlet]

    # Forward direction: inlet (upstream) -> outlet (downstream).
    Q_fwd = _directional_capacity(cv, g, D, A, inlet_submerged, hi, z_in,
                                                outlet_submerged, ho, z_out)
    allow_reverse || return Q_fwd       # downhill-only: never reverse

    # Reverse direction: outlet (upstream) -> inlet (downstream).  At most one of
    # the two directions has a positive driving head (the other's outlet-control
    # head clamps to 0), so the signed difference is the net flow: + inlet->outlet,
    # - outlet->inlet, and it passes continuously through 0 at the crossover.
    Q_rev = _directional_capacity(cv, g, D, A, outlet_submerged, ho, z_out,
                                                inlet_submerged, hi, z_in)
    return Q_fwd - Q_rev
end
