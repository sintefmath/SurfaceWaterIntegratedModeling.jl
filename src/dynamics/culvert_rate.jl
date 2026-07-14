# Culvert flow-rate hydraulics.
#
# Computes the volumetric flow (m^3/s, SI) through a `DynCulvert` given the water
# level at each end.  Simplified FHWA HDS-5 model:
#
#   * Inlet control: a hard switch on the upstream end's submergence -- weir below,
#     orifice once submerged.  No weir<->orifice transition blend.
#   * Outlet control: full-barrel pressure flow.  The tailwater is the downstream
#     pool surface when that end holds water, or a free-outfall critical-depth
#     tailwater when it is dry.
#   * The governing (more restrictive, i.e. smaller) of the two rates is returned.
#
# Flow runs from the higher water surface to the lower.  A culvert is hydraulically
# symmetric, so reverse flow (outlet pool higher than inlet) is the same model with
# the ends swapped, returned as a negative rate -- enabled by `allow_reverse`.
#
# All heads are water level above the inlet/outlet cell elevation (the invert).

export culvert_rate

# Non-negative capacity for a single flow direction, from the upstream end to the
# downstream end.  Each end is (submerged?, head above invert, invert elevation),
# but only the *upstream* submergence matters -- it selects weir vs orifice at the
# entrance.  The entrance-loss coefficient cv.Ke is used symmetrically.
function _directional_capacity(cv::DynCulvert,
                               up_submerged::Bool, up_head, up_z,
                               down_head, down_z)
    g = 9.81
    D = 2 * cv.r              # barrel diameter
    A = pi * cv.r^2           # full cross-sectional area

    # --- inlet control at the upstream end: weir below submergence, orifice above
    # (hard switch on submergence; no transition blend; weir width = D).
    Q_inlet = up_submerged ?
        cv.Cd * A * sqrt(2g * up_head) :    # orifice (submerged entrance)
        cv.Cw * D * up_head^1.5             # weir    (free entrance)

    # --- outlet control (full-barrel pressure flow) ---------------------------
    up_elev = up_z + up_head                # upstream water-surface elevation
    if down_head > 0
        down_elev = down_z + down_head      # real downstream pool surface
    else
        # Downstream end is dry -> free outfall; effective tailwater from the
        # critical depth at the outfall.
        # @@@ Approximate: yc uses the rough HDS-5 §5 formula, evaluated one-shot
        #     at Q_inlet (no Q<->yc iteration).  Acceptable because this branch
        #     rarely governs; use the exact iterative circular solution if more
        #     free-outfall accuracy is ever needed.
        yc = D * 0.325 * (Q_inlet / D^2.5)^(2/3) + 0.083 * D
        down_elev = down_z + max(yc, (yc + D) / 2)
    end
    dH = max(up_elev - down_elev, 0.0)      # net driving head (0 if not driving)
    Q_outlet = A * sqrt(2g * dH) / sqrt(1 + cv.Ke + cv.Kf)

    # The slope/Froude applicability gate (HDS-5 §3) is intentionally omitted: the
    # min() below already lets inlet control win when the barrel cannot fill.
    return min(Q_inlet, Q_outlet)
end

"""
    culvert_rate(cv;
                 inlet_submerged,  inlet_head,  inlet_invert,
                 outlet_submerged, outlet_head, outlet_invert,
                 allow_reverse = false) -> Float64

Volumetric flow through culvert `cv` in m^3/s (SI), given the submergence state, water head
(metres above the invert), and invert elevation at each end.

A positive result is flow from the inlet to the outlet.  By default
(`allow_reverse = false`) the culvert is downhill-only: if the outlet pool stands
higher than the inlet pool the forward driving head is zero and `0` is returned
("drowned").  With `allow_reverse = true`, when the outlet surface is the higher of
the two the flow reverses and a **negative** (outlet->inlet) rate is returned.

In the chosen direction the rate is `min(Q_inlet_control, Q_outlet_control)` -- the
bottleneck that physically governs.  See the module header for the assumptions.

!!! note
    @@@ The network solver does not yet consume a negative (reverse) rate -- the
    routing layer still assumes downhill flow (see `CLAUDE.md`/`AGENTS.md`).  Use
    `allow_reverse = true` only where the caller is prepared to handle it.
"""
function culvert_rate(cv::DynCulvert;
                      inlet_submerged::Bool,  inlet_head::Real,  inlet_invert::Real,
                      outlet_submerged::Bool, outlet_head::Real, outlet_invert::Real,
                      allow_reverse::Bool = false)
    hi = max(float(inlet_head),  0.0)   # head above inlet invert
    ho = max(float(outlet_head), 0.0)   # head above outlet invert
    z_in  = float(inlet_invert)
    z_out = float(outlet_invert)

    # Flow runs from the higher water surface to the lower one.  The downhill-only
    # default always treats the inlet as upstream; if the outlet pool is actually
    # higher the forward dH clamps to 0 and the culvert reads as drowned.  With
    # allow_reverse, an outlet surface above the inlet's genuinely reverses the
    # flow: compute it directly in the outlet->inlet direction and negate it.
    if allow_reverse && (z_out + ho) > (z_in + hi)
        return -_directional_capacity(cv, outlet_submerged, ho, z_out, hi, z_in)
    end
    return _directional_capacity(cv, inlet_submerged, hi, z_in, ho, z_out)
end
