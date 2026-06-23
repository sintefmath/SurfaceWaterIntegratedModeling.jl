# Culvert Hydraulics Reference

A self-contained specification for computing culvert capacity / flow rate given heads
on each side. Written to be served to an implementing agent (e.g. Claude Code) alongside
existing project scaffolding. Authoritative source for production work is **FHWA HDS-5**
(*Hydraulic Design of Highway Culverts*, 3rd ed.) and the free **HY-8** software; the
simplified forms below are conceptual limiting cases and are documented as such.

> Units note: friction-loss constants below are written for **English (US customary)**
> units (ft, ft³/s, g = 32.2 ft/s²). SI equivalents are flagged inline. Pick one unit
> system and keep it consistent throughout the implementation.

---

## 1. Core Concept

A culvert's capacity is set by whichever of two conditions is more restrictive:

- **Inlet control** — the entrance is the bottleneck. The barrel could carry more than
  the inlet will admit. Governed by inlet geometry, headwater depth, and culvert shape.
  Barrel length, roughness, slope, and tailwater do *not* matter.
- **Outlet control** — the barrel + outlet are the bottleneck. Governed by an energy
  balance across the whole barrel: friction, entrance loss, exit loss, tailwater.

**Procedure:** compute headwater (HW) under *both* conditions for the flow of interest;
the **higher HW (= lower capacity) governs**. Equivalently, for a given HW, the lower of
the two Q values governs.

Physical tendencies:
- Inlet control tends to govern for **short, steep** culverts or **poor entrances**
  (projecting pipe, no headwall).
- Outlet control tends to govern for **long, flat, rough** culverts or **high tailwater**.

---

## 2. Inlet Control Equations

The inlet behaves as a **weir** at low head and an **orifice** at high head, with a
transition zone between. HW/D is headwater depth above the inlet invert divided by
barrel rise/diameter D.

### Weir (unsubmerged inlet, roughly HW/D < 1.0)

    Q = Cw * L * H^(3/2)

- `Cw` weir coefficient: ~3.0–3.3 (English), ~1.7 (SI)
- `L`  width of opening (= D for circular, span for box)
- `H`  headwater depth measured from the **inlet invert** (channel floor)

### Orifice (submerged inlet, roughly HW/D > 1.2)

    Q = Cd * A * sqrt(2 * g * H)

- `Cd` discharge coefficient: ~0.6 (sharp/square edge) up to ~0.95 (well-rounded)
- `A`  full cross-sectional area of opening (circular: A = pi*D^2/4)
- `g`  32.2 ft/s² (English) or 9.81 m/s² (SI)
- `H`  headwater depth measured from the **centroid** of the opening
       (circular: H = HW_invert - D/2)

> DATUM WARNING: H is measured from the **invert** in the weir equation and from the
> **centroid** in the orifice equation. Do not mix them.

### Transition zone (roughly 1.0 < HW/D < 1.2)

The weir form (Q ∝ H^1.5) and orifice form (Q ∝ H^0.5) generally do **not** agree at the
switch point, so a hard switch produces a discontinuity. HDS-5 instead draws a curve
**tangent to both** across the transition. Implement an interpolation/blend here rather
than an `if` switch (see §6).

### Production accuracy

HDS-5 does not use the bare weir/orifice forms. It uses fitted **regression equations**
with coefficients specific to culvert shape and inlet-edge condition, plus a slope
correction term:
- Unsubmerged: coefficients **K, M** (with an additional slope-correction term).
- Submerged:   coefficients **c, Y**.
For real design, pull these from HDS-5 Appendix A tables (per shape/edge) or call HY-8.

---

## 3. Outlet Control Equation (full-barrel pressure flow)

Valid **only** when the barrel actually flows full under pressure.

    Q = A * sqrt(2 * g * dH) / sqrt(1 + Ke + (29 * n^2 * L) / D^(4/3))

(English units, circular pipe.) Terms:

- `A`   full cross-section, pi*D^2/4
- `dH`  difference between headwater elevation and effective tailwater elevation
- `Ke`  entrance loss coefficient:
          ~0.2 grooved end / headwall,
          ~0.5 square edge,
          ~0.9 projecting
- `n`   Manning roughness: ~0.012 concrete, ~0.024 corrugated metal
- `L`   barrel length
- `D`   diameter
- friction term `29 * n^2 * L / D^(4/3)` is Manning's equation recast as a loss
  coefficient. **SI:** replace the constant 29 with **19.6** (i.e. 2g/Ku² handling);
  verify against HDS-5 for the exact SI grouping before trusting it.

An exit-loss coefficient (~1.0) is sometimes added explicitly inside the denominator;
HDS-5's full energy form lists entrance, friction, and exit losses separately.

### On slope / inclination

Slope does **not** appear as an explicit parameter in the full-barrel formula — it is
absorbed into `dH` via the elevations of the two water surfaces. Slope matters mainly as
a **prerequisite check**: it determines *whether the barrel flows full at all*. A steep
slope drives supercritical free-surface flow so the barrel may never pressurize, in which
case this formula does not apply and inlet control governs. So treat slope as a gate on
applicability, not as a term in the equation.

---

## 4. Free Outfall (Outlet Hanging in Open Air, No Tailwater)

With no downstream water surface, the flow passes through **critical depth** at the
outfall and the outlet exerts no back-pressure. Outlet control almost never governs in
this case — but to run the check anyway, substitute an **effective tailwater**:

    TW = max( yc , (yc + D) / 2 )

where `yc` is critical depth at the outlet (see §5). Use that TW to set `dH` in the
outlet-control formula. The result will nearly always be a lower HW than inlet control,
confirming **inlet control governs**.

Practical workflow:
1. Compute inlet-control HW (weir or orifice per HW/D).
2. Compute critical depth `yc` at the outlet.
3. Compute outlet-control HW using `TW = max(yc, (yc+D)/2)`.
4. Take the higher HW (almost always inlet control for a free outfall).
5. Verify with HDS-5 charts / HY-8 for production.

---

## 5. Critical Depth (yc)

Critical depth is the depth at which **specific energy is minimum** for a given Q, i.e.
where the **Froude number = 1**. It is the threshold between subcritical (deep/slow,
Fr<1) and supercritical (shallow/fast, Fr>1) flow, and it is the depth flow naturally
seeks at a free outfall.

> Conceptual note for anyone reading the curve: the specific-energy diagram (E vs y at
> fixed Q) is a *menu of possible states at one cross-section*, not energy varying along
> the pipe. Real flow along the barrel only loses energy to friction. The "minimum" is
> the least energy *required* to pass Q at all — flow does not descend into it and climb
> back out.

### Rough approximation (circular)

    yc ~= D * 0.325 * (Q / D^2.5)^(2/3) + 0.083 * D

Approximate only — prefer charts or the exact iterative solution below.

### Exact (iterative) for a partly-full circle

Solve for depth `y` where:

    Q^2 * T / (g * A^3) = 1

with circular-segment geometry at depth y, using
`theta = 2 * arccos(1 - 2*y/D)` (radians, the angle subtended by the surface):

    A = (D^2 / 8) * (theta - sin(theta))     # flow area
    T = D * sin(theta / 2)                    # top width

Iterate y in (0, D) (bisection or Newton) until `Q^2 * T / (g * A^3) = 1`. That y is `yc`.
Also useful: hydraulic depth `Dh = A / T`, and `Fr = V / sqrt(g * Dh)` with `V = Q/A`.

### Slope regime check (feeds §3 applicability)

- Mild slope: normal depth > yc -> subcritical -> outlet control possible.
- Steep slope: normal depth < yc -> supercritical -> inlet control governs.

---

## 6. Discontinuities to Expect in a Capacity-vs-Head Curve

When sweeping inlet/outlet head, the implementation can show jumps or kinks. Categorize:

**Artificial (minimize these):**
- *Weir↔orifice hard switch*: value jump unless the §2 tangent transition is implemented.
- *Free-outfall `max()` tailwater*: a kink where one branch of `max(yc, (yc+D)/2)`
  overtakes the other.

**Inherent / physically meaningful (acceptable):**
- *Inlet↔outlet control crossover*: value is continuous (curves meet) but **slope** is
  generally not — a kink at the flow where the bottleneck genuinely shifts. Present in
  HDS-5's own performance curves; honest to keep.
- *Inlet priming / barrel filling*: real culverts show somewhat abrupt behavior as the
  entrance seals or the barrel transitions from part-full to full (turbulence, surging).

**Implementation guidance:**
- Replace the weir/orifice `if` with the HDS-5 transition-zone blend for value continuity.
- Accept the control-crossover kink; do not "smooth" it away — it is real.
- Expect (and document) minor kinks from the tailwater `max()`.
- Within any single regime the curve should be smooth.

---

## 7. Suggested Implementation Outline

```
inputs:  D, L (barrel length), shape, inlet_edge, n, Cw, Cd, Ke,
         HW (or sweep range), TW or free_outfall flag, units
solve_inlet_control(HW):
    compute HW/D
    if HW/D < 1.0:   Q = weir(HW)
    elif HW/D > 1.2: Q = orifice(HW)
    else:            Q = transition_blend(HW)   # tangent to both, not a hard switch
    return Q  (and the implied HW for a target Q)
solve_outlet_control(HW):
    if free_outfall: TW = max(yc(Q), (yc(Q)+D)/2)   # yc via §5 iteration
    dH = HW_elev - TW_elev
    Q  = A * sqrt(2 g dH) / sqrt(1 + Ke + 29 n^2 L / D^(4/3))
    return Q
govern:
    HW_inlet  = inlet HW for target Q
    HW_outlet = outlet HW for target Q
    governing = max(HW_inlet, HW_outlet)   # higher HW / lower capacity wins
checks:
    slope/Froude regime gate on outlet-control applicability (§3, §5)
    verify coefficients against HDS-5 tables for the chosen shape/edge
```

---

## 8. References

- FHWA **HDS-5**, *Hydraulic Design of Highway Culverts*, 3rd ed. — authoritative;
  free PDF from fhwa.dot.gov. Appendix A has the inlet-control regression equations and
  per-shape/edge coefficients; appendices cover hydraulic resistance and design charts.
- FHWA **HY-8** — free culvert analysis software implementing HDS-5 (both control
  conditions, performance curves).
- HDS-5 transition-zone definition: empirical curve drawn tangent to the unsubmerged and
  submerged equations between roughly HW/D = 1.0 and 1.2.
