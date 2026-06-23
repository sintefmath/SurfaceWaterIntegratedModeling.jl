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

> Symbol convention (matches HDS-5 eq. 3.4b / 3.5):
> - `Ku` = the **assembled friction constant** = **29** (English) / **19.63** (SI).
>   This is HDS-5's symbol. It is NOT the Manning unit coefficient.
> - `km` = the **Manning unit coefficient** inside Manning's velocity equation =
>   **1.486** (English) / **1.0** (SI). HDS-5 never names this separately; it is already
>   folded into Ku. Some texts write it as Km.
> - Bridge between them: **Ku = 2g / km^2** (derived below). Do not confuse the two.

HDS-5 builds the required head H as the sum of three velocity-head multiples
(H = Hv + He + Ho, with friction Hf inside; see eq. 3.1–3.5):

    H = [ 1 + Ke + Ku * n^2 * L / R^(1.33) ] * V^2 / (2g)      (HDS-5 eq. 3.5)

with V = Q/A the mean barrel velocity. Solving for Q (since H = dH for full flow):

    Q = A * sqrt(2 * g * dH) / sqrt(1 + Ke + Ku * n^2 * L / R^(1.33))

Terms:

- `A`   full cross-section, pi*D^2/4 for a circular pipe
- `dH`  difference between headwater elevation and effective tailwater elevation
- `Ke`  entrance loss coefficient:
          ~0.2 grooved end / headwall,
          ~0.5 square edge,
          ~0.9 projecting
- `n`   Manning roughness: ~0.012 concrete, ~0.024 corrugated metal
- `L`   barrel length
- `R`   hydraulic radius = A / p (p = wetted perimeter). For a circular pipe flowing
        **full**, R = D/4.
- `Ku`  friction constant = 29 (English) / 19.63 (SI). See bridge formula above.
- exponent `R^(1.33)`: HDS-5 writes 1.33; this is just 4/3 = 1.333... rounded. Use
  `R^(4/3)` in code for full precision — identical quantity.

> CRITICAL: the friction constant Ku pairs with **R^(1.33)** (hydraulic radius), NOT
> D^(1.33). An earlier draft of this file wrote `29 n^2 L / D^(4/3)`, which is wrong —
> mixing the R-based constant with a D-based denominator. Use the R-based form above for
> any cross-section; it is the safe default and matches the manual one-to-one.

### Circular-pipe shortcut (optional)

If you want the term written directly in D for a full circular pipe, substitute R = D/4.
Because (D/4)^(4/3) = D^(4/3) / 4^(4/3) and 4^(4/3) ≈ 6.35, the constant scales up by
~6.35:

    friction term (circular, full) = (Ku * 6.35) * n^2 * L / D^(4/3)
      English: ~184 * n^2 * L / D^(4/3)
      SI:      ~124 * n^2 * L / D^(4/3)

These D-based forms are numerically identical to the R-based form for a full circle —
just pre-substituted. Prefer the R-based form in code to match HDS-5 and avoid confusion.

### Derivation: why Ku = 2g / km^2 = 29 (English) / 19.63 (SI)

Manning's velocity equation uses the unit coefficient km (1.486 English, 1.0 SI):

    V = (km / n) * R^(2/3) * Sf^(1/2)

Solve for the friction slope Sf and multiply by L to get the friction head loss:

    Hf = Sf * L = n^2 * V^2 * L / (km^2 * R^(4/3))

Express Hf as a multiple of the velocity head V^2/(2g) by multiplying and dividing by 2g:

    Hf = [ 2g * n^2 * L / (km^2 * R^(4/3)) ] * V^2/(2g)
       = [ Ku * n^2 * L / R^(4/3) ] * V^2/(2g)     where  Ku = 2g / km^2

Evaluating Ku:
      English: 2 * 32.2 / 1.486^2 = 64.4 / 2.208 = 29.0   ✓ (matches HDS-5)
      SI:      2 * 9.81  / 1.0^2   = 19.62 ≈ 19.63          ✓ (matches HDS-5)

So HDS-5's Ku is exactly this assembled constant; km (the 1.486 / 1.0) lives *inside* it
and is never shown separately in the manual.

> Unit reminder: g = 32.2 ft/s^2 (English) or 9.81 m/s^2 (SI) — same acceleration, two
> unit systems. SI coincidence worth noting: Ku = 19.63 happens to equal the velocity-head
> denominator 2g = 19.62, because km = 1.0 drops out. In English they differ (Ku = 29 vs
> 2g = 64.4) because of the 1.486^2 factor. This coincidence is purely numerical, not a
> shared meaning.


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
    Q  = A * sqrt(2 g dH) / sqrt(1 + Ke + Ku * n^2 * L / R^(4/3))
         # Ku = 29 (English) / 19.63 (SI) = assembled friction constant (HDS-5 eq 3.5)
         # R = A/p (= D/4 for a full circle); R^(4/3) == HDS-5's R^1.33
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
