using Printf

export NBSLayer, NBSSystem
export mantillaRRmodel, elhadiGreenRoof, puddle
export compute_outflow

# ----------------------------------------------------------------------------
"""
    NBSLayer

A single storage layer in an NBS (Nature-Based Solution) system.

Storage triggers outflow (`Kout`, `nout`) when `S > Smax`, and infiltration to
the next layer (`Kinf`, `ninf`) when `S > Smin`. Evapotranspiration scales
linearly with storage up to the threshold `EVS11`, and at full rate above it.

# Fields
- `Smin`: lower storage threshold — infiltration activates above this (mm)
- `Smax`: upper storage threshold — outflow activates above this (mm)
- `Kout`, `nout`: outflow rate as `Kout * max(S - Smax, ε)^nout`
- `Kinf`, `ninf`: infiltration rate as `Kinf * max(S - Smin, ε)^ninf`
- `EVCoeff`: evapotranspiration coefficient (scalar or `Function(t)`)
- `EVS11`: storage level at which ET reaches full rate
- `A`: horizontal area of the layer (m²)
- `layer_name`: display name used in printing
"""
mutable struct NBSLayer
    Smin    # minimum capacity limit
    Smax    # maximum capacity limit
    Kout    # outflow rate coefficient
    Kinf    # infiltration rate coefficient
    nout    # outflow exponent
    ninf    # infiltration exponent
    EVCoeff # evapotranspiration coefficient
    EVS11   # soil moisture threshold
    A       # area of the layer
    layer_name::String
end

# ----------------------------------------------------------------------------
function Base.show(io::IO, layer::NBSLayer)
    println(io, "\nNBSLayer: '" * layer.layer_name * "'\n")

    println(io, "  - Storage params: " *
        "Smin=$(Printf.@sprintf("%.4g", layer.Smin)), " *
        "Smax=$(Printf.@sprintf("%.4g", layer.Smax))")
    if layer.Kout > 0.0
        println(io, "  - Outflow params: " *
                "K=$(Printf.@sprintf("%.4g", layer.Kout)), " *
                "n=$(Printf.@sprintf("%.4g", layer.nout))")
    else
        println(io, "  - There is no outflow from this layer.")
    end
    if layer.Kinf > 0.0
        println(io, "  - Infiltration params: " *
            "K=$(Printf.@sprintf("%.4g", layer.Kinf)), " *
            "n=$(Printf.@sprintf("%.4g", layer.ninf))")
    else
        println(io, "  - There is no infiltration from this layer.")
    end
    if layer.EVCoeff isa Function
        println(io, "  - Evapotranspiration params: " *
            "EVCoeff=function, " *
            "EVS11=$(Printf.@sprintf("%.4g", layer.EVS11))")
    elseif layer.EVCoeff > 0.0
        println(io, "  - Evapotranspiration params: " *
            "EVCoeff=$(Printf.@sprintf("%.4g", layer.EVCoeff)), " *
            "EVS11=$(Printf.@sprintf("%.4g", layer.EVS11))")
    else
        println(io, "  - There is no evapotranspiration from this layer.")
    end

    println(io, "  - Area: " * "$(Printf.@sprintf("%.4g", layer.A))")
end

# ----------------------------------------------------------------------------
"""
    NBSSystem

An NBS model consisting of vertically stacked [`NBSLayer`](@ref)s. Infiltration
from layer `i` feeds layer `i+1` automatically; no explicit connection is needed
within the same system.

# Fields
- `layers`: ordered vector of layers from top (surface) to bottom
- `name`: display name used in printing
"""
mutable struct NBSSystem
    layers::Vector{NBSLayer} # stack of NBS layers, from upper to lower
    name::String # String with model name
end

# ----------------------------------------------------------------------------
function Base.show(io::IO, model::NBSSystem)
    N = length(model.layers)
    println(io, "\nNBSSystem with name '" * model.name * "' and $N layers:\n")
    for i in 1:N
        layer = model.layers[i]
        println(io, " Layer $i:")
        show(io, layer)
        println(io, "")
    end
end

# ----------------------------------------------------------------------------
"""
    mantillaRRmodel(Smax, SM, Sdrain, kQ, kINF, kPER, nQ, nINF, nPER, EVC, EVS11;
                    area=1.0, name="Mantilla RR Model") -> NBSSystem

Three-layer rainfall-runoff model from Mantilla PhD thesis (2025), §3.7.3 / Figure 13.

Layers (top to bottom): **Surface** → **Soil** → **Drainage**. Surface overflows
to the environment; soil percolates to drainage; drainage discharges to the outlet.
Evapotranspiration acts on both surface and soil layers.

# Arguments
- `Smax`: surface layer storage capacity (mm)
- `SM`: soil layer holding capacity — water is held until `S > SM` (mm)
- `Sdrain`: drainage layer capacity (mm)
- `kQ`: outflow coefficient for surface and drainage layers
- `kINF`: infiltration coefficient from surface to soil
- `kPER`: percolation coefficient from soil to drainage
- `nQ`, `nINF`, `nPER`: power-law exponents for the respective fluxes
- `EVC`: evapotranspiration coefficient (applied to surface and soil)
- `EVS11`: storage threshold at which ET reaches full rate
"""
function mantillaRRmodel(Smax, SM, Sdrain,
                         kQ, kINF, kPER,
                         nQ, nINF, nPER,
                         EVC, EVS11; area = 1.0, name="Mantilla RR Model")
    toplayer = NBSLayer(
        0.0,            # Smin
        Smax,          # Smax
        kQ,            # Kout (overflow)
        kINF,          # Kinf (infiltration)
        nQ,            # nout
        nINF,          # ninf
        EVC,           # EVCoeff
        EVS11,         # EVS11
        area,          # layer area
        "Surface storage"
    )
    midlayer = NBSLayer(
        SM,            # Smin
        0.0,          # Smax
        0.0,          # Kout (overflow)
        kPER,         # Kinf (infiltration)
        0.0,          # nout
        nPER,         # ninf
        EVC,          # EVCoeff
        EVS11,        # EVS11
        area,         # layer area
        "Soil storage"
    )
    bottomlayer = NBSLayer(
        0.0,          # Smin
        Sdrain,       # Smax
        kQ,           # Kout (overflow)
        0.0,          # Kinf (infiltration)
        nQ,           # nout
        0.0,          # ninf
        0.0,          # EVCoeff
        1.0,          # EVS11
        area,         # layer area
        "Drainage storage"
    )
    NBSSystem([toplayer, midlayer, bottomlayer], name)
end

# ----------------------------------------------------------------------------
"""
    elhadiGreenRoof(Smin_soil, Smin_drain, kINF, kDRAIN, nINF, nDRAIN, EVC, EVS11;
                    area=1.0, name="Elhadi Green Roof Model") -> NBSSystem

Two-layer green roof model: **Soil** → **Drainage**.

Water is retained in the soil layer until `S > Smin_soil`, then infiltrates to
drainage at rate `kINF * (S - Smin_soil)^nINF`. Drainage discharges to the
outlet at rate `kDRAIN * (S - Smin_drain)^nDRAIN`.

# Arguments
- `Smin_soil`: soil field capacity — water held before drainage begins (mm)
- `Smin_drain`: drainage layer capacity (mm)
- `kINF`, `nINF`: infiltration coefficient and exponent (soil → drainage)
- `kDRAIN`, `nDRAIN`: outflow coefficient and exponent (drainage → outlet)
- `EVC`: evapotranspiration coefficient
- `EVS11`: storage threshold at which ET reaches full rate
"""
function elhadiGreenRoof(Smin_soil, Smin_drain,
                         kINF, kDRAIN,
                         nINF, nDRAIN,
                         EVC, EVS11;
                         area = 1.0,
                         name="Elhadi Green Roof Model")
    toplayer = NBSLayer(
        Smin_soil,     # Smin
        Inf,           # Smax
        0.0,           # Kout (overflow)
        kINF,          # Kinf
        0.0,           # nout
        nINF,          # ninf
        EVC,           # EVCoeff
        EVS11,         # EVS11
        area,          # layer area
        "Soil layer"
    )
    bottomlayer = NBSLayer(
        0.0,           # Smin
        Smin_drain,    # Smax
        kDRAIN,        # Kout (overflow)
        0.0,           # Kinf
        nDRAIN,        # nout
        0.0,           # ninf
        EVC,           # EVCoeff
        1.0,           # EVS11
        area,          # layer area
        "Drainage layer"
    )
    NBSSystem([toplayer, bottomlayer], name)
end

# ----------------------------------------------------------------------------
"""
    puddle(Smax; kOUT=1.0, nOUT=1.0, kINF=0.0, nINF=0.0,
           EVC=0.0, EVS11=1.0, area=1.0, name="Simple Puddle Model") -> NBSSystem

Single-layer surface storage. Overflows at rate `kOUT * max(S - Smax, ε)^nOUT`
once storage exceeds `Smax`. Defaults produce a pure overflow model with no
infiltration or evapotranspiration.
"""
function puddle(Smax;
                kOUT=1.0, nOUT=1.0, kINF=0.0, nINF=0.0, EVC=0.0, EVS11=1.0,
                area = 1.0,
                name="Simple Puddle Model")
    NBSSystem(
       [NBSLayer(
           0.0,           # Smin
           Smax,          # Smax
           kOUT,          # Kout
           kINF,          # Kinf
           nOUT,          # nout
           nINF,          # ninf
           EVC,           # EVCoeff
           EVS11,         # EVS11
           area,          # layer area
           "Puddle layer"
       )],
        name
    )
end

# ----------------------------------------------------------------------------
"""
    compute_outflow(K, n, Smax, S; mollifier=0.0) -> Float64

Power-law outflow: `K * max(S - Smax, ε)^n`.

`mollifier > 0` smooths the threshold over a fraction of `Smax`, which improves
gradient conditioning when this function is differentiated (e.g. via ForwardDiff,
as done for calibration in the standalone NBS package).
"""
function compute_outflow(K, n, Smax, S; mollifier=0.0)
    tmp = max(S - Smax, eps())

    if mollifier > 0
        mu = mollifier * Smax
        tmp += mu * (  (S < Smax) ? S/Smax : max(2*Smax - S, 0.0)/Smax )
    end
    return K * tmp^n
end

# This file holds only the layer-storage model, so it can be included before elements.jl;
# the placement type (`DynNBSPlacement`) lives with the other Dyn* elements in elements.jl.
