import Graphs

export DynObject, DynFlowPath, DynTrap, DynCulvert, DynNBSPlacement, DynNetwork, setup_network

# Make generic baseclass for dynamic objects
abstract type DynObject end

"""
        DynFlowPath(cells, departure_point, target_trap, culvert_inlets, culvert_outlets, nbs_outlets, merges)

Represent a flow path over the terrain as a sequence of grid cells.  The path may
lead into a trap (`target_trap > 0`), terminate in an intersecting flow path, or
flow out of the domain (`target_trap == 0` for both).

Culverts and NBS outlets along the way subtract or add water; the infiltration
capacity of each cell is represented externally.  Convenience constructors default
`departure_point` to the head cell and the culvert/nbs/merge lists to empty.

"""
struct DynFlowPath <: DynObject
    cells::Vector{CartesianIndex{2}} # cells along the flow path

    # Cell the path departed from; valid even when `cells` is empty (zero-length connector)
    # or the tail was truncated. Authoritative source — do not infer it from `cells[1]`.
    departure_point::CartesianIndex{2}

    # Target trap index (0 for out-of-domain or intersection with another flow path)
    target_trap::Int

    # culvert inlets and outlets on this path: each (culvert_index, cell_position),
    # where cell_position is the 1-based index of the culvert's inlet/outlet cell in
    # *this* path's `cells`.  The position lets routing charge the infiltration up to
    # the abstraction/addition point, exactly like a `merges` junction.
    culvert_inlets::Vector{Tuple{Int,Int}}
    culvert_outlets::Vector{Tuple{Int,Int}}

    # nbs outlets represent other external sources
    nbs_outlets::Vector{Tuple{Int,Int}}

    # tributary paths that merge into this one: (tributary_path_index, junction_cell_index)
    # where junction_cell_index is the 1-based index of the junction cell in *this* path's cells.
    merges::Vector{Tuple{Int,Int}}
end

# content-first form: departure_point defaults to the head cell (requires non-empty cells;
# a zero-length connector uses the full constructor with an explicit departure_point).
DynFlowPath(cells, target_trap, cin, cout, nbs, mg) =
    DynFlowPath(cells, first(cells), target_trap, cin, cout, nbs, mg)
DynFlowPath(cells, target_trap) =
    DynFlowPath(cells, target_trap, Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[])
DynFlowPath(cells) =
    DynFlowPath(cells, 0, Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[])

"""
        DynTrap(trap_ix, spill_path, culvert_inlets, culvert_outlets)

Represent a trap in the terrain, identified by its index in the spillanalysis
structure.  The `spill_path` field is the flow path that water would take when the
trap's capacity is exceeded:
- `spill_path > 0`: the (network-local) index of that flow path;
- `spill_path == 0`: the trap is not (yet) full, so it has no spill path;
- `spill_path == -1`: the trap is full and spills straight out of the domain
  (no in-network successor).  The routing and ordering layers treat any
  `spill_path <= 0` as "no successor"; the sentinel only distinguishes a full
  domain-exiting trap from a not-yet-full one for the solver's state validation.

The trap may also have culvert inlets or outlets within its footprint, which
would add or subtract water from the trap respectively, depending on the current
water level in the trap.  The infiltration capacity of each cell in the trap is
represented externally.

"""
mutable struct DynTrap <: DynObject
    trap_ix::Int # index of the trap in the spillanalysis structure

    # spill path
    spill_path::Int

    # culvert inlets and outlets, each represented by a culvert ID
    culvert_inlets::Vector{Int}
    culvert_outlets::Vector{Int}

    # Live incoming flow connections to dynamic sources (reachability count); the trap is
    # dynamic while in_count > 0.  Set by init_in_counts! and maintained on grow/shrink.
    in_count::Int

    DynTrap(trap_ix, spill_path, culvert_inlets, culvert_outlets) =
        new(trap_ix, spill_path, culvert_inlets, culvert_outlets, 0)
end

DynTrap(trap_ix, spill_path) = DynTrap(trap_ix, spill_path, Int[], Int[])


"""
       DynCulvert(inlet, outlet, r, Cd, Ke, Kf, Cw)

Represent a culvert connecting two points in the terrain, represented by the
grid cell indices of the inlet and outlet.  The culvert has a limited capacity,
which is determined by its internal parameters as well as the dynamic water levels
at the inlet and outlet.

Internal parameters include the radius (r), discharge coefficient (Cd), entrance
loss coefficient (Ke), friction loss coefficient (Kf), and weir coefficient
(Cw).  Which parameters are used in the computation of the flow through the
culvert depends on the water levels at the inlet and outlet.

"""
struct DynCulvert <: DynObject
    inlet::CartesianIndex{2}  # grid cell index of the culvert inlet
    outlet::CartesianIndex{2} # grid cell index of the culvert outlet

    r::Float64  # radius of the culvert
    Cd::Float64 # discharge coefficient of the culvert
    Ke::Float64 # entrance loss coefficient
    Kf::Float64 # friction loss coefficient
    Cw::Float64 # weir coefficient (for overtopping flow)
end

"""
    DynCulvert(tstruct, inlet, outlet; r, n=0.013, Cd=0.6, Ke=0.5, Cw=1.7)

Convenience constructor that builds a [`DynCulvert`](@ref) from a minimum of
physical data, filling in SI default hydraulic coefficients and deriving the
friction-loss coefficient `Kf` from Manning's roughness `n` and the barrel
length.

The barrel length is the straight-line distance between the `inlet` and `outlet`
cells: horizontal extent plus the elevation drop read from `tstruct.topography`.
All quantities are SI (metres, m^3/s).

Defaults assume a concrete pipe (`n`) with a square-edged entrance (`Cd`, `Ke`)
and the SI weir coefficient (`Cw`); override any of them when better data exists.
"""
function DynCulvert(tstruct, inlet::CartesianIndex{2}, outlet::CartesianIndex{2};
                    r::Real,
                    n::Real  = 0.013,   # Manning roughness (~concrete)
                    Cd::Real = 0.6,     # orifice discharge coef (square edge)
                    Ke::Real = 0.5,     # entrance loss coef (square edge)
                    Cw::Real = 1.7)     # weir coefficient, SI
    D = 2r
    di, dj = Tuple(outlet - inlet)
    # @@@ grid resolution assumed 1 m/cell; use the real cell size once available
    horiz = hypot(float(di), float(dj))
    drop  = tstruct.topography[inlet] - tstruct.topography[outlet]
    L = hypot(horiz, drop)              # full barrel length (m)
    # Manning friction as a dimensionless loss coefficient (HDS-5 eq. 3.5, SI):
    # Kf = Ku * n^2 * L / R^(4/3), with assembled friction constant Ku = 19.63
    # (SI; = 2g/km^2 with Manning unit coef km = 1) and hydraulic radius R = D/4
    # for a full circular pipe (R-based form, not D-based).
    # @@@ R = D/4, not D/2: hydraulic radius is A/P = pi*r^2 / (2*pi*r) = r/2 = D/4.
    R  = D / 4                          # hydraulic radius, full circular pipe
    Kf = 19.63 * n^2 * L / R^(4/3)
    return DynCulvert(inlet, outlet, float(r), float(Cd), float(Ke), float(Kf), float(Cw))
end

# ----------------------------------------------------------------------------
"""
       DynNBSPlacement(system, footprint, n_terrain, outlets)
       DynNBSPlacement(system, footprint, outlets)          # n_terrain defaults to 1

Represent a Nature-Based Solution (NBS) installation placed on the terrain,
governed by the layered storage model `system` (an [`NBSSystem`](@ref)).

`footprint` is the set of grid cells (linear indices, matching
`TrapStructure.footprints`'s convention) covered by the installation.

`n_terrain` is the number of *topmost* layers whose overflow re-emits at the
footprint's lower-edge exit boundary.  Every outflowing layer (`Kout > 0`) *below*
the top `n_terrain` requires one explicit piped `outlet`, in top-to-bottom layer
order.  `outlets` are grid cells (geometry — matching `DynCulvert`'s
`inlet`/`outlet` convention) where that piped discharge re-enters the terrain; the
terrain-re-emit targets of the top `n_terrain` layers are derived from the
footprint geometry at network-build time.

The `footprint_*`/`internal_accumulation_cells` and `id` fields are filled in at
network-build time (`setup_network`), when the `TrapStructure` is available.
"""
mutable struct DynNBSPlacement <: DynObject
    system::NBSSystem
    footprint::Vector{Int}             # grid cells covered (linear indices)
    n_terrain::Int                     # topmost layers re-emitting at the exit boundary
    outlets::Vector{CartesianIndex{2}} # piped-outlet cells for the outflowing layers below the top n

    # filled in at network-build time, when the TrapStructure is available
    footprint_inflow_cells::Vector{CartesianIndex{2}}      # cells feeding into the footprint from outside
    footprint_outflow_cells::Vector{CartesianIndex{2}}     # external cells the footprint drains to
    internal_accumulation_cells::Vector{CartesianIndex{2}} # footprint cells that pond (trap bottoms)

    # @@@ stable identity for cross-event layer-state persistence (the key for restoring
    #     each placement's layer state across rebuilds).  Assigned by the caller/distributor;
    #     0 until then.
    id::Int

    function DynNBSPlacement(system::NBSSystem, footprint::Vector{Int}, n_terrain::Integer,
                             outlets::Vector{CartesianIndex{2}})
        nlayer = length(system.layers)
        0 <= n_terrain <= nlayer ||
            error("DynNBSPlacement: n_terrain must be in 0:$nlayer (layer count), got $n_terrain")
        # every outflowing layer (Kout > 0) below the top n_terrain needs one piped outlet
        npiped = count(l -> system.layers[l].Kout > 0.0, (n_terrain + 1):nlayer)
        length(outlets) == npiped ||
            error("DynNBSPlacement: expected $npiped piped outlet(s) for the outflowing layer(s) " *
                  "below the top n_terrain=$n_terrain, got $(length(outlets))")
        new(system, footprint, n_terrain, outlets, Vector{CartesianIndex{2}}(),
            Vector{CartesianIndex{2}}(), Vector{CartesianIndex{2}}(), 0)
    end
end

# Backward-compatible tail: no explicit n_terrain defaults to 1 (topmost layer re-emits at terrain).
DynNBSPlacement(system::NBSSystem, footprint::Vector{Int}, outlets::Vector{CartesianIndex{2}}) =
    DynNBSPlacement(system, footprint, 1, outlets)

# Network of dynamic objects
"""
         DynNetwork(flow_paths, traps, culverts[, nbs])

Represent the dynamic elements of the terrain as a network of flow paths, traps,
and culverts.  Each flow path may lead into a trap, and each trap has a spill
path that represents the flow path that water would take when the trap's
capacity is exceeded.  Culverts can be present along flow paths or within traps,
and would add or subtract water from the flow path or trap depending on their
location and the current water levels.

The network is tree-like (flow paths can merge but not split, and multiple flow
paths may flow into the same trap, but each flow path can only lead into one
trap and each trap only has one spill path out of it when it overflows).  Loops
are avoided by requiring that water always flows downstreams (whether through
flow paths or culverts).  Overlapping flow paths are not allowed; if multiple
flow paths share cells, one would be truncated as it merges into the other, and
the trap at the end of the truncated path would be connected to the remaining
path instead.
"""
mutable struct DynNetwork
    flow_paths::Vector{DynFlowPath}
    traps::Vector{DynTrap}
    culverts::Vector{DynCulvert}
    nbs::Vector{DynNBSPlacement}
end

# Backward-compatible constructor: a network with no NBS overlay elements.  Keeps
# every `DynNetwork(paths, traps, culverts)` construction site working unchanged
# while the NBS overlay wiring is layered on (populated only where NBS are present).
DynNetwork(flow_paths::Vector{DynFlowPath}, traps::Vector{DynTrap},
           culverts::Vector{DynCulvert}) =
    DynNetwork(flow_paths, traps, culverts, DynNBSPlacement[])

DynNetwork(culverts::Vector{DynCulvert}, nbs::Vector{DynNBSPlacement}) =
    DynNetwork(DynFlowPath[], DynTrap[], culverts, nbs)

DynNetwork() = DynNetwork(DynFlowPath[], DynTrap[], DynCulvert[])
