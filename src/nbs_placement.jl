# ============================================================================
# Static-analysis support for Nature-Based-Solution (NBS) placements.
#
# These helpers let `spillanalysis` treat each `NBSPlacement` footprint as an
# artificial trap: the terrain under the footprint is dug down to a common low
# level so it forms a single, positive-volume trap, and each system's per-layer
# outlets are validated and resolved against the computed spill regions.  The
# dynamic solver (see `dynamics/`) later drives these traps as rate-limited
# elements; from the static side we only need the footprint to become one
# identifiable trap and the outlets to be resolved.
#
# See `agent/NBS_INTEGRATION_PLAN.md` (§2) and `agent/prompts/nbs_instructions.org`
# for the settled design decisions referenced here as [Qn] / [§n].
# ============================================================================

# Cell size in grid units.  Layer areas and volumes are computed as
# (cell count) * NBS_CELL_AREA.
# @@@ hardcoded to 1.0 (one area-unit per cell), matching the same assumption in
# the `DynCulvert` constructor (elements.jl); read the real grid resolution once
# it is threaded through.
const NBS_CELL_AREA = 1.0

# How far below the global terrain minimum an NBS footprint is dug [Q10].  Any
# positive value guarantees the footprint forms a trap of positive volume (so the
# flat-trap elimination in `spillregions`, which is gated on *zero* volume, never
# removes it).  Kept small so the dug cells stay close to real elevations and are
# easy to recognise when inspecting the surface.
const NBS_DIG_DROP = 1.0

# ----------------------------------------------------------------------------
# 8-connected contiguity test over a footprint (linear indices) on a grid of
# size `dims`.  Used only to warn about disconnected footprints.
function _nbs_footprint_contiguous(footprint::Vector{Int}, dims::Tuple{Int,Int})
    isempty(footprint) && return true
    fset  = Set(footprint)
    CI    = CartesianIndices(dims)
    LI    = LinearIndices(dims)
    seen  = Set{Int}()
    stack = [footprint[1]]
    while !isempty(stack)
        c = pop!(stack)
        (c in seen) && continue
        push!(seen, c)
        r, col = Tuple(CI[c])
        for dr in -1:1, dc in -1:1
            (dr == 0 && dc == 0) && continue
            rr, cc = r + dr, col + dc
            (1 <= rr <= dims[1] && 1 <= cc <= dims[2]) || continue
            n = LI[rr, cc]
            (n in fset && !(n in seen)) && push!(stack, n)
        end
    end
    return length(seen) == length(fset)
end

# ----------------------------------------------------------------------------
"""
    _prepare_nbs!(placements, grid)

Validate and normalise a vector of [`NBSPlacement`](@ref)s against terrain `grid`,
before the spill analysis runs.  Mutates each placement's layer areas in place.

For every placement this checks that the footprint is non-empty, in-bounds and
free of duplicates (warning if it is not 8-connected), that there is exactly one
outlet per layer, and overwrites every layer's area `A` with
`length(footprint) * NBS_CELL_AREA` — the footprint is the single source of truth
for the area [Q6],[Q9].
"""
function _prepare_nbs!(placements, grid::AbstractMatrix)
    dims = size(grid)
    n    = length(grid)
    for (pi, p) in enumerate(placements)
        isempty(p.footprint) && error("NBS placement $pi has an empty footprint")
        for c in p.footprint
            (1 <= c <= n) ||
                error("NBS placement $pi: footprint cell $c is out of bounds (1:$n)")
        end
        (length(unique(p.footprint)) == length(p.footprint)) ||
            error("NBS placement $pi: footprint contains duplicate cells")
        _nbs_footprint_contiguous(p.footprint, dims) ||
            @warn "NBS placement $pi: footprint is not 8-connected; it may form more than one trap"

        nlayers = length(p.system.layers)
        (length(p.outlets) == nlayers) ||
            error("NBS placement $pi: has $(length(p.outlets)) outlets but the " *
                  "system has $nlayers layers (one outlet per layer required)")

        # Footprint is the single source of truth for layer area [Q6],[Q9].
        area = length(p.footprint) * NBS_CELL_AREA
        for layer in p.system.layers
            layer.A = area
        end
    end
    return placements
end

# ----------------------------------------------------------------------------
"""
    _dig_nbs_traps!(grid, placements) -> level

Lower every NBS footprint cell in `grid` to a common fixed elevation so each
footprint forms a single artificial trap [Q10].  The level is
`minimum(grid) - NBS_DIG_DROP`, i.e. strictly below all real terrain, which
guarantees a positive trap volume (a *uniform* flat bottom is fine — the
flat-trap elimination in `spillregions` triggers only on zero volume, not on
flatness).  `grid` is mutated in place; returns the dig level.

Assumes distinct placements have disjoint, non-adjacent footprints; adjacent
footprints dug to the same level would merge into one trap.
"""
function _dig_nbs_traps!(grid::AbstractMatrix, placements)
    isempty(placements) && return nothing
    # Note: `grid` here is the modelling copy (`gridcpy`) in `spillanalysis`, so
    # the caller's original terrain is untouched.  Height unit assumed metres if
    # ambiguous [Q3]; NBS_DIG_DROP is in that same unit.
    level = minimum(grid) - NBS_DIG_DROP
    for p in placements
        for c in p.footprint
            grid[c] = level
        end
    end
    return level
end

# ----------------------------------------------------------------------------
"""
    _resolve_nbs!(placements, regions, spillpoints, topography)

After the spill analysis has run on the dug terrain, resolve each placement's
per-layer outlets against the computed `regions`/`spillpoints`, mutating the
`outlets` vectors in place.

For each placement:
- The NBS trap index is looked up from the footprint's region [Q5]; every
  footprint cell must share one region (otherwise the footprint failed to form a
  single trap and an error is raised).
- Each *specified* outlet must lie outside the NBS trap's own spill region [§25]
  and *strictly below* the trap's natural spillpoint elevation [§15] — the
  latter guaranteeing water cannot route back into the trap (the cycle guard,
  see plan §4).
- An *unspecified* outlet `CartesianIndex(0,0)` [Q7] is backfilled: the lowermost
  layer's takes the trap's natural discharge cell (the *downstream* side of the
  trap's spillpoint — where the trap overflows to, not the trap-side cell); every
  other layer's takes the (resolved) lowermost outlet.  Each backfill emits a
  warning [Q7a],[Q4].
- Finally, no outlet may share a supertrap with the NBS trap: if the outlet's
  region and the NBS trap have a common ancestor in the sub/supertrap hierarchy,
  a later fill of that supertrap would merge the two and route water back into
  the NBS — a cycle.  This is rejected (see [`_check_nbs_no_shared_supertrap`]).

Layers are ordered top→bottom (`system.layers`), so the lowermost layer — the one
whose outlet defines the trap's discharge point — is the last entry, and
`outlets[end]` is its outlet.

Note: this does *not* rewrite the computed `Spillpoint` of the NBS trap.  Per
[Q3b] the artificial trap's static spillpoint and volume are superseded by the
dynamic solver, which routes discharge through the resolved `outlets`; forcing
the static spillpoint would only desync the already-computed trap volumes and
hierarchy.  @@@ If a future use-case needs the static spillgraph to reflect the
specified bottom outlet, revisit here.
"""
function _resolve_nbs!(placements, regions::AbstractMatrix{Int},
                       spillpoints, topography::AbstractMatrix, supertraps_of)
    dims = size(regions)
    CI   = CartesianIndices(dims)
    unspecified = CartesianIndex(0, 0)
    for (pi, p) in enumerate(placements)
        trap_ix = regions[p.footprint[1]]
        trap_ix > 0 ||
            error("NBS placement $pi: footprint did not form a valid trap " *
                  "(region $trap_ix); check the footprint and dig level")
        all(regions[c] == trap_ix for c in p.footprint) ||
            error("NBS placement $pi: footprint spans more than one region " *
                  "(expected a single trap $trap_ix); footprints must be contiguous")

        sp      = spillpoints[trap_ix]
        sp_elev = sp.elevation
        # Natural discharge cell: the downstream side of the trap's spillpoint,
        # i.e. the cell the trap overflows *into* (not `current_region_cell`,
        # which lies inside the NBS region).  `nothing` when the trap spills out
        # of the domain (no downstream cell exists).
        nat_cell = sp.downstream_region_cell == 0 ? nothing : CI[sp.downstream_region_cell]
        nlayers  = length(p.outlets)

        # Resolve the lowermost outlet first so upper layers can default to it.
        if p.outlets[nlayers] == unspecified
            nat_cell === nothing &&
                error("NBS placement $pi: lowermost outlet unspecified and the " *
                      "trap spills out of the domain; specify an explicit outlet")
            p.outlets[nlayers] = nat_cell
            @warn "NBS placement $pi: lowermost outlet unspecified; set to the " *
                  "trap's natural discharge cell at $(Tuple(nat_cell))"
        else
            _validate_nbs_outlet(pi, "lowermost", p.outlets[nlayers], trap_ix,
                                  regions, topography, sp_elev)
            if nat_cell !== nothing && p.outlets[nlayers] != nat_cell
                @warn "NBS placement $pi: specified lowermost outlet " *
                      "$(Tuple(p.outlets[nlayers])) differs from the trap's " *
                      "natural discharge cell $(Tuple(nat_cell))"
            end
        end
        bottom = p.outlets[nlayers]

        for l in 1:(nlayers - 1)
            if p.outlets[l] == unspecified
                p.outlets[l] = bottom
                @warn "NBS placement $pi: outlet for layer $l unspecified; set to " *
                      "the lowermost outlet at $(Tuple(bottom))"
            else
                _validate_nbs_outlet(pi, "layer $l", p.outlets[l], trap_ix,
                                     regions, topography, sp_elev)
            end
        end

        # Cycle guard: no outlet may share a supertrap with the NBS trap.
        for (l, outlet) in enumerate(p.outlets)
            _check_nbs_no_shared_supertrap(pi, l, outlet, trap_ix, regions, supertraps_of)
        end
    end
    return placements
end

# Reject an outlet that shares a supertrap with the NBS trap: if the outlet's
# region and the NBS trap have a common ancestor in the sub/supertrap hierarchy,
# filling that supertrap would submerge both and route water back into the NBS.
# An out-of-domain / unclassified outlet region (<= 0) shares nothing and passes.
function _check_nbs_no_shared_supertrap(pi, layer, outlet::CartesianIndex{2},
                                        trap_ix::Int, regions::AbstractMatrix{Int},
                                        supertraps_of)
    checkbounds(Bool, regions, outlet) || return nothing
    r = regions[outlet]
    (r > 0 && r <= length(supertraps_of)) || return nothing
    shared = intersect(supertraps_of[trap_ix], supertraps_of[r])
    isempty(shared) ||
        error("NBS placement $pi: outlet for layer $layer (region $r) shares " *
              "supertrap(s) $(collect(shared)) with the NBS trap (region $trap_ix); " *
              "water would route back once that supertrap fills [cycle guard]")
    return nothing
end

# Validate a single specified outlet cell against the NBS trap.
function _validate_nbs_outlet(pi, label, outlet::CartesianIndex{2}, trap_ix::Int,
                              regions::AbstractMatrix{Int}, topography::AbstractMatrix,
                              sp_elev)
    (checkbounds(Bool, regions, outlet)) ||
        error("NBS placement $pi: $label outlet $(Tuple(outlet)) is out of bounds")
    regions[outlet] != trap_ix ||
        error("NBS placement $pi: $label outlet $(Tuple(outlet)) lies inside the " *
              "NBS trap's own spill region — it must discharge elsewhere [§25]")
    topography[outlet] < sp_elev ||
        error("NBS placement $pi: $label outlet $(Tuple(outlet)) is not strictly " *
              "below the trap's natural spillpoint elevation ($sp_elev); water " *
              "could route back into the trap [§15]")
    return nothing
end
