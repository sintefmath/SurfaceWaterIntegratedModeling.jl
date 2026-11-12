# ============================================================================
# How an NBS placement couples to the terrain.  The link (cells) depends only on the footprint
# and the TrapStructure; the baseline (throughput and its split over the endpoints) is read off
# the NBS-oblivious runoff field.  Both the solver and `network_watercourses` use these, so the
# two cannot disagree.  Invariants: agent/NBS_INVARIANTS.md.
# ============================================================================

# Below this the endpoint weights count as zero and the throughput is split evenly instead.
const _NBS_O0_EPS = 1e-12

"""
    _prepare_nbs!(nbs, tstruct) -> nbs

Resolve the placements against the terrain: set each one's `id` (its position in `nbs`, the
key of its persistent layer storage) and its [`NBSLink`](@ref).  Called by `setup_network`;
idempotent.  This is the only place a placement is mutated.
"""
function _prepare_nbs!(nbs, tstruct)
    for (i, nb) in enumerate(nbs)
        nb.id   = i
        nb.link = _nbs_link(tstruct, nb.footprint)
    end
    return nbs
end

"""
    _nbs_link(tstruct, footprint) -> NBSLink

The cells through which `footprint` exchanges water with the terrain: the outside cells draining
into it, the outside cells it drains to, and the trap bottoms inside it.  Read off the flow
graph's own in/out edges, so the cost is proportional to the footprint.
"""
function _nbs_link(tstruct, footprint::Vector{Int})
    g      = tstruct.flowgraph
    CI     = CartesianIndices(tstruct.topography)
    LI     = LinearIndices(tstruct.topography)
    inside = Set(footprint)
    inflow, exits = Set{Int}(), Set{Int}()
    for c in footprint
        for u in Graphs.inneighbors(g, c)
            u in inside || push!(inflow, u)
        end
        ds, terminus = _downstream_cell(g, c)
        terminus || ds in inside || push!(exits, ds)
    end
    # a closed basin is legal only if its layers drain through piped outlets
    isempty(exits) && @warn "NBS placement footprint has no outflow cells."
    ponds = [b for b in tstruct.trap_bottoms if LI[b] in inside]
    return NBSLink(CI[sort!(collect(inflow))], CI[sort!(collect(exits))], ponds)
end

"""
    _nbs_baseline(nb, tstruct, runoff, precipitation) -> (O_0, exit_shares, pond_shares)

The placement's oblivious throughput `O_0` and how it splits over its endpoints.

`O_0` is the rain on the footprint plus the oblivious flow entering it (exact, as the footprint
has no infiltration).  The shares, aligned with `nb.link.exits` and `nb.link.ponds`, are each
endpoint's part of that throughput: for an exit, the flow the footprint delivers to it (not the
exit cell's own runoff, which also carries water that never crossed the footprint); for a pond,
the flow arriving at it.  They sum to 1, or are equal when `O_0 ≈ 0`.
"""
function _nbs_baseline(nb::DynNBSPlacement, tstruct, runoff, precipitation)
    link = nb.link
    O_0  = _footprint_rain(precipitation, nb.footprint) +
           sum(max(Float64(runoff[c]), 0.0) for c in link.inflow_cells; init = 0.0)
    w_exit = _exit_throughput(tstruct, nb.footprint, link.exits, runoff)
    w_pond = Float64[max(Float64(runoff[c]), 0.0) for c in link.ponds]
    W = sum(w_exit; init = 0.0) + sum(w_pond; init = 0.0)
    @assert isapprox(W, O_0; rtol = 1e-6, atol = 1e-9) || (W < _NBS_O0_EPS && O_0 < _NBS_O0_EPS) """
        NBS endpoint weights ($W) must sum to the throughput O_0 ($O_0)"""
    n = length(w_exit) + length(w_pond)
    share(w) = W > _NBS_O0_EPS ? w / W : 1.0 / n
    return O_0, share.(w_exit), share.(w_pond)
end

# Oblivious flow the footprint hands to each exit (aligned with `exits`): the runoff of the
# footprint cells draining straight out to it.
function _exit_throughput(tstruct, footprint::Vector{Int}, exits, runoff)
    inside = Set(footprint)
    LI = LinearIndices(tstruct.topography)
    w  = Dict{Int,Float64}()
    for c in footprint
        ds, terminus = _downstream_cell(tstruct.flowgraph, c)
        (terminus || ds in inside) && continue
        w[ds] = get(w, ds, 0.0) + max(Float64(runoff[c]), 0.0)
    end
    return Float64[get(w, LI[e], 0.0) for e in exits]
end

"""
    _footprint_rain(precipitation, footprint) -> Float64

Rain landing on `footprint` (scalar or per-cell rate), in the units of the runoff grid.
"""
_footprint_rain(precipitation::Real, footprint::Vector{Int}) =
    Float64(precipitation) * length(footprint)
_footprint_rain(precipitation::AbstractMatrix{<:Real}, footprint::Vector{Int}) =
    sum(Float64(precipitation[f]) for f in footprint; init = 0.0)

"""
    _pond_trap(net, tstruct, cell) -> Int

Local index of the trap in `net` that a pond cell (a trap bottom inside a footprint) drains
into: the single network trap among the supertraps of the cell's region.
"""
function _pond_trap(net::DynNetwork, tstruct, cell::CartesianIndex{2})
    r   = tstruct.regions[cell]
    cov = [l for (l, t) in enumerate(net.traps) if t.trap_ix in tstruct.supertraps_of[r]]
    @assert length(cov) == 1 "NBS pond in region $r maps to $(length(cov)) network traps (expected 1)"
    return cov[1]
end
