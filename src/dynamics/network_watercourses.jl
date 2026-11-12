export network_watercourses

# ----------------------------------------------------------------------------
"""
    network_watercourses(tstruct, full_traps, trap_volumes, nbs_layers;
                         dyn_traps, culverts, nbs, precipitation, infiltration)

Compute the flow-intensity field across the terrain at ONE timepoint, taking the dynamic
network (culverts and NBS) into account.

[`watercourses`](@ref) computes the flow field routed purely by terrain — it is deliberately
network-oblivious, and is the baseline the dynamic solver corrects.  This function produces the
network-aware field for a single instant: it starts from that oblivious baseline and rewrites it
where the network diverts flow — a culvert draws water at its inlet cell and delivers it at its
outlet, and an NBS placement intercepts the inflow reaching its footprint and re-emits its
layers' current outflow at its exits.

Surface flow is undefined inside an NBS footprint (the water is held by the layered-storage
model, not routed across the terrain), so footprint cells carry **no** surface flow: their field
value is set to `NaN`, marking them as undefined rather than zero.  `NaN` also masks cleanly in
log-scale flow-intensity plots, where a zero would be `-Inf`.  The intercepted flow reappears
only at the placement's exit cells, from where it propagates downstream as usual.

It is a *snapshot*: the culvert and NBS flows depend on the network's state at the timepoint,
which the caller supplies.  That state is exactly what the timepoint-query functions return —
`full_traps` / `trap_volumes` from [`all_states_at_timepoints`](@ref) and `nbs_layers` from
[`network_states_at_timepoints`](@ref) — evaluated for the same `dyn_traps` / `culverts` / `nbs`.

# Arguments
- `tstruct::TrapStructure{<:Real}`: the terrain trapping structure from `spillanalysis`.
- `full_traps::Vector{Bool}`: which traps are filled at the timepoint (one entry per trap).
- `trap_volumes::Dict{Int,Float64}`: network-covered trap index -> water volume at the timepoint,
  used to read the inlet/outlet water levels a culvert's rate depends on.
- `nbs_layers::Dict{Int,Vector{Float64}}`: NBS placement `id` -> per-layer water volumes at the
  timepoint (upper layer first), driving each placement's re-emission.
- `dyn_traps::Vector{Int}`: dynamic-trap seeds, as passed to [`fill_sequence`](@ref) — together
  with the culvert/NBS seeds these define the network membership (default: empty).
- `culverts::Vector{DynCulvert}`: the culverts, as passed to [`fill_sequence`](@ref)
  (default: empty).
- `nbs::Vector{DynNBSPlacement}`: the NBS placements, as passed to [`fill_sequence`](@ref)
  (default: empty).
- `precipitation::Union{Matrix{<:Real}, Real}`: precipitation rate per grid cell (default: 1.0).
- `infiltration::Union{Matrix{<:Real}, Real}`: maximum infiltration rate per grid cell
  (default: 0.0).
- `zvt`: optional precomputed z-volume tables (`_compute_z_vol_tables(tstruct)`). This depends
  only on `tstruct`; pass it to avoid rebuilding it on every call when sampling many timepoints
  (default: `nothing`, built on demand).

# Returns
- `runoff::Matrix{Float64}`: the network-aware flow-intensity field — infiltration-excess runoff
  rate (positive) or remaining infiltration capacity (negative) per grid cell, with the culvert
  and NBS diversions applied.

See also [`watercourses`](@ref), [`network_states_at_timepoints`](@ref),
[`all_states_at_timepoints`](@ref).
"""
function network_watercourses(tstruct::TrapStructure{<:Real},
                              full_traps::Vector{Bool},
                              trap_volumes::Dict{Int,Float64},
                              nbs_layers::Dict{Int,Vector{Float64}};
                              dyn_traps::Vector{Int}=Int[],
                              culverts::Vector{DynCulvert}=DynCulvert[],
                              nbs::Vector{DynNBSPlacement}=DynNBSPlacement[],
                              precipitation::Union{Matrix{<:Real}, Real}=1.0,
                              infiltration::Union{Matrix{<:Real}, Real}=0.0,
                              zvt=nothing)

    gridres = size(tstruct.topography)
    LI      = LinearIndices(gridres)

    # expand precipitation/infiltration to grids; NBS footprints never infiltrate (the layer
    # models own that water) — the same rule fill_sequence applies
    precip = precipitation isa Real ? fill(Float64(precipitation), gridres) : Float64.(precipitation)
    infil  = infiltration  isa Real ? fill(Float64(infiltration),  gridres) : Float64.(infiltration)
    isempty(nbs) || (infil[vcat([n.footprint for n in nbs]...)] .= 0.0)

    # the network-oblivious baseline field, fully routed by terrain. `runoff_live`
    # accumulates the network corrections; `runoff_orig` keeps the untouched
    # oblivious field. Culvert draws read the *live* flow (what is actually there,
    # matching the solver's min(capacity, current)); NBS corrections are measured
    # against the *oblivious* field (matching the solver's precomputed O_0_total),
    # so a downstream placement is never contaminated by an upstream correction.
    runoff_live, = watercourses(tstruct, full_traps; precipitation = precip, infiltration = infil)
    runoff_orig  = copy(runoff_live)

    # build the network components: topology, cached NBS cell lists, culvert owners.
    # `z_vol_tables` depends only on `tstruct`, so a caller sampling many timepoints
    # should compute it once and pass it via `zvt`; otherwise it is built here.
    z_vol_tables = zvt !== nothing ? zvt : _compute_z_vol_tables(tstruct)
    seeds = _dyn_seeds(tstruct, dyn_traps, DynCulvert[])
    comps = setup_network(tstruct, findall(full_traps);
                          dyn_coords = seeds, culverts = culverts, nbs = nbs)

    footprint_set  = Set{Int}(vcat([n.footprint for n in nbs]...))
    full_traps_int = findall(full_traps)

    for net in comps
        # --- culverts: draw at the inlet, deliver at the outlet (mass-conserving) ------------
        if !isempty(net.culverts)
            cvplan = _build_culvert_plan(net, tstruct)
            geom   = _build_trap_geometry(tstruct, net, infil; zvt = z_vol_tables)
            trap_level = Float64[water_level(geom[i], get(trap_volumes, net.traps[i].trap_ix, 0.0))
                                 for i in eachindex(net.traps)]
            for ci in eachindex(net.culverts)
                inlet  = LI[net.culverts[ci].inlet]
                outlet = LI[net.culverts[ci].outlet]
                # capacity capped by the LIVE flow actually reaching the inlet — never draw more
                # than is available (mass conservation), including any upstream network corrections
                q = min(_culvert_flow(cvplan, net, ci, trap_level), max(runoff_live[inlet], 0.0))
                q <= 0.0 && continue
                _propagate_field_delta!(runoff_live, footprint_set, tstruct, full_traps_int, inlet,  -q)
                _propagate_field_delta!(runoff_live, footprint_set, tstruct, full_traps_int, outlet,  q)
            end
        end

        # --- NBS: the layers' current output replaces the oblivious throughput (as in the solver)
        for nb in net.nbs
            A      = Float64(length(nb.footprint))
            layers = get(nbs_layers, nb.id, zeros(Float64, length(nb.system.layers)))
            # measured on the untouched oblivious field, so stacked placements stay independent
            O_0, exit_shares, _ = _nbs_baseline(nb, tstruct, runoff_orig, precip)
            diff = _nbs_terrain_outflow(nb.system, layers, A, nb.n_terrain) - O_0
            for (c, sh) in zip(nb.link.exits, exit_shares)
                _propagate_field_delta!(runoff_live, footprint_set, tstruct, full_traps_int,
                                        LI[c], diff * sh)
            end
            # @@@ pond shares stay in the trap under the footprint (field undefined there); if
            #     that trap is full they should continue along its spill path
            for (o, l) in zip(nb.outlets, _piped_layers(nb.system, nb.n_terrain))
                _propagate_field_delta!(runoff_live, footprint_set, tstruct, full_traps_int,
                                        LI[o], first(_layer_fluxes(nb.system.layers[l], layers[l], A)))
            end
        end
    end

    # surface flow is undefined inside NBS footprints — mark them, after all propagation
    for c in footprint_set
        runoff_live[c] = NaN
    end
    return runoff_live
end

# Carry a flow change `amount` downstream from `start_lin` with the solver's per-cell rule (as
# `_track_flow!`): each cell's value changes by what arrives, and what passes on is
# `max(V+d,0) - max(V,0)` — spare infiltration absorbs added water, a removal stops at zero flow.
# Stops at an NBS footprint (it intercepts surface flow), when nothing is left, or where the
# route ends.  The route (`flow_path_from`) passes through full traps via their spill paths.
function _propagate_field_delta!(field, footprint_set, tstruct, full_traps_int::Vector{Int},
                                 start_lin::Int, amount::Float64)
    amount == 0.0 && return
    paths, _ = flow_path_from(tstruct, start_lin; full_traps = full_traps_int)
    for seg in paths, c in seg
        (c in footprint_set) && return
        V = field[c]
        field[c] = V + amount
        amount = max(V + amount, 0.0) - max(V, 0.0)
        amount == 0.0 && return
    end
    return
end
