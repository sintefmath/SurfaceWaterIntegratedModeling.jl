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

    # expand precipitation/infiltration to grids, prepared as fill_sequence does
    precip = precipitation isa Real ? fill(Float64(precipitation), gridres) : Float64.(precipitation)
    infil  = infiltration  isa Real ? fill(Float64(infiltration),  gridres) : Float64.(infiltration)
    _nbs_zero_infiltration!(infil, nbs)

    # the network-oblivious baseline field, fully routed by terrain. `runoff_live`
    # accumulates the network corrections; `runoff_orig` keeps the untouched
    # oblivious field. Culvert draws read the *live* flow (what is actually there,
    # matching the solver's min(capacity, current)); what an NBS removes at its exits
    # is measured on the *oblivious* field, so a downstream placement is never
    # contaminated by an upstream correction.
    runoff_live, = watercourses(tstruct, full_traps; precipitation = precip, infiltration = infil)
    runoff_orig  = copy(runoff_live)

    # build the network components: topology, cached NBS cell lists, culvert owners.
    # `z_vol_tables` depends only on `tstruct`, so a caller sampling many timepoints
    # should compute it once and pass it via `zvt`; otherwise it is built here.
    z_vol_tables = zvt !== nothing ? zvt : _compute_z_vol_tables(tstruct)
    seeds = _dyn_seeds(tstruct, dyn_traps, DynCulvert[])
    comps = setup_network(tstruct, findall(full_traps);
                          dyn_coords = seeds, culverts = culverts, nbs = nbs)

    footprint_set = Set{Int}(vcat([n.footprint for n in nbs]...))
    # what every flow change propagates through: the live field, the unchanged one (a full
    # trap's own output is read from it), each crossed full trap's net output so far
    fc = (; field = runoff_live, field0 = runoff_orig, trap_out = Dict{Int,Float64}(), footprint_set,
          tstruct, full = findall(full_traps))

    # the solver's baseline, which sets how an NBS splits its outflow over its exits: there,
    # network-covered full traps do not spill over the terrain (their spill reaches an NBS as draw)
    runoff_solver = nothing
    if !isempty(nbs)
        covered = _covered_of(comps, tstruct)
        solver_full = [full_traps[t] && !(t in covered) for t in eachindex(full_traps)]
        runoff_solver, = watercourses(tstruct, solver_full; precipitation = precip,
                                      infiltration = infil)
    end

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
                _propagate_field_delta!(fc, inlet,  -q)
                _propagate_field_delta!(fc, outlet,  q)
            end
        end

        # --- NBS: the layers' current output replaces the oblivious throughput (as in the solver)
        for nb in net.nbs
            A      = Float64(length(nb.footprint))
            layers = get(nbs_layers, nb.id, zeros(Float64, length(nb.system.layers)))
            # at each exit, replace what this field's footprint delivers there (measured on the
            # untouched field, so stacked placements stay independent) by the solver's share of
            # the layers' outflow
            O_t = _nbs_terrain_outflow(nb.system, layers, A, nb.n_terrain)
            _, exit_shares, _ = _nbs_baseline(nb, tstruct, runoff_solver, precip)
            delivered = _exit_throughput(nb.link, runoff_orig)
            for (c, sh, w) in zip(nb.link.exits, exit_shares, delivered)
                _propagate_field_delta!(fc, LI[c], O_t * sh - w)
            end
            # @@@ pond shares stay in the trap under the footprint (field undefined there); if
            #     that trap is full they should continue along its spill path
            for (o, l) in zip(nb.outlets, _piped_layers(nb.system, nb.n_terrain))
                _propagate_field_delta!(fc, LI[o], first(_layer_fluxes(nb.system.layers[l], layers[l], A)))
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
# The route (`flow_path_from`) jumps each full trap it meets to its spill point; the trap passes
# on change by the same rule, with its net output as V (`_through_trap!`).  Stops at an NBS
# footprint (it intercepts surface flow), when nothing is left, or where the route ends.  `fc`
# holds the field and what the propagation needs (see `network_watercourses`).
function _propagate_field_delta!(fc, start_lin::Int, amount::Float64)
    amount == 0.0 && return
    paths, fulls = flow_path_from(fc.tstruct, start_lin; full_traps = fc.full)
    traps = Iterators.Stateful(fulls)
    # a start inside the first full trap has no leading segment: that trap comes first
    if !isempty(fulls) && start_lin in fc.tstruct.footprints[fulls[1]]
        amount = _through_trap!(fc, popfirst!(traps), amount)
    end
    for seg in paths
        for c in seg
            (amount == 0.0 || c in fc.footprint_set) && return
            V = fc.field[c]
            fc.field[c] = V + amount
            amount = max(V + amount, 0.0) - max(V, 0.0)
        end
        isempty(traps) || (amount = _through_trap!(fc, popfirst!(traps), amount))
    end
    return
end

# Pass change `d` through full trap `t`, which acts as one cell whose value is its net output D:
# it passes on `max(D+d,0) - max(D,0)`, as the solver's full trap spills `inflow - footprint loss`
# only when positive.  A trap short of water (D < 0, e.g. just after the rain) takes added water
# first.  D starts from the unchanged field and accumulates over the changes through `t`.
function _through_trap!(fc, t, d)
    D = get!(() -> _trap_net_output(fc.field0, fc.tstruct, t), fc.trap_out, t)
    fc.trap_out[t] = D + d
    return max(D + d, 0.0) - max(D, 0.0)
end

# Full trap `t`'s net output in `field`: what reaches its bottom cells less the infiltration room
# left on its submerged footprint, as base SWIM's `_compute_trap_output`, but not clamped at 0.
function _trap_net_output(field, tstruct, t)
    sp, leaves = tstruct.spillpoints[t].elevation, tstruct.lowest_subtraps_for[t]
    net = 0.0
    for c in tstruct.footprints[t]
        tstruct.topography[c] < sp && (net += min(field[c], 0.0))
        isempty(Graphs.outneighbors(tstruct.flowgraph, c)) && tstruct.regions[c] in leaves &&
            (net += max(field[c], 0.0))
    end
    return net
end
