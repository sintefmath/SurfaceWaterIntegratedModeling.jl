# Traces the dynamic network over the terrain: from each seed cell, follow the flow
# path to the next trap and repeat, so seeds grow into a set of connected components.
#
# @@@ Contract: NBS footprints should have zero infiltration (avoids special-casing in
#     fill_sequence).

# ----------------------------------------------------------------------------
"""
    setup_network(tstruct, full_traps; dyn_coords, culverts, nbs)

Trace the dynamic-flow network over the terrain and return it as a
`Vector{DynNetwork}` (one per connected component).

Seeds are the `dyn_coords`, both endpoints of every culvert, and the NBS outflow
cells; from each seed the flow path is followed trap-to-trap until it leaves the
domain or reaches a trap already in the network.

# Arguments
- `tstruct`: the `TrapStructure` holding terrain, flow graph, traps and spillpoints.
- `full_traps`: indices of traps currently at capacity (they spill onward).
- `dyn_coords`, `culverts`, `nbs`: the dynamic seeds to trace from.

# Returns
`Vector{DynNetwork}` — the traced network split into connected components.  `culverts` and
`nbs` are held by reference on the components, not copied.

!!! note
    Despite the name, this **mutates `nbs`**: [`_prepare_nbs!`](@ref) sets each placement's
    `id` and terrain `link`.  `tstruct`, `full_traps` and `dyn_coords` are read only.
"""
function setup_network(tstruct, full_traps;
                       dyn_coords::Vector{CartesianIndex{2}}=CartesianIndex{2}[],
                       culverts::Vector{DynCulvert}=DynCulvert[],
                       nbs::Vector{DynNBSPlacement}=DynNBSPlacement[])

    # Input validation (no culverts or dyn_coords inside NBS footprints)
    _validate_network_inputs(tstruct, dyn_coords, culverts, nbs)

    _prepare_nbs!(nbs, tstruct)

    # seeds for tracking dynamic flows are the dyn_coords, BOTH culvert endpoints (the outlet
    # so its downstream chain is traced, the inlet so its host trap is an evolving node the
    # culvert can draw from — not just fill statically), and every cell an NBS emits water at
    seeds = union(Set(dyn_coords),
                  Set([culvert.inlet  for culvert in culverts]),
                  Set([culvert.outlet for culvert in culverts]),
                  Set(vcat([n.outlets for n in nbs]...)),
                  Set(vcat([n.link.ponds for n in nbs]...)),
                  Set(vcat([n.link.exits for n in nbs]...)))

    # Generate map to keep track of all existing network elements on terrain
    pathmap = Dict{Int, Int}() # terrain cell -> flow path index

    # Track all seeds and build up corresponding network of paths and traps.  Seed order is
    # immaterial: a path's source is its `departure_point` (not `cells[1]`), so a seed whose
    # first cell is already claimed just becomes a zero-length connector, and a chain reaching
    # an already-traced trap stops there (keeping that trap's single spill path).
    monolithic_network = DynNetwork(culverts, nbs)
    foreach(s -> _grow_network_from_seed!(monolithic_network, pathmap, s, tstruct, full_traps),
            seeds)

    # Reachability counts for dynamic-membership tracking; the split copies them per component.
    init_in_counts!(monolithic_network)

    # Split network in connected components
    networks = split_network_into_connected_components(monolithic_network, tstruct)

    return networks
end

# ----------------------------------------------------------------------------
"""
    _trace_to_next_trap(tstruct, start_cell, full_traps) -> (path, trap_ix)

Follow the flow path from `start_cell` to the trap it ends in.  `trap_ix == 0` means the flow
exits the domain.

The end trap is the highest *full* supertrap over the terminal region — a full trap pools with
its siblings, so water arriving there is really arriving at the merged parent — or the
lowest-level trap when none are full.  The target's own footprint cells are stripped from
`path`, since they are trap, not flow path.  Nothing is mutated.
"""
function _trace_to_next_trap(tstruct, start_cell, full_traps)

    path = _trace_path(tstruct, start_cell)
    end_reg = tstruct.regions[path[end]]

    if end_reg <= 0
        return path, 0
    end

    # Determine the end trap.  This is the same as the uppermost full trap,
    # unless all its siblings are also full, in which case we go up to the next
    # supertrap.
    supertraps = tstruct.supertraps_of[end_reg]
    highest_full_supertrap = findlast(t -> t in full_traps, supertraps)
    if isnothing(highest_full_supertrap)
        # no full traps in the supertrap hierarchy, so we return the lowest-level trap
        return path, end_reg
    end
    
    # if this is not the highest level trap, and all its siblings are full, return the
    # trap above it.  Otherwise, return the highest full trap.
    target_trap =
        (highest_full_supertrap < length(supertraps) &&
        all(s -> s in full_traps, subtrapsof(tstruct, supertraps[highest_full_supertrap + 1]))) ?
        supertraps[highest_full_supertrap + 1] :
        supertraps[highest_full_supertrap]

    # remove trap footprint cells from path, since they are not part of the flow path. 
    path = filter(c -> !(c in tstruct.footprints[target_trap]), path)
        
    return path, target_trap
end

# ----------------------------------------------------------------------------
"""
    _grow_network_from_seed!(network, pathmap, seed, tstruct, full_traps;
                             departing_trap_ix = 0) -> nothing

Trace from `seed` and attach what it finds to `network`: follow the flow path to the next trap,
add that trap and the connecting path, and repeat from its spillpoint until the flow leaves the
domain or reaches a trap already present.  Stopping at a present trap is deliberate — its
downstream is already represented, so only the incoming connector is attached rather than
re-tracing and overwriting that trap's single spill path.

# Arguments
- `network`: **mutated** — grows by the traced traps, flow paths, merges and `spill_path` links.
- `pathmap`: **mutated** — the terrain cell -> flow path index map, extended for each new path.
- `seed`: terrain cell to trace from.
- `tstruct`, `full_traps`: read-only terrain/trap state driving the trace.
- `departing_trap_ix`: the trap whose `spill_path` the first connector sets.  Used by
  [`grow_spill!`](@ref), where the seed is an existing trap's spillpoint; 0 at build, where
  seeds are `dyn_coords` / culvert / NBS cells.
"""
function _grow_network_from_seed!(network, pathmap, seed::CartesianIndex{2},
                                 tstruct, full_traps; departing_trap_ix::Int=0)
    LI = LinearIndices(tstruct.topography)
    CI = CartesianIndices(tstruct.topography)
    terminus = CartesianIndex(0, 0)
    while seed != terminus
        path, trap_ix = _trace_to_next_trap(tstruct, LI[seed], full_traps)

        trap_local = 0
        was_present = false
        # registering trap if it exists
        if trap_ix > 0
            # switch 'trap' from index into tstruct to index into network.traps
            trap_local = findfirst(dtrap -> dtrap.trap_ix == trap_ix, network.traps)
            was_present = !isnothing(trap_local)
            if !was_present
                # trap_ix is not yet in network.traps, so add it
                cv_in, cv_out = _intersecting_culverts(CI[tstruct.footprints[trap_ix]],
                                                       network.culverts)
                push!(network.traps, DynTrap(trap_ix, 0, cv_in, cv_out)) # @@ handle spill_path index later
                trap_local = length(network.traps)
            end
        end
        # Always register a connecting path (possibly zero cells) from the departing element
        # to the segment end; `seed` is its departure_point, valid even when `cells` is empty.
        isect_path, isect_cell =
            _update_pathmap!(pathmap, path, length(network.flow_paths)+1) # may empty path

        cv_in, cv_out, nbs_in = _intersecting_on_path(CI[path], network.culverts, network.nbs)
        target_trap = isect_path > 0 ? 0 : trap_local
        push!(network.flow_paths,
              DynFlowPath(CI[path], seed, target_trap, cv_in, cv_out, nbs_in, Tuple{Int,Int}[]))
        new_ix = length(network.flow_paths)

        (departing_trap_ix > 0) && (network.traps[departing_trap_ix].spill_path = new_ix)

        if isect_path > 0
            # from the intersection cell on, this path merges into the owning one
            isect_pt = findfirst(network.flow_paths[isect_path].cells .== CI[isect_cell])
            @assert isect_pt > 0 "Intersection point should be in the other path's cells."
            push!(network.flow_paths[isect_path].merges, (new_ix, isect_pt))
        end

        # if path ended in a full trap, it is spilling and we should keep on
        # tracing.  Otherwise, we are done.
        departing_trap_ix = trap_local

        if trap_ix == 0 || was_present
            seed = terminus
        else
            spoint = tstruct.spillpoints[trap_ix]
            if (trap_ix ∈ full_traps) && spoint.downstream_region_cell != spoint.current_region_cell
                seed = CI[spoint.downstream_region_cell]   # full trap spilling on: keep tracing
            else
                # Stop here.  A FULL trap with no distinct downstream cell spills straight out of
                # the domain: mark its spill_path with the out-of-domain sentinel (-1), otherwise
                # it reads as an unfilled transitory frontier (spill_path 0) and violates the
                # solver's three-state contract.  An unfilled frontier keeps spill_path 0.
                (trap_ix ∈ full_traps) && (network.traps[trap_local].spill_path = -1)
                seed = terminus
            end
        end

    end
end

# `_grow_network_from_seed!`'s `departing_trap_ix` kwarg above is used by `grow_spill!`
# (in network_updating.jl), the live-grow counterpart of `detach_spill!`.

# ----------------------------------------------------------------------------
"""
    _update_pathmap!(pathmap, path, path_ix) -> (isect_path, isect_cell)

Claim `path`'s cells in `pathmap` for `path_ix`, truncating `path` where it meets a cell an
earlier path already owns — flow is deterministic, so from that cell on the two coincide and
the owner carries it.

**Mutates both** `pathmap` (cells registered) and `path` (truncated in place; truncating at the
first cell empties it, and the caller then makes a zero-length connector, its source preserved
in `departure_point`).  Returns the meeting point as `(isect_path, isect_cell)`, both `0` if
the path met nothing.
"""
function _update_pathmap!(pathmap, path::Vector{Int}, path_ix)
    # check if path intersects with any existing NBS or culvert
    isect_path, isect_cell = 0, 0
    for (ix, cell) in enumerate(path)
        if haskey(pathmap, cell)
            # Truncate where it meets the owning path (flow is deterministic); ix==1 empties
            # it -> caller makes a zero-length connector (source kept in departure_point).
            splice!(path, ix:lastindex(path))
            isect_path, isect_cell = pathmap[cell], cell
        end
    end
    # register path in pathmap
    for cell in path
        pathmap[cell] = path_ix
    end
    
    return isect_path, isect_cell
end

# ----------------------------------------------------------------------------
"""
    _intersecting_culverts(footprint, culverts) -> (cv_in, cv_out)

The culverts whose inlet / outlet lies inside a trap `footprint`, as id lists — a trap is a pool,
so there is no along-path position to record (contrast `_intersecting_on_path`).
"""
function _intersecting_culverts(footprint, culverts)
    cv_in  = [ix for (ix, cv) in enumerate(culverts) if cv.inlet  in footprint]
    cv_out = [ix for (ix, cv) in enumerate(culverts) if cv.outlet in footprint]
    return cv_in, cv_out
end

# ----------------------------------------------------------------------------
"""
    _intersecting_on_path(path, culverts, nbs) -> (cv_in, cv_out, nbs_in)

The culvert inlets/outlets and NBS inflow cells on a flow `path`, each as `(id, pos)` with `pos`
the 1-based index of the cell within `path`.  The position matters: routing charges
infiltration up to that cell before applying the element, like a `merges` junction.  Every NBS
inflow cell crossed is registered; the first draws the flow and later ones see zero.
"""
function _intersecting_on_path(path, culverts, nbs)
    cv_in  = Tuple{Int,Int}[]
    cv_out = Tuple{Int,Int}[]
    nbs_in = Tuple{Int,Int}[]
    for (ix, culvert) in enumerate(culverts)
        p = findfirst(==(culvert.inlet),  path); p !== nothing && push!(cv_in,  (ix, p))
        p = findfirst(==(culvert.outlet), path); p !== nothing && push!(cv_out, (ix, p))
    end
    for (ix, n) in enumerate(nbs), (pos, cell) in enumerate(path)
        cell in n.link.inflow_cells && push!(nbs_in, (ix, pos))
    end
    return cv_in, cv_out, nbs_in
end

# ----------------------------------------------------------------------------
"""
    _validate_network_inputs(tstruct, dyn_coords, culverts, nbs) -> nothing

Error if an NBS footprint is empty or overlaps another (invariant I1), or if a culvert endpoint,
a `dyn_coord`, or a sink lies inside one (I2, agent/NBS_INVARIANTS.md): the layer model owns
all water there.  Nothing is mutated.
"""
function _validate_network_inputs(tstruct, dyn_coords, culverts, nbs)
    isempty(nbs) && return nothing
    nbs_footprints = Set(vcat([n.footprint for n in nbs]...))
    any(n -> isempty(n.footprint), nbs) && error("An NBS footprint is empty.")
    length(nbs_footprints) == sum(n -> length(unique(n.footprint)), nbs) ||
        error("NBS footprints overlap.")
    LI = LinearIndices(tstruct.topography)
    for culvert in culverts
        if LI[culvert.inlet] in nbs_footprints || LI[culvert.outlet] in nbs_footprints
            error("Culvert inlet or outlet is inside an NBS footprint.")
        end
    end
    for coord in dyn_coords
        if LI[coord] in nbs_footprints
            error("Dynamic flow coordinate is inside an NBS footprint.")
        end
    end
    for sink in tstruct.sinks
        if LI[sink] in nbs_footprints
            error("Sink cell $sink is inside an NBS footprint.  An NBS footprint must not " *
                  "contain a sink: water terminating there is invisible to the placement's " *
                  "ponding endpoints (its link's ponds cover trap bottoms only).")
        end
    end
    return nothing
end
