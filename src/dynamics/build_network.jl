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
`Vector{DynNetwork}` — the traced network split into connected components.
"""
function setup_network(tstruct, full_traps;
                       dyn_coords::Vector{CartesianIndex{2}}=CartesianIndex{2}[],
                       culverts::Vector{DynCulvert}=DynCulvert[],
                       nbs::Vector{DynNBSPlacement}=DynNBSPlacement[])

    # Input validation (no culverts or dyn_coords inside NBS footprints)
    _validate_network_inputs(tstruct, dyn_coords, culverts, nbs)

    # Stamp each placement's stable id = its position in this nbs vector, so its layer
    # storage persists across events (keyed by `nb.id` in the driver's nbs_state store).
    for (i, nb) in enumerate(nbs); nb.id = i; end

    # Compute inflow and outflow cells for each NBS footprint
    _compute_nbs_inflow_outflow_cells!(nbs, tstruct)
    
    # seeds for tracking dynamic flows are the dyn_coords, BOTH culvert endpoints (the outlet
    # so its downstream chain is traced, the inlet so its host trap is an evolving node the
    # culvert can draw from — not just fill statically), and all NBS outflow points
    seeds = union(Set(dyn_coords),
                  Set([culvert.inlet  for culvert in culverts]),
                  Set([culvert.outlet for culvert in culverts]),
                  Set(vcat([n.outlets for n in nbs]...)),
                  Set(vcat([n.internal_accumulation_cells for n in nbs]...)),
                  Set(vcat([n.footprint_outflow_cells for n in nbs]...)))

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
# Invert the flow graph: map each cell to the list of cells that drain into it.
function _map_inverse_flow(flowgraph)
    inv_flow = Dict{Int, Vector{Int}}()
    for cell in 1:Graphs.nv(flowgraph)
        ds = Graphs.outneighbors(flowgraph, cell)
        for d ∈ ds
            inv_flow[d] = get(inv_flow, d, Int[]) # initialize if not present
            push!(inv_flow[d], cell)
        end
    end
    return inv_flow
end

# ----------------------------------------------------------------------------
# Fill each placement's inflow / outflow / internal-accumulation cell lists from its footprint.
function _compute_nbs_inflow_outflow_cells!(nbs, tstruct)
    inv_flow = _map_inverse_flow(tstruct.flowgraph)
    CI = CartesianIndices(tstruct.topography)
    LI = LinearIndices(tstruct.topography)
    tbottoms = Set(LI[tstruct.trap_bottoms])
    
    for n in nbs
        n.footprint_outflow_cells = CI[_footprint_outflow_cells(tstruct, n.footprint)]
        n.footprint_inflow_cells = CI[_footprint_inflow_cells(inv_flow, n.footprint)]
        n.internal_accumulation_cells = CI[collect(intersect(Set(n.footprint), tbottoms))]
    end
end

# ----------------------------------------------------------------------------
# Follow the flow path from `start_cell` to the trap it ends in; return
# `(path, trap_ix)`, where `trap_ix == 0` means it exits the domain.  The end trap
# is the highest full supertrap over the terminal region (or the lowest-level trap
# if none are full); its footprint cells are stripped from `path`.
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
# `departing_trap_ix` seeds the trap whose spill_path the first connector sets (used by
# grow, where the seed is an existing trap's spillpoint; 0 at build, where seeds are
# dyn_coords / culvert / NBS cells).  The trace stops when it reaches a trap already in the
# network: its downstream is already represented, so we attach the incoming connector and
# stop, rather than re-tracing and overwriting that trap's single spill path.
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
                cv_in, cv_out, nbs_out =
                    _intersecting_culverts_and_nbs_outlets(CI[tstruct.footprints[trap_ix]],
                                                           network.culverts, network.nbs)
                push!(network.traps, DynTrap(trap_ix, 0, cv_in, cv_out)) # @@ handle spill_path index later
                trap_local = length(network.traps)
            end
        end
        # Always register a connecting path (possibly zero cells) from the departing element
        # to the segment end; `seed` is its departure_point, valid even when `cells` is empty.
        isect_path, isect_cell =
            _update_pathmap!(pathmap, path, length(network.flow_paths)+1) # may empty path

        cv_in, cv_out, nbs_in, nbs_out = _intersecting_on_path(CI[path], network.culverts, network.nbs)
        target_trap = isect_path > 0 ? 0 : trap_local
        push!(network.flow_paths,
              DynFlowPath(CI[path], seed, target_trap, cv_in, cv_out, nbs_in, nbs_out, Tuple{Int,Int}[]))
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
# Register `path`'s cells in `pathmap` under `path_ix`, truncating it where it first
# meets an already-claimed cell.  Returns `(isect_path, isect_cell)` of that meeting
# point (both 0 if none).
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
# The culvert inlets/outlets and NBS outlets falling inside a trap `footprint`,
# returned as bare id lists (a trap has no along-path position).
function _intersecting_culverts_and_nbs_outlets(footprint, culverts, nbs)
    cv_in = Int[]
    cv_out = Int[]
    nbs_out = Int[]
    for (ix, culvert) in enumerate(culverts)
        if culvert.inlet in footprint
            push!(cv_in, ix)
        end
        if culvert.outlet in footprint
            push!(cv_out, ix)
        end
    end
    for (ix, n) in enumerate(nbs)
        if any(outlet in footprint for outlet in n.outlets)
            push!(nbs_out, ix)
        end
    end

    return cv_in, cv_out, nbs_out
end

# ----------------------------------------------------------------------------
# For a flow path, the culvert inlets/outlets and NBS in/outlets that fall on it, each
# paired with the 1-based position of the matching cell within `path` (routing charges
# infiltration up to that cell, like a `merges` junction).  NBS inlets are the
# footprint-inflow boundary cells crossed by the path (every crossing is registered; the
# first draws the flow, later ones then see zero).  The trap form above stores bare ids
# since a trap has no along-path position.
function _intersecting_on_path(path, culverts, nbs)
    cv_in   = Tuple{Int,Int}[]
    cv_out  = Tuple{Int,Int}[]
    nbs_in  = Tuple{Int,Int}[]
    nbs_out = Tuple{Int,Int}[]
    for (ix, culvert) in enumerate(culverts)
        p = findfirst(==(culvert.inlet),  path); p !== nothing && push!(cv_in,  (ix, p))
        p = findfirst(==(culvert.outlet), path); p !== nothing && push!(cv_out, (ix, p))
    end
    for (ix, n) in enumerate(nbs)
        for (pos, cell) in enumerate(path)
            cell in n.footprint_inflow_cells && push!(nbs_in, (ix, pos))
        end
        for outlet in n.outlets
            p = findfirst(==(outlet), path)
            p !== nothing && (push!(nbs_out, (ix, p)); break)
        end
    end
    return cv_in, cv_out, nbs_in, nbs_out
end

# ----------------------------------------------------------------------------
# Error if any culvert endpoint or dyn_coord lies inside an NBS footprint.
function _validate_network_inputs(tstruct, dyn_coords, culverts, nbs)
    nbs_footprints = vcat([n.footprint for n in nbs]...)
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
end

# ----------------------------------------------------------------------------
# External cells the footprint drains into (downstream neighbours outside it).
function _footprint_outflow_cells(tstruct, footprint::Vector{Int})
    flowgraph = tstruct.flowgraph
    outflow_cells = Set{Int}()
    for cell in footprint
        ds, terminus = _downstream_cell(flowgraph, cell)   # utils.jl: (downstream cell, is-terminus?)
        if !terminus && !(ds in footprint)
            push!(outflow_cells, ds)
        end
    end
    # if there are no outflow cells, that is OK, but if upper NBS layers do not
    # have other, explicitly given outlets, this NBS system will become
    # submerged with no runoff.
    if isempty(outflow_cells)
        @warn "NBS placement footprint has no outflow cells."
    end
    return collect(outflow_cells)
end

# ----------------------------------------------------------------------------
# External cells that drain into the footprint (upstream neighbours outside it).
function _footprint_inflow_cells(inv_flow::Dict{Int, Vector{Int}},
                                 footprint::Vector{Int})
    all_inflow_cells = vcat([get(inv_flow, cell, Int[]) for cell in footprint]...)
    return collect(setdiff(Set(all_inflow_cells), Set(footprint)))
end

