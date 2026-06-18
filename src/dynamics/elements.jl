import Graphs

export DynObject, DynFlowPath, DynTrap, DynCulvert, DynNetwork, setup_network

# Make generic baseclass for dynamic objects
abstract type DynObject end

"""
        DynFlowPath(cells, target_trap, culvert_inlets, culvert_outlets, merges)

Represent a flow path over the terrain, represented by a sequence of grid cells.
The flow path may lead into a trap (target_trap), terminate in an intersecting
flow path (target_trap=0), or flow out of the domain (target_trap=0).

The flow path may also include culverts along the way.  Culvert inlets would
subtract water from the flow path, whereas culvert outlets would add water to it.
The infiltration capacity of each cell in the path is represented externally.

"""
struct DynFlowPath <: DynObject
    cells::Vector{CartesianIndex{2}} # cells along the flow path

    # Target trap index (0 for out-of-domain or intersection with another flow path)
    target_trap::Int

    # culvert inlets and outlets, each represented by a culvert ID
    culvert_inlets::Vector{Int}
    culvert_outlets::Vector{Int}

    # other flow paths that merge into this one
    merges::Vector{Int}
end

DynFlowPath(cells, target_trap) = DynFlowPath(cells, target_trap, Int[], Int[], Int[])
DynFlowPath(cells) = DynFlowPath(cells, 0, Int[], Int[], Int[])

"""
        DynTrap(trap_ix, spill_path, culvert_inlets, culvert_outlets)

Represent a trap in the terrain, identified by its index in the spillanalysis
structure.  Every trap must have a spill path (spill_path) that represents the
flow path that water would take when the trap's capacity is exceeded.

The trap may also have culvert inlets or outlets within its footprint, which
would add or subtract water from the trap respectively, depending on the current
water level in the trap.  The infiltration capacity of each cell in the trap is
represented externally.

"""
struct DynTrap <: DynObject
    trap_ix::Int # index of the trap in the spillanalysis structure

    # spill path
    spill_path::Int

    # culvert inlets and outlets, each represented by a culvert ID
    culvert_inlets::Vector{Int}
    culvert_outlets::Vector{Int}
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

# Network of dynamic objects
"""
         DynNetwork(flow_paths, traps, culverts)

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
struct DynNetwork
    flow_paths::Vector{DynFlowPath}
    traps::Vector{DynTrap}
    culverts::Vector{DynCulvert}
end

DynNetwork() = DynNetwork(DynFlowPath[], DynTrap[], DynCulvert[])

# todo: inclusion of culverts as separate argument
"""
      setup_network(tstruct, dyn_coords, full_traps)

Given a spillanalysis structure, a set of coordinates representing dynamic points
in the terrain, and a set of traps that are currently filled, build a network
of dynamic flow paths and traps that represent the flow of water through the terrain.  The flow paths are built starting from the dynamic points, and are connected
to traps as they lead into them.

TODO: For the moment, culverts are not supported. They will be included in a future update.

"""
function setup_network(tstruct, dyn_coords, full_traps)
    _merge_networks([_subnetwork(tstruct, c, full_traps) for c in dyn_coords])
end

"""
        _subnetwork(tstruct, coord, full_traps)

Create a simple line-like network associated with the complete flow path from a
point in the terrain.
"""
function _subnetwork(tstruct, coord::CartesianIndex, full_traps)
    LI = LinearIndices(tstruct.topography)

    paths, ftraps = flow_path_from(tstruct, LI[coord]; full_traps=full_traps)
    @assert abs(length(paths) - length(ftraps)) <= 1

    # The result is a single chain that strictly alternates between path segments
    # and traps in downstream order.  It begins with a trap when the start point
    # lies inside the first full trap (which then spills into the first path); this
    # is detected by checking whether the first path segment emanates from that
    # trap's spill point.
    starts_with_trap = !isempty(ftraps) &&
        (isempty(paths) ||
         paths[1][1] == tstruct.spillpoints[ftraps[1]].downstream_region_cell)

    # When the chain ends on a path, that path may still discharge into an unfilled
    # trap downstream; include it as the terminating trap.
    traps = collect(ftraps)
    ends_with_path = length(paths) > length(ftraps) ||
                     (starts_with_trap && length(paths) == length(ftraps))
    if ends_with_path
        tix = _unfilled_trap_at(tstruct, paths[end][end], full_traps)
        tix > 0 && push!(traps, tix)
    end

    return _build_network(paths, traps, starts_with_trap, tstruct)
end

# Return the index of the uppermost unfilled trap that `cell` drains into, or 0 if
# the cell drains out of the domain or only into already-full traps.
function _unfilled_trap_at(tstruct, cell::Int, full_traps)
    cur_reg = tstruct.regions[cell]
    cur_reg <= 0 && return 0
    unfilled = filter(x -> x ∉ full_traps, tstruct.supertraps_of[cur_reg])
    isempty(unfilled) ? 0 : minimum(unfilled)
end

# Build the DynNetwork from the alternating path/trap chain.  Both vectors are
# already in downstream order, so their local indices are simply 1:length, and
# `starts_with_trap` fixes the offset between a path and the trap it flows into:
# when the chain starts with a trap, path i flows into trap i+1 and trap i spills
# into path i; otherwise path i flows into trap i and trap i spills into path i+1.
function _build_network(paths, traps, starts_with_trap, tstruct)
    CI = CartesianIndices(tstruct.topography)
    link(i, off, n) = i + off <= n ? i + off : 0
    pt = starts_with_trap ? 1 : 0  # offset from a path to the trap it flows into

    dyn_paths = [DynFlowPath([CI[k] for k in paths[i]], link(i, pt, length(traps)))
                 for i in eachindex(paths)]
    dyn_traps = [DynTrap(traps[i], link(i, 1 - pt, length(paths)))
                 for i in eachindex(traps)]

    return DynNetwork(dyn_paths, dyn_traps, DynCulvert[])
end

"""
        _merge_networks(networks)

Given a vector of flow path networks, merge any that overlap into a single one.
Return a vector of disjoint networks.  Where flow paths share grid cells, the
later path is truncated at the shared cell and registered as a tributary (a
"merge") of the path it runs into, so that no grid cell is shared by multiple
flow paths.

"""
function _merge_networks(networks::Vector{DynNetwork})
    isempty(networks) && return networks

    # flatten into one pool with globally unique indices
    all_paths, all_traps, all_culverts = _combine_networks(networks)

    # truncate paths that share cells, registering each truncated path as a tributary
    _resolve_cell_overlaps!(all_paths, all_culverts)

    # _build_component: reconstruct each component as a self-contained DynNetwork
    # _path_components: group paths into disjoint connected components
    return [_build_component(all_paths, all_traps, all_culverts, ids)
            for ids in _path_components(all_paths, all_traps)]     
end

# Merge all networks into a single flat pool with globally unique indices
function _combine_networks(networks::Vector{DynNetwork})
    path_offsets    = [0; cumsum([length(n.flow_paths) for n in networks])[1:end-1]]
    trap_offsets    = [0; cumsum([length(n.traps)      for n in networks])[1:end-1]]
    culvert_offsets = [0; cumsum([length(n.culverts)   for n in networks])[1:end-1]]

    remap(idx, off)   = idx == 0 ? 0 : idx + off
    remap_vec(v, off) = isempty(v) ? copy(v) : v .+ off

    all_paths    = DynFlowPath[]
    all_traps    = DynTrap[]
    all_culverts = DynCulvert[]

    for (ni, net) in enumerate(networks)
        poff, toff, coff = path_offsets[ni], trap_offsets[ni], culvert_offsets[ni]
        for p in net.flow_paths
            push!(all_paths, DynFlowPath(copy(p.cells),
                                         remap(p.target_trap, toff),
                                         remap_vec(p.culvert_inlets, coff),
                                         remap_vec(p.culvert_outlets, coff),
                                         remap_vec(p.merges, poff)))
        end
        for t in net.traps
            push!(all_traps, DynTrap(t.trap_ix,
                                     remap(t.spill_path, poff),
                                     remap_vec(t.culvert_inlets, coff),
                                     remap_vec(t.culvert_outlets, coff)))
        end
        append!(all_culverts, net.culverts)
    end

    return all_paths, all_traps, all_culverts
end

# Truncate any flow path whose cells overlap with a previously-processed path.
# The truncated path is registered as a tributary (a "merge") of the primary
# (earlier) path that owns the shared cell, and culverts beyond the truncation
# point are dropped.  Trap connectivity is left untouched.
function _resolve_cell_overlaps!(all_paths, all_culverts)
    cell_owner = Dict{CartesianIndex{2}, Int}()  # grid cell → owning path index

    for pi in 1:length(all_paths)
        path = all_paths[pi]
        merge_pos = findfirst(cell -> haskey(cell_owner, cell), path.cells)

        if merge_pos !== nothing
            merge_into = cell_owner[path.cells[merge_pos]]
            kept       = path.cells[1:merge_pos-1]
            kept_set   = Set(kept)

            # Drop culverts whose cell falls beyond the truncation point
            inlets  = filter(c -> all_culverts[c].inlet  ∈ kept_set, path.culvert_inlets)
            outlets = filter(c -> all_culverts[c].outlet ∈ kept_set, path.culvert_outlets)

            all_paths[pi] = DynFlowPath(kept, 0, inlets, outlets, path.merges)

            primary = all_paths[merge_into]
            all_paths[merge_into] = DynFlowPath(primary.cells, primary.target_trap,
                primary.culvert_inlets, primary.culvert_outlets, [primary.merges; pi])

            for cell in kept
                cell_owner[cell] = pi
            end
        else
            for cell in path.cells
                cell_owner[cell] = pi
            end
        end
    end
end

# Group path indices into disjoint connected components using union-find.
# Two paths are in the same component if one merges into the other, or if
# they are connected through a shared trap (path → trap → spill_path).
function _path_components(all_paths, all_traps)
    n      = length(all_paths)
    parent = collect(1:n)

    find_root(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    function unite!(x, y)
        rx, ry = find_root(x), find_root(y)
        rx != ry && (parent[rx] = ry)
    end

    for pi in 1:n
        for mp in all_paths[pi].merges
            unite!(pi, mp)
        end
        t = all_paths[pi].target_trap
        if t > 0
            sp = all_traps[t].spill_path
            sp > 0 && unite!(pi, sp)
        end
    end

    components = Dict{Int, Vector{Int}}()
    for pi in 1:n
        push!(get!(components, find_root(pi), Int[]), pi)
    end
    return collect(values(components))
end

# Return (sorted_path_ids, sorted_trap_ids) in upstream-to-downstream order,
# using a joint topological sort of the path/trap DAG via Graphs.
function _topological_order(global_path_ids, global_trap_ids, all_paths, all_traps)
    np = length(global_path_ids)
    path_node = Dict(gpi => i      for (i, gpi) in enumerate(global_path_ids))
    trap_node = Dict(gti => np + i for (i, gti) in enumerate(global_trap_ids))

    # Every path/trap referenced from within this component is itself part of the
    # component: _path_components unions paths via their merges and via
    # path→trap→spill_path, and _build_component collects every target trap.  So a
    # nonzero reference is guaranteed to resolve here (asserted, not tested).
    g = Graphs.SimpleDiGraph(np + length(global_trap_ids))
    for (li, gpi) in enumerate(global_path_ids)
        p = all_paths[gpi]
        if p.target_trap > 0
            @assert haskey(trap_node, p.target_trap)
            Graphs.add_edge!(g, li, trap_node[p.target_trap])
        end
        for m in p.merges
            @assert haskey(path_node, m)
            Graphs.add_edge!(g, path_node[m], li)
        end
    end
    for (li, gti) in enumerate(global_trap_ids)
        sp = all_traps[gti].spill_path
        if sp > 0
            @assert haskey(path_node, sp)
            Graphs.add_edge!(g, np + li, path_node[sp])
        end
    end

    # The flow graph is assumed acyclic (water flows strictly downstream); a cycle
    # is a programming error, and topological_sort_by_dfs throws if one is present.
    order = Graphs.topological_sort_by_dfs(g)
    return ([global_path_ids[i]      for i in order if i <= np],
            [global_trap_ids[i - np] for i in order if i >  np])
end

# Reconstruct one DynNetwork from a set of global path indices, remapping all
# internal references to local (1-based) indices.
function _build_component(all_paths, all_traps, all_culverts, global_path_ids)
    path_set = Set(global_path_ids)

    # Collect traps targeted by paths in this component (the normal case), plus
    # leaf traps that no path flows into but which spill into a path here —
    # this arises when the starting point falls inside a full trap's footprint,
    # leaving that trap with no upstream flow path but a spill_path into the network.
    global_trap_ids = Int[]
    for gpi in global_path_ids
        t = all_paths[gpi].target_trap
        t > 0 && t ∉ global_trap_ids && push!(global_trap_ids, t)
    end
    for (ti, trap) in enumerate(all_traps)
        ti ∉ global_trap_ids && trap.spill_path ∈ path_set && push!(global_trap_ids, ti)
    end

    global_path_ids, global_trap_ids =
        _topological_order(global_path_ids, global_trap_ids, all_paths, all_traps)

    global_culvert_ids = unique(vcat(
        [c for gpi in global_path_ids
           for c in vcat(all_paths[gpi].culvert_inlets, all_paths[gpi].culvert_outlets)],
        [c for gti in global_trap_ids
           for c in vcat(all_traps[gti].culvert_inlets,  all_traps[gti].culvert_outlets)]))

    path_map    = Dict(gpi => lpi for (lpi, gpi) in enumerate(global_path_ids))
    trap_map    = Dict(gti => lti for (lti, gti) in enumerate(global_trap_ids))
    culvert_map = Dict(gci => lci for (lci, gci) in enumerate(global_culvert_ids))

    local_paths = [DynFlowPath(
        all_paths[gpi].cells,
        get(trap_map, all_paths[gpi].target_trap, 0),
        [culvert_map[c] for c in all_paths[gpi].culvert_inlets],
        [culvert_map[c] for c in all_paths[gpi].culvert_outlets],
        [path_map[m] for m in all_paths[gpi].merges if m ∈ path_set]
    ) for gpi in global_path_ids]

    local_traps = [DynTrap(
        all_traps[gti].trap_ix,
        get(path_map, all_traps[gti].spill_path, 0),
        [culvert_map[c] for c in all_traps[gti].culvert_inlets],
        [culvert_map[c] for c in all_traps[gti].culvert_outlets]
    ) for gti in global_trap_ids]

    return DynNetwork(local_paths, local_traps, [all_culverts[gci] for gci in global_culvert_ids])
end
