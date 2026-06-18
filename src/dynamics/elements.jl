export DynObject, DynFlowPath, DynTrap, DynCulvert

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
struct
DynFlowPath <: DynObject cells::Vector{CartesianIndex} # cells along the flow
path

    # Target dynamic index (0 for out-of-domain or intersection with another flow path)
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
struct DynTrap <: DynObject trap_ix::Int # index of
the trap in the spillanalysis structure

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
cell indices of the inlet and outlet.  The culvert has a limited capacity, which
is determined by its internal parameters as well as the dynamic water levels at
the inlet and outlet.

Internal parameters include the radius (r), discharge coefficient (Cd), entrance
loss coefficient (Ke), friction loss coefficient (Kf), and weir coefficient
(Cw).  Which parameters are used in the computation of the flow through the
culvert depends on the water levels at the inlet and outlet.

""" 
struct DynCulvert <: DynObject
    inlet::CartesianIndex # cell index of the culvert inlet
    outlet::CartesianIndex # cell index of the culvert outlet
    
    r::Float64 # radius of the culvert
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
    LI = LinearIndices(tstruct.topography)
    _merge_networks([_subnetwork(tstruct, LI(c)) for c in dyn_coords])
end

"""
        _subnetwork(tstruct, pt_ix, full_traps)

Create a simple line-like network associated with the complete flow path from a
point in the terrain.
"""
function _subnetwork(tstruct, pt_ix, full_traps)
    dyn_paths = Vector{DynFlowPath}()
    dyn_traps = Vector{DynTrap}()

    # compute basic flowpath
    paths, downstream_full_traps = flow_path_from(tstruct, coord, full_traps)
    @assert length(paths) == length(downstream_full_traps) ||
            length(paths) == length(downstream_full_traps) + 1
    
    for i in 1:length(paths)

        if length(paths[i]) > 0
            push!(dyn_paths, DynFlowPath(paths[i], length(dyn_paths) + 1))
        else
            @assert i == 0 "Only first path from flow_path_from can be empty"
        end
        
        if i <= length(downstream_full_traps)
            next_path_ix = i < length(paths) ? length(dyn_paths) + 1 : 0
            push!(dyn_traps, DynTrap(downstream_full_traps[i], next_path_ix))
        else
            @assert i == length(paths) "Only last path from flow_path_from can be without a trap"
            # the last path did not lead into a filled trap.  But perhaps it lead into an
            # unfilled one?
            cur_reg = tstruct.regions[paths[end][end]]
            if cur_reg > 0
                unfilled_supertraps = filter(x -> x ∉ full_traps, tstruct.supertraps_of[cur_reg])
                if !isempty(unfilled_supertraps)
                    push!(dyn_traps, DynTrap(minimum(unfilled_supertraps), 0))
                end
            end
        end
    end
end

"""
        _merge_networks(networks)

Given a vector of flow path networks, merge any that overlap into a single one.
Return a vector of disjoint networks.  Intersecting flow paths should be
properly merged and truncated, so that no grid cell is shared by multiple flow
paths, and traps at the end of truncated paths are connected to the remaining
path instead.

Returns a vector of disjoint networks, where each network is a tree-like
structure of flow paths and traps.  The networks are disjoint in the sense that
no flow path in one network shares any grid cell with a flow path in another
network.

"""
function _merge_networks(networks::Vector{DynNetwork})
   # merge all networks that overlap into a single one.  Return vector of disjoint networks.
end
