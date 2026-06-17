export DynObject, DynFlowPath, DynTrap, DynCulvert

# Make generic baseclass for dynamic objects
abstract type DynObject end

# Dynamic flow path
struct DynFlowPath <: DynObject
    cells::Vector{CartesianIndex} # cells along the flow path

    # Target dynamic index (0 for out-of-domain)
    target_trap::Int

    # culvert inlets and outlets, each represented by a culvert ID
    culvert_inlets::Vector{Int}
    culvert_outlets::Vector{Int}
    
    # other flow paths that merge into this one
    merges::Vector{Int} 
end

# Dynamic trap
struct DynTrap <: DynObject
    trap_ix::Int # index of the trap in the spillanalysis structure

    # spill path
    spill_path::Int

    # culvert inlets and outlets, each represented by a culvert ID
    culvert_inlets::Vector{Int}
    culvert_outlets::Vector{Int}
end

# Dynamic culvert
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
struct DynNetwork
    flow_paths::Vector{DynFlowPath}
    traps::Vector{DynTrap}
    culverts::Vector{DynCulvert}
end

# todo: inclusion of culverts as separate argument
function setup_network(tstruct, dyn_coords, full_traps)
    LI = LinearIndices(tstruct.topography)
    _merge_networks([_subnetwork(tstruct, LI(c)) for c in dyn_coords])
end

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


function _merge_networks(networks::Vector{DynNetwork})
   # merge all networks that overlap into a single one.  Return vector of disjoint networks.
end
