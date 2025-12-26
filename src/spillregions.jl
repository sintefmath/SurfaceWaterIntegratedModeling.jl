import Graphs
import OffsetArrays: OffsetArray

export spillregions, update_spillregions!

_SPILLFIELD_INACTIVE_FLAG = -10; # could in principle be any integer outside the
                                 # range [-2:7]
# ----------------------------------------------------------------------------
"""
    spillregions(spillfield)

Identify all spill regions to be derived from a given [`spillfield`](@ref).

A spill region is a group of cells that flow into the same trap (here: the
lowest-level traps).

Returns a Matrix{Int} of the same shape as `spillfield`, where all cells with
the same integer value belong to the same spill region.  Spill regions that 
exit the domain are assigned negative region numbers.

# Arguments
- `spillfield::Matrix{Int}`: the matrix containing the spillfield already 
                             computed from the [`spillfield`](@ref) function.
- `usediags` : if `true`, diagonal connections between cells will also be considered
- `tiling::Union{Tuple{Int, Int}, Nothing}`: 
      tuple specifying number of 'tiles' to subdivide surface in for parallel 
      processing.  Default is (1,1), which means the whole surface is treated
      as a single tile (no parallel processing).
- `cut_edges::Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}}`:
      dictionary specifying edges that should be cut (i.e., no flow allowed
      across these edges).  Keys are CartesianIndices of grid cells, and values
      are Vectors of CartesianIndices of neighboring grid cells to which flow
      is blocked.  This dict should also have been used when generating the spillfield.
- `culverts::Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}}`:
      vector specifying culverts as pairs of CartesianIndices of grid cells
      that are connected by the culvert.  The grid cells should be ordered so that
      the first cell in the tuple is at higher elevation than the second cell.
See also [`spillfield`](@ref), [`update_spillregions!`](@ref).
"""
function spillregions(spillfield::Matrix{Int8};
                      usediags::Bool=true,
                      tiling=nothing,
                      cut_edges::Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}}=
                          Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}}(),
                      culverts::Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}}=
                          Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}}())
    
    regions = similar(spillfield, Int64)
    edges = Set{Tuple{Int, Int}}()
    domain = Domain2D(1:size(spillfield,1), 1:size(spillfield,2))
    cut_edges_bidir = _prepare_cut_edges_bidir(cut_edges, size(spillfield))
    
    if tiling == nothing
        _process_domain!(regions, edges, spillfield,
                         domain=domain,
                         usediags=usediags,
                         cut_edges=cut_edges_bidir,
                         culverts=culverts)
    else
        @assert false "proper parallel implementation of spillregions not yet complete."
        # @@ TODO What needs to be done is ensure that edges crossing tile boundaries
        # are registered (which is not the case for the code below).  Moreover,
        # when fixing boundary seams, we need to make sure that cut_edges are considered.
        # Also, there's issues around sharing of `edges` across threads.
        
        tiles, splits = tiledomain(domain, tiling[1], tiling[2])
        
        Threads.@threads for i = 1:length(tiles)
            _process_domain!(regions, edges, spillfield,
                             domain=tiles[i],
                             usediags=usediags,
                             cut_edges=cut_edges_bidir,
                             culverts=culverts)
        end

        xsplits = splits[1]
        ysplits = splits[2]
        # @@ TODO: needs to process `edges` too
        _fix_boundary_seams!(regions, spillfield, usediags,
                             xsplits[2:end-1], ysplits[2:end-1])
    end
    # Identifying exit nodes (nodes whose stream direction points out of the
    # domain, or bottom nodes lying right on the boundary).  
    enodes = _find_exit_nodes(spillfield)

    # Determine the spill regions of exit nodes
    enoderegs = unique(regions[enodes])

    # give consistent numbering to regions: from 1 and upwards for regions
    # associated with a trap, and from -1 downwards for regions spilling out of
    # the domain.
    _renumerate_regions!(regions, exitregions=enoderegs)

    # construct spillfield graph

    # remove self-loops and edges connecting nodes within the same region
    edge_filt = filter(a->(a[1] != a[2]), edges)
    
    g = Graphs.SimpleDiGraph([Graphs.SimpleEdge{Int64}((e[1], e[2])) for e in edge_filt])
    
    return regions, g
    
end

# ----------------------------------------------------------------------------
"""
    update_spillregions!(regions, spillfield, domain; usediags=true, 
                         return_region_reindex=false)

Update a previously computed spillregion field inside a limited domain (where
the spillfield has presumably changed).

# Arguments
- `regions::Matrix{Int}`        : the previously computed spillregion field
- `spillfield::Matrix{Int8}`    : the locally updated spill field
- `domain::Domain2D`            : domain within which the spillfield has changed, 
                                  and thus regions need to be updated
- `usediags{Bool}`              : if `true`, diagonal connections between cells 
                                  will also be considered
- `return_region_reindex{Bool}' : compute and return correspondence between old 
                                  and new region indices as a 2-column matrix, 
                                  first column represents old region numbers and 
                                  second column the new ones.  Note that there 
                                  may be both one-to-one, one-to-many and
                                  many-to-one correspondences.

See also [`spillregions`](@ref).
"""
function update_spillregions!(regions::Matrix{Int64}, spillfield::Matrix{Int8},
                              domain::Domain2D;
                              usediags=true, return_region_reindex=false)
    @assert false "This function needs to be updated to handle cut_edges properly, as well as also updating the flowgraph."
    
    # expand domain by one gridcell in all direction, since the gridcell closest
    # to the modified area may also have changed spill direction (which depends
    # on the immediate neighbors)
    domain2 = Domain2D(
        max(domain.xrange[1]-1, 1):min(domain.xrange[end]+1, size(regions, 1)),
        max(domain.yrange[1]-1, 1):min(domain.yrange[end]+1, size(regions, 2)))

    reindex = _process_extended_domain!(regions, spillfield, domain2, usediags,
                                        return_region_reindex)

    # Identifying exit node and renumbering
    enodes = _find_exit_nodes(spillfield)
    enoderegs = unique(regions[enodes])
    _renumerate_regions!(regions, exitregions=enoderegs, reindex=reindex)

    return reindex
end

# ----------------------------------------------------------------------------
function _process_extended_domain!(regions::Matrix{Int64},
                                   spillfield::Matrix{Int8},
                                   domain::Domain2D,
                                   usediags::Bool,
                                   return_region_reindex::Bool)
    (xmin, xmax) = (domain.xrange[1], domain.xrange[end])
    (ymin, ymax) = (domain.yrange[1], domain.yrange[end])
    boundary_regs = unique(vcat(regions[[xmin, xmax], ymin:ymax][:],
                                regions[xmin:xmax, [ymin, ymax]][:]))

    mask = fill(false, maximum(regions))
    (min_reg, max_reg) = extrema(regions)
    mask = OffsetArray(fill(false, max_reg - min_reg + 1), min_reg-1)

    mask[boundary_regs] .= true
    mask = mask[regions]; # mask now has same shape as 'regions', where regions
                          # touched by domain boundary are flagged as 'false'
    mask[xmin:xmax, ymin:ymax] .= true; # flag inside of domain as 'false' too.

    spillfield_modif = copy(spillfield)
    spillfield_modif[.!mask] .= _SPILLFIELD_INACTIVE_FLAG

    edges = Vector{Tuple{Int, Int}}()

    # Identify grid edges connecting neighbor 'trap bottom' cells.
    _flat_zone_connecting_edges!(edges, spillfield_modif, usediags)

    # Identify grid edges constituting "streamlines"
    _spillfield_flow_edges!(edges, spillfield_modif)

    # Identify connected components
    components = _determine_connected_components(edges, length(regions))

    # filling in spill regions
    max_reg = maximum(regions)
    
    reindex = []

    if return_region_reindex
        # find first cell index of each unique region
        fix = unique(z->regions[z], 1:length(regions))
        oldreg = regions[fix]; # equivalent to unique(regions)
    end
    
    for i in 1:length(components)

        new_regnum = i + max_reg
        if return_region_reindex
            regs = unique(regions[components[i]])
            push!(reindex, hcat(regs, fill(new_regnum, length(regs))))
        end
        
        regions[components[i]] .= new_regnum
    end

    if return_region_reindex

        newreg = regions[fix]; # region numbers after updating
        
        return unique(vcat(hcat(oldreg, newreg), reindex...), dims=1)
    end
end

# ----------------------------------------------------------------------------
function _determine_connected_components(edges, num_nodes)

    # make temporary map to avoid SimpleGraph filling in additional nodes in-between

    actives = BitArray(undef, num_nodes)
    actives .= false; # 'undef' seems to set 'actives' to false above, but just
                      # in case, we do it explicitly here
    for e in edges
        actives[e[1]] = true
        actives[e[2]] = true
    end
    unique_reg_ixs = findall(actives)

    # The following line does the same as the code block above, but is much
    # slower due to the call to 'unique'.  The code line is left in place as a
    # clarification.
    #unique_reg_ixs = unique(collect(Iterators.flatten(edges)))

    map = zeros(Int, maximum(unique_reg_ixs))
    map[unique_reg_ixs] = 1:length(unique_reg_ixs)
    g  = Graphs.SimpleGraph([Graphs.SimpleEdge{Int64}((map[e[1]], map[e[2]]))
                             for e in edges])
    components = Graphs.connected_components(g)

    # map back to the original indices used for edges

    for i in 1:length(components)
        components[i] = unique_reg_ixs[components[i]]
    end
    
    return components
end

# ----------------------------------------------------------------------------
function _process_domain!(regions::Matrix{Int64},
                          edges::Set{Tuple{Int, Int}},
                          spillfield::Matrix{Int8};
                          usediags::Bool=true,
                          domain=nothing,
                          cut_edges::Set{Tuple{Int, Int}},
                          culverts::Vector{Tuple{CartesianIndex{2}, CartesianIndex{2}}})

    spilldomain = view(spillfield, domain.xrange, domain.yrange)
    regionsdomain = view(regions, domain.xrange, domain.yrange)

    # Identify grid edges connecting neighbor 'trap bottom' cells.
    all_connections = Set{Tuple{Int, Int}}()
    _flat_zone_connecting_edges!(all_connections, spilldomain, usediags)

    # eliminate the edges where barriers split flat zones into separate regions
    setdiff!(all_connections, cut_edges)
    
    # Identify grid edges constituting "streamlines"
    _spillfield_flow_edges!(edges, spilldomain)

    # Redirect flow across culverts: If the culvert is between p1 to p2, insert
    # an edge between these two points, directed from highest (p1) to lowest
    # (p2) elevation, and remove any edges flowing out from p1.
    LI = LinearIndices(size(spilldomain))
    edges_to_remove = Set{Tuple{Int, Int}}()
    for c in culverts
        p1, p2 = LI[c[1]], LI[c[2]]
        # identify any edges flowing out from p1
        push!(edges_to_remove, filter(e->e[1] == p1, edges)...)
        # add edge from p1 to p2
        push!(edges, (p1, p2))
    end
    setdiff!(edges, edges_to_remove)

    # identify regions connected by streamlines
    union!(all_connections, edges)
    components = _determine_connected_components(all_connections, length(regions))
    
    # filling in spill regions
    regionsdomain .= 0 # masked areas do not have any spill regions, and will be
                       # left with region value 0 (masked, inactive)
    for i in 1:length(components)
        regionsdomain[components[i]] .= i
    end
end

# ----------------------------------------------------------------------------
function _prepare_cut_edges_bidir(cut_edges::Dict{CartesianIndex{2}, Vector{CartesianIndex{2}}}, 
                                  gridsize::Tuple{Int, Int})
    cut_edges_bidir = Set{Tuple{Int, Int}}()
    for (k, v) in cut_edges
        for neigh in v
            push!(cut_edges_bidir, (LinearIndices(gridsize)[k],
                                    LinearIndices(gridsize)[neigh]))
        end
    end
    return cut_edges_bidir
end


# ----------------------------------------------------------------------------
function _remap!(regions::AbstractArray{Int64}, fromto::Array{Int64, 2})

    # @@ This routine may consume more memory than necessary if there is a large
    # spread in region numbers.
    
    (min_reg1, max_reg1) = extrema(regions)
    (min_reg2, max_reg2) = extrema(fromto[:,1])
    min_reg = min(min_reg1, min_reg2)
    max_reg = max(max_reg1, max_reg2)

    mapping = OffsetArray(collect(min_reg:max_reg), min_reg-1)
    mapping[fromto[:,1]] .= fromto[:,2]

    regions .= mapping[regions]
end

# ----------------------------------------------------------------------------
function _renumerate_regions!(regions::Matrix{Int64}; exitregions=[],
                              reindex=nothing)
    
    regnums = unique(regions)
    regnums = setdiff(regnums, [0]); # remove 0, which represents masked areas
    
    new_numbering = zeros(Int64, length(regnums))

    insidereg = .!in.(regnums, Ref(exitregions))
    outsidereg = .!insidereg
    new_numbering[insidereg] = 1:sum(insidereg)
    new_numbering[outsidereg] = -(1:sum(outsidereg))
    
    _remap!(regions, [regnums new_numbering])

    if reindex != nothing
        #local_region = view(regions, tiles[i].xrange, tiles[i].yrange)
        _remap!(view(reindex, :, 2),
                [regnums new_numbering])
    end
end


# ----------------------------------------------------------------------------
function _update_correspondences!(regions::Matrix{Int64},
                                  bcorr::Vector{Tuple{Int64, Int64}})

    
    # Note: the Graphs library doesn't support nodes with nonpositive indices.
    # Before calling this function, make sure all regions have positive indices.
    @assert(minimum(regions) > 0)

    if isempty(bcorr)
        return
    end

    components = _determine_connected_components(bcorr, length(regions))

    reg_ix = Vector{Int64}()
    
    for c in components
        ix = minimum(c)
        append!(reg_ix, fill(ix, length(c)))
    end
    mapping = hcat(vcat(components...), reg_ix)

    
    # give each grouping the region number corresponding to the lowest region
    # number of its members.
    _remap!(regions, mapping)
    
end


# ----------------------------------------------------------------------------
function _compute_region_correspondences!(correspondences::Vector{Tuple{Int64, Int64}},
                                          spillfield::Matrix{Int8}, 
                                          regions::Matrix{Int64},
                                          split::Int64, dir::Symbol)

    if split < 1
        # this does not divide the domain into two separate parts.  Nothing to do.
        return
    end

    insidereg = (dir == :x)  ? regions[split, :]   : regions[:, split]
    outsidereg = (dir == :x) ? regions[split-1, :] : regions[:, split-1]

    insideflowdir  = (dir == :x) ? spillfield[split, :]   : spillfield[:, split]
    outsideflowdir = (dir == :x) ? spillfield[split-1, :] : spillfield[:, split-1]

    icorr = (dir == :x) ? [4 0 7] : [4 2 6]
    ocorr = (dir == :x) ? [5 1 6] : [5 3 7]
    
    for shift in -1:1
        rstart = 1 - min(shift, 0)
        rend   = size(regions, dir == :x ? 2 : 1) - max(shift, 0)

        ic = icorr[shift + 2]
        oc = ocorr[shift + 2]
        
        for i in rstart:rend
            # check if one cell spills into the other, or alternatively, if both
            # cells are (adjacent) bottom cells
            if ((insideflowdir[i] == ic) || (outsideflowdir[i + shift] == oc)) ||
               ((insideflowdir[i] == -1) && (outsideflowdir[i + shift] == -1))
                push!(correspondences, (insidereg[i], outsidereg[i + shift]))
            end
        end
    end
end
# ----------------------------------------------------------------------------
function _ensure_unique_region_numbering!(regions, tiles)

    # ensuring unique global region numbers across tiles
    (min_reg, max_reg) = extrema(regions)
    @assert(min_reg > 0)
    offset = max_reg

    Threads.@threads for i = 1:prod(size(tiles))
        local_region = view(regions, tiles[i].xrange, tiles[i].yrange)
        local_region .+= (i-1) * offset
    end
    
end

# ----------------------------------------------------------------------------
function _fix_boundary_seams!(regions::Matrix{Int64}, spillfield::Matrix{Int8},
                              usediags::Bool, xsplits::Vector{Int}, ysplits::Vector{Int})

    # ensuring unique global region numbers across tiles    
    _ensure_unique_region_numbering!(regions,
                                     tiles([1, xsplits'..., size(regions, 1)+1],
                                           [1, ysplits'..., size(regions, 2)+1]))

    # map which regions corresponds across tiles
    correspondences = Vector{Tuple{Int64, Int64}}()

    for xs in xsplits
        _compute_region_correspondences!(correspondences, spillfield, regions, xs, :x)
    end
    
    for ys in ysplits
        _compute_region_correspondences!(correspondences, spillfield, regions, ys, :y)
    end
                
    # update region numbers based on the detected correspondences
    _update_correspondences!(regions, correspondences)


end

# ----------------------------------------------------------------------------
function _find_exit_nodes(spillfield)

    sz = size(spillfield)
    enodes = Vector{Int}()

    linear = LinearIndices(sz)

    # searching boundary j=1 for exit points
    for I in CartesianIndex(1,1):CartesianIndex(sz[1], 1)
        dir = spillfield[I]
        if dir == -1 || dir == 2 || dir == 4 || dir == 6
            push!(enodes, linear[I])
        end
    end

    # searching boundary j=jmax for exit points
    for I in CartesianIndex(1, sz[2]):CartesianIndex(sz)
        dir = spillfield[I]
        if dir == -1 || dir == 3 || dir == 5 || dir == 7
            push!(enodes, linear[I])
        end
    end

    # searching boundary i=imin for exit points
    for I in CartesianIndex(1,1):CartesianIndex(1, sz[2])
        dir = spillfield[I]
        if dir == -1 || dir == 0 || dir == 4 || dir == 7
            push!(enodes, linear[I])
        end
    end

    # searching boundary i=imax for exit points
    for I in CartesianIndex(sz[1], 1):CartesianIndex(sz)
        dir = spillfield[I]

        if dir == -1 || dir == 1 || dir == 5 || dir == 6
            push!(enodes, linear[I])
        end
    end

    # add sink nodes, if any
    sinkcells = linear[findall(spillfield .== -3)]
    enodes = [enodes; sinkcells]
    
    return enodes
end

# ----------------------------------------------------------------------------
# Identify all edges that constitutes "streamlines"
function _spillfield_flow_edges!(edges, spillfield)

    # * -3 : sink (any passing streamline is terminated here)
    # * -2 : covered by building / clipped away
    # * -1 : no downward slope (gridcell is a trap)
    # *  0 : steepest slope towards (i-1, j)
    # *  1 : steepest slope towards (i+1, j)
    # *  2 : steepest slope towards (i, j-1)
    # *  3 : steepest slope towards (i, j+1)
    # *  4 : steepest slope towards (i-1, j-1)
    # *  5 : steepest slope towards (i+1, j+1)
    # *  6 : steepest slope towards (i+1, j-1)
    # *  7 : steepest slope towards (i-1, j+1)

    linear = LinearIndices(size(spillfield)); # @@ does this consume memory?
    
    for cur_node in CartesianIndices(spillfield)

        dir = spillfield[cur_node]
        n1 = linear[cur_node]

        if dir == _SPILLFIELD_INACTIVE_FLAG
            continue
        elseif dir in [-3, -1]
            # this is a sink or the bottom of a trap, no downward flow from
            # here.  Add an edge of this cell pointing to itself, to ensure it
            # will be not ignored by the graph algorithm identifying components,
            # in the rare case that no other cell spills into this one.
            push!(edges, (n1, n1))
            continue
        elseif dir == -2 # masked gridcell (buildings, etc.).  No flow here
            continue
        end
        
        ishift = (dir == 0 || dir == 4 || dir == 7) ? -1 :
                 (dir == 1 || dir == 5 || dir == 6) ?  1 : 0
        jshift = (dir == 2 || dir == 4 || dir == 6) ? -1 :
                 (dir == 3 || dir == 5 || dir == 7) ?  1 : 0

        neigh_node = cur_node + CartesianIndex(ishift, jshift)

        if neigh_node[1] > 0 && neigh_node[1] <= size(spillfield, 1) &&
           neigh_node[2] > 0 && neigh_node[2] <= size(spillfield, 2)
            
            n2 = linear[neigh_node]
            push!(edges, (n1, n2))
        else
            # this node points out of the domain.  In case no other node
            # connects to it, add a self-connection, to ensure it will not be
            # ignored by the graph algorithm identifying components
            push!(edges, (n1, n1))
        end
    end
end

# ----------------------------------------------------------------------------
# Identify all grid edges that connects two neighbor trap bottoms.  If two trap
# bottoms are connected as neighbors, they have the same z-value, and are part
# of the same 'bottom point' of a trap.
function _flat_zone_connecting_edges!(edges, M, usediags::Bool=true)

    nodeixs = findall(M .== -1); # locate all cells marked as trap bottoms

    num_nodes = length(nodeixs)

    linear = LinearIndices(size(M))
    
    IMin = CartesianIndex(1, 1)
    IMax = CartesianIndex(size(M))

    shift = usediags ?
        [CartesianIndex(1, 0); CartesianIndex(-1, 1); CartesianIndex(0, 1); CartesianIndex(1,1)] :
        [CartesianIndex(1, 0); CartesianIndex(0, 1)]

    for ix_startnode in nodeixs

        for s in shift

            ix_endnode = ix_startnode + s

            if ix_endnode[1] <= IMax[1] && ix_endnode[2] <= IMax[2] &&
               ix_endnode[1] >= IMin[1] && ix_endnode[2] >= IMin[2]
                # we found a neighbor to this node, add the edge to edgelist
                if M[ix_endnode] == -1
                    n1 = linear[ix_startnode]
                    n2 = linear[ix_endnode]
                    push!(edges, (n1, n2))
                end
            end
        end
    end
end
