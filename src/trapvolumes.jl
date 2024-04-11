export trapvolumes


"""
    trapvolumes(grid, spillregions, spillpoints, lowest_regions)

Compute the volume of all traps (at all levels) identified for a topography grid.

The necessary input to this function can be provided by [`spillregions`](@ref), 
[`spillpoints`](@ref) and [`sshierarchy!`](@ref).  The result is returned in form
of a vector of floating point numbers, with one entry per trap, giving the volume 
of that trap.

# Arguments
    - `grid::Matrix{<:Real}`: the topographical grid that has been analysed 
    - `spillregions::Matrix{Int}`: of same size as `grid`, where the integer at entry
                                   (i, j) gives the region number of cell (i, j) in `grid`.
    - `spillpoints::Vector{Spillpoint}`: vector with the identified spillpoints, as
                                         provided by the [`spillpoints`](@ref) function.
                                         There is one spillpoint per trap.
    - `lowest_regions::Vector{Vector{Int}}`: each entry in this vector corresponds to 
                                             a trap, and list all the low-level regions
                                             that combined constitute the spill region
                                             of this trap.

See also [`spillregions`](@ref), [`spillpoints`](@ref) and [`sshierarchy!`](@ref).  
"""
function trapvolumes(grid::Matrix{<:Real},
                     spillregions::Matrix{Int},
                     spillpoints::Vector{Spillpoint},
                     lowest_regions::Vector{Vector{Int}})

    num_lowlevel_regions = maximum(spillregions)
    num_all_regions = length(lowest_regions)
    @assert(length(spillpoints) == num_all_regions)
    if num_all_regions == 0
        return Vector{Int64}()
    end

    # for each lowest-level spill region, determine all supertraps by inverting
    # 'lowest_regions'
    supertraps = [[r] for r in 1:num_lowlevel_regions]

    for superreg_ix in num_lowlevel_regions+1:num_all_regions
        for reg in lowest_regions[superreg_ix]
            push!(supertraps[reg], superreg_ix)
        end
    end

    elevations = [sp.elevation for sp in spillpoints]

    num_threads = 1 #@@Threads.nthreads()
    chunksize = Int(ceil(length(spillregions)/num_threads))

    function process_chunk(first, last)
        resvec = fill(0.0, num_all_regions)
        for ix in first:last
            lowreg = spillregions[ix]
            zval = grid[ix]
            if lowreg > 0
                allregs = supertraps[lowreg]
                if length(allregs) == 1
                    reg = allregs[1]; # only a single region, so no need to build vector
                    dz = elevations[reg] - zval
                    if dz > 0
                        resvec[reg] += dz
                    end
                else
                    dzvals = [max(elevations[r] - zval, 0) for r in allregs]
                    resvec[allregs] .+= dzvals
                end
            end
        end
        return resvec
    end

    # We work on separate accumulating vectors to prevent locking issues.  In
    # the end, we sum up all the partial results go get the real answer.
    res = Vector{Any}(undef, num_threads)
    for i in 1:num_threads
        res[i] = Threads.@spawn process_chunk((i-1)*chunksize+1,
                                              min(i*chunksize, length(spillregions)))
    end
    volumes = fetch(res[1])
    for i in 2:num_threads
        volumes += fetch(res[i])
    end

    return volumes
end

                     
