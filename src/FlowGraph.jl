export FlowGraph, downstream_cells, has_downstream, topological_order, upstream_dfs

"""
    mutable struct FlowGraph

Compact forward-adjacency flow graph replacing `Graphs.SimpleDiGraph` as the
`flowgraph` field of `TrapStructure`.

Stores downstream adjacency as a flat `Vector{Int32}` (`forward[i] = j` means
cell i flows to cell j; 0 means no downstream neighbour).  Cells with more than
one downstream neighbour (rare flat-zone configurations, future multi-direction
flow) are recorded in `multi`.

The topological order is precomputed once at construction time via Kahn's
algorithm, so `watercourses` no longer needs to re-sort on every call.

The backward adjacency (`_backward`) is built lazily on the first call to
`upstream_dfs`; it is used only by the UI-side `current_upstream_area` function.

Memory for 80M cells: ~320 MB (vs ~5 GB for `SimpleDiGraph`).
See `docs/flowgraph-redesign.md` for the full design rationale.
"""
mutable struct FlowGraph
    forward::Vector{Int32}              # forward[i]=j: cell i → j; 0=no downstream
    multi::Dict{Int32, Vector{Int32}}   # cells with >1 downstream neighbour
    topo_order::Vector{Int32}           # precomputed topological order, sources first
    _backward::Union{Nothing, Vector{Vector{Int32}}}  # lazily built reverse adjacency
end

"""
    FlowGraph(edges, n_cells)

Build a `FlowGraph` from an iterable of `(i, j)` edge pairs.

Self-loops (`i == j`) are skipped.  Cells with multiple outgoing edges in
`edges` are stored in `multi` with `forward[i] = 0` as a sentinel.

`n_cells` must equal the total number of grid cells (vertices).
"""
function FlowGraph(edges, n_cells::Int)
    forward = zeros(Int32, n_cells)
    multi   = Dict{Int32, Vector{Int32}}()

    for (i, j) in edges
        i == j && continue
        i32, j32 = Int32(i), Int32(j)
        if haskey(multi, i32)
            push!(multi[i32], j32)
        elseif forward[i] != 0
            # second outgoing edge for this cell — promote to multi
            multi[i32] = [forward[i], j32]
            forward[i] = Int32(0)
        else
            forward[i] = j32
        end
    end

    topo = _compute_topo_order(forward, multi, n_cells)
    return FlowGraph(forward, multi, topo, nothing)
end

# ----------------------------------------------------------------------------
"""
    downstream_cells(g, cell) → Vector{Int32}

Return the downstream neighbour(s) of `cell`.  Returns an empty vector for
trap bottoms and domain exits.  Almost always returns a single-element vector.
"""
@inline function downstream_cells(g::FlowGraph, cell::Int)
    i32 = Int32(cell)
    if haskey(g.multi, i32)
        return g.multi[i32]
    end
    j = g.forward[cell]
    return j == 0 ? Int32[] : Int32[j]
end

# ----------------------------------------------------------------------------
"""
    has_downstream(g, cell) → Bool

Return `true` if `cell` has at least one downstream neighbour (i.e. is not a
trap bottom or domain exit).
"""
@inline function has_downstream(g::FlowGraph, cell::Int)
    return g.forward[cell] != 0 || haskey(g.multi, Int32(cell))
end

# ----------------------------------------------------------------------------
"""
    topological_order(g) → Vector{Int32}

Return the precomputed topological order of all cells (sources before sinks).
"""
@inline topological_order(g::FlowGraph) = g.topo_order

# ----------------------------------------------------------------------------
"""
    upstream_dfs(g, start) → Vector{Int}

Return all cells that drain (directly or indirectly) into `start`, found via
DFS on the backward adjacency.  The backward adjacency is built lazily on the
first call.

`start` itself is excluded from the result, matching the behaviour of
`Graphs.dfs_parents(g, start, dir=:in) .> 0`.
"""
function upstream_dfs(g::FlowGraph, start::Int)
    isnothing(g._backward) && _build_backward!(g)
    n       = length(g.forward)
    visited = falses(n)
    result  = Int[]
    stack   = Int[start]
    visited[start] = true
    while !isempty(stack)
        cell = pop!(stack)
        cell != start && push!(result, cell)
        for up in g._backward[cell]
            up_int = Int(up)
            if !visited[up_int]
                visited[up_int] = true
                push!(stack, up_int)
            end
        end
    end
    return result
end

# ----------------------------------------------------------------------------
# Internal helpers

function _compute_topo_order(forward::Vector{Int32},
                             multi::Dict{Int32, Vector{Int32}},
                             n::Int)
    # Kahn's algorithm: BFS from zero-in-degree nodes.
    # The result vector doubles as the queue (head pointer avoids O(N) dequeue).
    indegree = zeros(Int32, n)
    for i in eachindex(forward)
        j = forward[i]
        j > 0 && (indegree[j] += 1)
    end
    for targets in values(multi)
        for j in targets
            0 < j <= n && (indegree[j] += 1)
        end
    end

    topo = Int32[]
    sizehint!(topo, n)
    for i in 1:n
        indegree[i] == 0 && push!(topo, Int32(i))
    end

    head = 1
    while head <= length(topo)
        cell = topo[head]; head += 1
        i32  = Int32(cell)
        if haskey(multi, i32)
            for j in multi[i32]
                0 < j <= n || continue
                indegree[j] -= 1
                indegree[j] == 0 && push!(topo, j)
            end
        else
            j = forward[cell]
            if j > 0
                indegree[j] -= 1
                indegree[j] == 0 && push!(topo, j)
            end
        end
    end

    # Graph should be a DAG for terrain flow; append any cycle members last
    if length(topo) < n
        in_topo = falses(n)
        @inbounds in_topo[topo] .= true
        for i in 1:n
            in_topo[i] || push!(topo, Int32(i))
        end
    end

    return topo
end

function _build_backward!(g::FlowGraph)
    n = length(g.forward)
    backward = [Int32[] for _ in 1:n]
    for i in eachindex(g.forward)
        j = g.forward[i]
        j > 0 && push!(backward[j], Int32(i))
    end
    for (i, targets) in g.multi
        for j in targets
            push!(backward[j], i)
        end
    end
    g._backward = backward
end
