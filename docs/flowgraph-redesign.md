# FlowGraph Redesign: Compact Forward-Adjacency Flow Graph

## Problem

`TrapStructure.flowgraph` is currently a `Graphs.SimpleDiGraph`, which stores
both forward and backward adjacency as `Vector{Vector{Int64}}`.  For a terrain
with 80 million cells each vector has 80M entries, and each entry is a
heap-allocated `Vector{Int64}`.  Due to per-vector metadata overhead (~24 bytes
each), the two adjacency structures alone cost roughly:

    80M × 2 × 24 bytes ≈ 3.8 GB  (vector metadata)
    + 2 × 80M × 8 bytes  ≈ 1.3 GB  (actual edge data)

Total: **~5 GB** just for the flow graph.

---

## Constraints

Three properties of the terrain flow graph prevent simpler representations:

1. **Culverts** create non-local edges (pipe connections under roads that link
   non-adjacent cells).  Any replacement must support arbitrary edges, not just
   8-connectivity.
2. **Barriers** cut some connections between adjacent trap-bottom cells in flat
   zones (cells that appear identical in the spillfield).  These cuts are
   incorporated into the edge set at build time; they cannot be recovered from
   the spillfield alone.
3. **Future multi-direction flow**: some algorithms may eventually emit cells
   with more than one downstream neighbour.  The representation must not
   structurally prevent this.

**Rejected alternative — spillfield only:** `Matrix{Int8}` cannot represent
culverts (non-adjacent edges) or barrier-cut flat-zone pairs, and does not
support multi-direction flow.

---

## Chosen Solution: Custom `FlowGraph` Struct

```julia
mutable struct FlowGraph
    forward::Vector{Int32}            # forward[i]=j: cell i flows to j; 0=none
    multi::Dict{Int32, Vector{Int32}} # cells with >1 downstream (future / rare)
    topo_order::Vector{Int32}         # precomputed topological order (sources first)
    _backward::Union{Nothing, Vector{Vector{Int32}}}  # built lazily
end
```

### Why This Works

- **`forward` encodes all edges**, including culverts and barrier cuts, which are
  already baked into the edge set during `spillanalysis` before the graph is built.
  No runtime overrides are needed for the single-downstream case (the common one).
- **`multi` handles the general case** without touching the hot path.  Only cells
  with ≥2 downstream neighbours appear here (rare flat-zone configurations and
  future multi-direction).
- **`topo_order` cached once at construction** — `watercourses` previously
  re-ran `Graphs.topological_sort_by_dfs` (O(N)) on every invocation.  Now it
  reads a precomputed vector at zero cost.
- **`_backward` is built lazily** — backward adjacency is only used by
  `current_upstream_area`, a UI-only function triggered by a user click.  The
  O(N) build cost on first call is acceptable; it is not part of steady-state
  memory and can be GC'd if needed.

### Memory Comparison

| Representation | Memory (80M cells) |
|---|---|
| `SimpleDiGraph` (current) | ~5 GB |
| `FlowGraph` forward only (`Int32`) | ~320 MB |
| `FlowGraph` + backward (when built) | ~640 MB |

---

## Public Interface

```julia
# Downstream neighbours — returns [] for trap bottom / domain exit
downstream_cells(g::FlowGraph, cell::Int) → Vector{Int32}

# Fast boolean test (no allocation)
has_downstream(g::FlowGraph, cell::Int) → Bool

# Precomputed topological order (sources first, used by watercourses)
topological_order(g::FlowGraph) → Vector{Int32}

# All cells upstream of `start` via DFS on lazily-built backward adjacency
upstream_dfs(g::FlowGraph, start::Int) → Vector{Int}
```

---

## Construction

In `spillregions.jl`, replace:
```julia
g = Graphs.SimpleDiGraph([Graphs.SimpleEdge{Int64}((e[1], e[2])) for e in edges
                              if e[1] != e[2]])
Graphs.add_vertices!(g, N - Graphs.nv(g))
return regions, g, bottomcells
```
With:
```julia
fg = FlowGraph(edges, N)
return regions, fg, bottomcells
```

The `FlowGraph(edges, N)` constructor:
1. Iterates `edges` once to populate `forward` and `multi`
2. Runs Kahn's algorithm (O(N+E)) to build `topo_order`
3. Returns the fully initialised struct (no backward adjacency yet)

### Topological Sort: Kahn's Algorithm

Uses the result vector itself as the BFS queue (no extra allocation):
```julia
# 1. Compute in-degrees in one O(N) pass over forward
indegree[forward[i]] += 1   # for all i where forward[i] > 0

# 2. Seed queue with zero-indegree cells
push all i where indegree[i] == 0

# 3. BFS: pop front cell, decrement downstream's in-degree, enqueue if zero
```

Produces the same semantics as DFS-based topological sort.

---

## Call Sites Updated

| File | Old call | New call |
|---|---|---|
| `watercourses.jl` | `Graphs.topological_sort_by_dfs(g)` | `topological_order(g)` |
| `watercourses.jl` | `Graphs.outneighbors(g, cur_node)` | `downstream_cells(g, cur_node)` |
| `watercourses.jl` `_update_runoff!` | `Graphs.outneighbors(dstream, target)` | `downstream_cells(dstream, target)` |
| `watercourses.jl` `_trace_path` | `Graphs.outneighbors(tstruct.flowgraph, …)` | `downstream_cells(…)` |
| `fill_sequence/flow.jl` `_is_trap_bottom` | `isempty(Graphs.outneighbors(…))` | `!has_downstream(…)` |
| `fill_sequence/flow.jl` `_track_flow!` | `Graphs.outneighbors(tstruct.flowgraph, cell)` | `downstream_cells(…)` |
| `utils.jl` `_downstream_cell` | `Graphs.outneighbors(flowgraph, lix)` | `downstream_cells(flowgraph, lix)` |
| `utils.jl` `upstream_area` | `Graphs.dfs_parents(tstruct.flowgraph, …)` | `upstream_dfs(…)` |
| `utils.jl` `reconstruct_spillfield` | `Graphs.outneighbors(tstruct.flowgraph, …)` | `downstream_cells(…)` |
| `TrapStructure.jl` field | `flowgraph::Graphs.SimpleDiGraph` | `flowgraph::FlowGraph` |

---

## Files Added / Changed

- **New:** `src/FlowGraph.jl` — struct definition and all interface functions
- **Modified:** `src/SurfaceWaterIntegratedModeling.jl` — `include("FlowGraph.jl")` before `spillregions.jl`
- **Modified:** `src/TrapStructure.jl` — field type change
- **Modified:** `src/spillregions.jl` — graph construction
- **Modified:** `src/watercourses.jl` — topological sort + outneighbors calls + `_update_runoff!` signature
- **Modified:** `src/fill_sequence/flow.jl` — two call sites
- **Modified:** `src/utils.jl` — three call sites

---

## Notes / Future Work

- The backward adjacency (`_backward`) is built once on first `upstream_dfs`
  call and then cached.  If memory pressure is severe it could be cleared after
  use; the lazy build cost (O(N) ≈ seconds for 80M cells) is acceptable for a
  UI callback.
- `multi` is populated during construction for any cell that genuinely has more
  than one outgoing edge in the input edge set.  In current terrain data this
  occurs for some flat-zone trap-bottom cells.  Future multi-direction flow
  algorithms would populate it deliberately.
- `Spillpoint.elevation::Real` (abstract type) is a separate, smaller
  inefficiency (type instability + boxing).  Parameterising it as
  `Spillpoint{T<:Real}` is straightforward and documented in
  `docs/performance_review.md`.
