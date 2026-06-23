# Performance Review

This document identifies inefficiencies in the codebase that may matter when processing
large terrain grids.  Issues are grouped by estimated impact.  No code has been changed;
this is a catalogue for prioritisation.

---

## High Priority

### 1. `Spillpoint.elevation` is an abstract type — pervasive type instability

**File:** `src/spillpoints.jl:22`

```julia
struct Spillpoint
    ...
    elevation::Real   # ← abstract
end
```

`Real` is an abstract type, so the struct is not fully concrete.  Every access to
`spillpoint.elevation` anywhere in the codebase triggers dynamic dispatch and heap
boxing.  Spillpoints are accessed in essentially every hot path (volume computation,
fill tracking, slope comparisons).

**Suggestion:** Parameterise the struct:

```julia
struct Spillpoint{T<:Real}
    downstream_region::Int
    current_region_cell::Int
    downstream_region_cell::Int
    elevation::T
end
```

---

### 2. `waterheight` allocates a large array inside the root-finder loop

**File:** `src/utils.jl:288`

```julia
volfun = z -> sum(max.(z .- tstruct.topography[trapcells], 0.0))
```

`tstruct.topography[trapcells]` performs an index-gather and allocates a new array
*on every evaluation of the closure*.  The root-finder (`Roots.find_zero`) calls this
closure many times per trap.

**Suggestion:** Capture the slice once before the closure:

```julia
bottom = tstruct.topography[trapcells]
volfun = z -> sum(max.(z .- bottom, 0.0))
```

---

### 3. `_setup_dvdt` allocates a `BitArray` on every ODE right-hand-side evaluation

**File:** `src/fill_sequence/fill_sequence.jl:524`

```julia
return dv[1] = inflow - infilfun(trap_bottom .<= z)
```

`trap_bottom .<= z` produces a fresh `BitArray` on each call.  This function is the
inner loop of an ODE solver and is evaluated hundreds of times per fill event.

**Suggestion:** Pre-sort `trap_bottom` (it is already sorted in `_compute_z_vol_tables`)
and replace the broadcast with a binary-search-based area count, which requires no
allocation and is also faster algorithmically.

---

### 4. Intermediate `Vector{SimpleEdge}` allocated for every graph construction

**Files:** `src/spillregions.jl:107`, `src/spillregions.jl:250`,
`src/utils.jl:348`, `src/utils.jl:419`, `src/utils.jl:1372`

Repeated pattern:

```julia
Graphs.SimpleDiGraph([Graphs.SimpleEdge{Int64}((e[1], e[2])) for e in edges
                          if e[1] != e[2]])
```

This materialises a full intermediate `Vector{SimpleEdge}` before handing it to the
graph constructor.  For millions of terrain cells this doubles peak memory and
allocation time.

**Suggestion:** Use `add_edge!` in a loop, or pass an edge iterator directly if the
Graphs version supports it.  At minimum, the filter `if e[1] != e[2]` should be
combined with the construction to avoid a two-pass scan.

---

### 5. `_is_sink` reconstructs `CartesianIndices` and does a linear search on every call

**File:** `src/fill_sequence/flow.jl:96–99`

```julia
function _is_sink(cell, tstruct)
    return cell > 0 && (tstruct.sinks !== nothing) &&
        CartesianIndices(size(tstruct.topography))[cell] in tstruct.sinks
end
```

Two problems, each compounding because this is called inside `_track_flow!` which is
itself in a loop:

- `CartesianIndices(size(tstruct.topography))` is reconstructed on every call.
- `in tstruct.sinks` is a `Vector` linear search — O(number of sinks) per call.

**Suggestion:** Store a pre-computed `Set{Int}` of linear sink indices in
`TrapStructure` at construction time, and query that instead.

---

## Medium Priority

### 6. `trapvolumes` allocates a vector for every multi-region cell

**File:** `src/trapvolumes.jl:69`

```julia
dzvals = [max(elevations[r] - zval, 0) for r in allregs]
resvec[allregs] .+= dzvals
```

A new `dzvals` vector is allocated for each grid cell that belongs to more than one
supertrap.  For large nested traps this is a significant number of allocations.

**Suggestion:** Replace with a scalar accumulation loop:

```julia
for r in allregs
    resvec[r] += max(elevations[r] - zval, 0)
end
```

---

### 7. `_prepare_cut_edges_bidir` recomputes `LinearIndices` inside the inner loop

**File:** `src/spillregions.jl:400–401`

```julia
push!(cut_edges_bidir, (LinearIndices(gridsize)[k],
                        LinearIndices(gridsize)[neigh]))
```

`LinearIndices(gridsize)` is constructed twice per neighbor per key.  It should be
computed once outside both loops.

---

### 8. Assertion in the weather-event loop materialises a `Bool` vector

**File:** `src/fill_sequence/fill_sequence.jl:65`

```julia
@assert(all([ca.time == cur_time for ca ∈ cur_amounts]))
```

The `[...]` forces allocation of a full boolean vector before `all` inspects it.
Use the generator form which short-circuits and allocates nothing:

```julia
@assert all(ca.time == cur_time for ca ∈ cur_amounts)
```

---

### 9. `_identify_next_status_change!` — `findall` over a materialised comprehension

**File:** `src/fill_sequence/fill_sequence.jl:283–285`

```julia
candidates = findall([all(filled_traps[subtrapsof(tstruct, ix)]) &&
                      changetimeest[ix].min < earliest_changetime
                      for ix ∈ 1:num_traps])
```

Materialises a `Bool` vector for the sole purpose of passing it to `findall`.
Use the predicate form instead:

```julia
candidates = findall(ix -> all(filled_traps[subtrapsof(tstruct, ix)]) &&
                           changetimeest[ix].min < earliest_changetime,
                     1:num_traps)
```

---

### 10. `_compute_initial_rateinfo` uses `ones(...)` multiplication instead of `fill`

**File:** `src/fill_sequence/flow.jl:174–175`

```julia
precipitation = precipitation .* ones(size(tstruct.regions))
```

Allocates a temporary ones-array and then multiplies element-wise.  The equivalent
`fill(precipitation, size(tstruct.regions))` is both simpler and avoids the
intermediate allocation.

---

### 11. `Set{Tuple{Int,Int}}` used for the large flow-edges collection

**Files:** `src/spillregions.jl:57`, `src/spillregions.jl:282`

The `edges` collection is a `Set`, paying hashing overhead on every `push!` for
what is often millions of edges.  The only duplicates that need to be suppressed are
self-loops, and those are rare.

**Suggestion:** Use a `Vector` with `sizehint!` based on the domain size, and filter
self-loops at graph-construction time.  This is already done correctly in
`compute_spillfield_graph` (`utils.jl:1040`); the same pattern could be applied in
`_process_domain!`.

---

### 12. `_find_trapcells` scans the full region matrix per region

**File:** `src/utils.jl:189`

```julia
regcells = findall(tstruct.regions[:] .== r)
```

`tstruct.regions[:]` creates a copy of the entire region matrix, and this is done
once per region inside the loop — O(N × number of regions) total work.  The
`_regionmap` helper in the same file already builds the inverse map efficiently in
O(N); `_find_trapcells` should use it.

---

## Lower Priority / Structural

### 13. `_active_region_spilltree` runs a full DFS once per filled trap

**File:** `src/utils.jl:1362`

```julia
subtraps = findall(Graphs.dfs_parents(tstruct.agglomerations, i, dir=:in) .> 0)
```

`dfs_parents` allocates a full-size array and performs a complete DFS for every
filled trap in the loop.  Total cost is O(T × (V + E)) where T is the number of
filled traps.

**Suggestion:** A single multi-source BFS from all filled traps simultaneously, or
pre-computed subtrap ancestry sets, reduces this to O(V + E) total.

---

### 14. `all_subtraps_of` mixes `Set` and `Vector` for `active_set`, allocating each BFS level

**File:** `src/utils.jl:481–490`

```julia
active_set = Set(trap_ixs)       # first iteration: Set
while !isempty(active_set)
    found = Vector{Int64}()       # allocates every level
    for i in active_set
        union!(found, Graphs.inneighbors(subtrap_graph, i))
    end
    union!(result, found)
    active_set = found            # subsequent iterations: Vector (no dedup)
end
```

`found` is a new `Vector` each BFS level, and after the first level `active_set` is
a `Vector` so nodes are not deduplicated before the next round.  Reusing two
pre-allocated sets and swapping them would eliminate allocations and ensure
correctness.

---

### 15. Multithreaded path in `trapvolumes` is disabled

**File:** `src/trapvolumes.jl:52`

```julia
num_threads = 1 #@@@Threads.nthreads()
```

The multithreaded chunked accumulation path is present and correctly avoids race
conditions (each thread writes to its own `resvec`), but has been commented out.
`trapvolumes` iterates every grid cell and is one of the most parallelisable
bottlenecks in the pipeline.  Re-enabling this with `Threads.nthreads()` would give
a near-linear speedup on multi-core machines.

---

### 16. `@inbounds` missing in inner loops with provably safe indexing

Several loops iterate `CartesianIndices` or linear ranges whose bounds are
structurally guaranteed, but lack `@inbounds`:

- `src/spillfield.jl:455–465` — `_compare_slopes` inner loop; indices are derived
  from `shape` which was just computed from the array size.
- `src/spillregions.jl:623–650` — `_spillfield_flow_edges!`; all indices are bounds-
  checked once by the `if` guards before use.
- `src/trapvolumes.jl:57–74` — the per-cell processing loop; `ix` is always a valid
  linear index of `spillregions`.

Adding `@inbounds` (with care) to these loops avoids redundant bounds checks on
every array access in what are essentially O(N) scans of the full terrain.
