# Int32 Migration: `regions` and `footprints`

## Motivation

For an 80M-cell terrain grid, `TrapStructure` holds two large `Int64` arrays:

- `regions::Matrix{Int64}` — one value per cell: **640 MB**
- `footprints::Vector{Vector{Int64}}` — collectively span all N cells: **640 MB**

Changing both to `Int32` halves these to 320 MB each, saving ~640 MB total.
No overflow risk: region count and linear cell indices are both bounded by N, and
80M << 2³¹ (~2.1B).

See also `flowgraph-redesign.md` (flowgraph: 5 GB → 320 MB) and
`performance_review.md` (catalogue of other inefficiencies).

---

## `regions::Matrix{Int32}`

### Hardcoded `Int64` signatures — mechanical updates required

`spillregions.jl` hard-codes `Int64` in approximately 8 places:

| Location | Change needed |
|---|---|
| `similar(spillfield, Int64)` (line 56) | `→ Int32` |
| `update_spillregions!(regions::Matrix{Int64}, ...)` | `→ Int32` |
| `_process_extended_domain!(regions::Matrix{Int64}, ...)` | `→ Int32` |
| `_process_domain!(regions::Matrix{Int64}, ...)` | `→ Int32` |
| `_remap!(regions::AbstractArray{Int64}, fromto::Array{Int64, 2})` | both `→ Int32` |
| `_renumerate_regions!(regions::Matrix{Int64}; ...)` | `→ Int32` |
| `new_numbering = zeros(Int64, ...)` inside `_renumerate_regions!` | `→ Int32` |
| `_update_correspondences!(regions::Matrix{Int64}, ...)` | `→ Int32` |
| `_fix_boundary_seams!(regions::Matrix{Int64}, ...)` | `→ Int32` |

All are mechanical — no algorithmic changes.

### `_remap!` and its `fromto` argument

`_renumerate_regions!` builds `new_numbering = zeros(Int64, ...)` and passes it as
`fromto` to `_remap!`.  Both must be updated together; the `fromto` array type must
match the `regions` element type.

### `OffsetArray` indexing — worth testing

In `_process_extended_domain!`:

```julia
(min_reg, max_reg) = extrema(regions)   # becomes Int32
mask = OffsetArray(fill(false, max_reg - min_reg + 1), min_reg-1)
mask[boundary_regs] .= true
```

`min_reg - 1` with `Int32` is fine.  `boundary_regs` come from `unique(regions[...])`
so also `Int32`.  OffsetArrays handles integer offsets generically, but verify in tests.

### `maximum(tstruct.regions)` as vector length

```julia
region_accum = zeros(Float64, max(maximum(tstruct.regions), 0))
```

`maximum` returns `Int32`; `max(Int32, 0)` promotes to `Int64` (literal `0` is
`Int64`), so `zeros(Float64, Int64)` — no change needed here.

---

## `footprints::Vector{Vector{Int32}}`

### **Critical:** `_compute_trap_footprints` pushes `Int64` indices into the vector

```julia
# spillanalysis.jl:218
footprints = [Vector{Int64}() for i in 1:num_traps]
...
for i in LinearIndices(regions)   # i is Int64 on 64-bit systems
    ...
    push!(footprints[tr], i)      # would throw InexactError with Vector{Int32}
```

`LinearIndices` returns `Int64` on 64-bit systems.  Pushing into a `Vector{Int32}`
without an explicit `Int32(i)` conversion throws `InexactError` at runtime.
**This is the one non-obvious correctness risk in the migration.**

Fix: change allocation to `Vector{Int32}()` and push site to `push!(footprints[tr], Int32(i))`.

### `setdiff` mixes `Int64` paths with `Int32` footprints

In `watercourses.jl:201`:

```julia
paths[end] = setdiff(paths[end], tstruct.footprints[largest_full_supertrap])
```

`paths[end]` is `Vector{Int64}` (path cells come from `Int64` linear indices).
`setdiff(Vector{Int64}, Vector{Int32})` works — Julia promotes for comparison and
returns `Vector{Int64}` — but is implicit type mixing.  Consider converting at the
call site: `Int.(tstruct.footprints[...])`, or making path vectors `Int32` too.

### Array indexing — no issue

All patterns like `runoff[tstruct.footprints[trap]]`,
`topography[tstruct.footprints[tix]]` work fine with `Int32` indices.  Julia
zero-extends `Int32` to `Int` for array subscripting; this is one instruction,
inlined and eliminated by the JIT.

---

## Runtime performance

| Aspect | Impact |
|---|---|
| Cache efficiency | Positive — smaller elements, more fit per cache line |
| SIMD scans | Positive — twice as many Int32 per register |
| Array subscript conversion `Int32 → Int` | ~1 instruction per access, negligible |
| OffsetArray / promotion edge cases | Test coverage recommended |

Expected outcome: neutral to slightly positive.

---

## Summary table

| Concern | Severity | Notes |
|---|---|---|
| ~9 hardcoded `Int64` signatures in `spillregions.jl` | Low — mechanical | All in internal helpers |
| `_compute_trap_footprints` pushes `Int64` into `Vector{Int32}` | **High — correctness** | Must add `Int32(i)` at push site |
| `setdiff` mixing `Int64` paths with `Int32` footprints | Low — works, untidy | Type mismatch across call boundary |
| OffsetArray indexing with `Int32` | Low — likely fine | Test to confirm |
| Overflow at any realistic terrain size | None | 80M << 2³¹ |
| Runtime regression | None expected | Likely slight improvement |
