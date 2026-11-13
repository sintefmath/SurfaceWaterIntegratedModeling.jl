# SWIM performance, memory and correctness review

- **Date:** 2026-11-11/12
- **Branch reviewed:** `SUrbArea` at `96f0b8b`
- **Scope:** whole `src/`. The focus is memory and CPU cost of the static spill-point analysis (`spillanalysis`) and of dynamic time sequencing (`fill_sequence` plus `dynamics/`) on large grids.
- **Status:** review only. **No code was changed.** Nothing was profiled; there is no Julia in the review sandbox. The cost estimates come from reading the code. Before starting any item, measure it (see §7).

**Notation:**
- N = grid cells; T = traps (all levels); R = lowest-level regions
- F = footprint size of one trap (up to ≈N for top-level traps)
- E = spill events; W = weather periods; P = sample timepoints; depth = nesting depth of the sub/supertrap hierarchy

**Verification marks:**
- **[V]** I re-read the code myself and confirmed the issue.
- **[R]** Found in a careful read by a review pass but not re-checked line by line. Re-read before acting.
- **[P]** Plausible, but needs a failing test to confirm.

Line numbers refer to commit `96f0b8b` and will drift. Search by function name if they no longer match.

---

## 0. Executive summary

| # | Item | Kind | Impact | Effort |
|---|------|------|--------|--------|
| D1 | Network ODE stores the full dense trajectory; only the endpoint is used | memory | very high (network runs) | trivial |
| F1 | Full runoff and rain grid copied into a `SpillEvent` every weather period | memory | very high (many periods) | medium |
| S1 | `flowgraph::SimpleDiGraph` over every cell (~130–150 B/cell) | memory | very high | large (layout change) |
| S2 | `spillregions` construction peak ~250–300 B/cell (edge Sets plus a graph) | memory | very high (peak) | large |
| F2 | `_track_flow!` → `_update_Smin_Smax!` rescans whole ancestor footprints | CPU | very high | medium |
| F3 | `fill_trap_until`: O(F) allocations, plus an ODE where a closed form exists | CPU | high | medium |
| D2 | Every network component committed and re-predicted twice per touched event | CPU | high | medium |
| D3 | Network RHS: O(F) per trap per call, about 8 allocations per call, untyped field | CPU | high | medium |
| S3 | `footprints` / `z_vol_tables` grow as O(N·depth) | memory | high (deep hierarchies) | large |
| S4 | Waterbody flattening is O(N·#lakes) | CPU | high when lakes present | trivial |
| S5 | `trapvolumes` makes two heap allocations per cell | CPU/GC | high | trivial |
| C1–C7 | Correctness bugs in the dynamic path (time origin, mass loss) | correctness | — | small–medium |
| C14 | Network commit used post-event runoff (NBS mass loss) — **fixed 2026-11-12** | correctness | high (NBS runs) | done |
| C8–C12 | Correctness bugs in static path / API | correctness | — | small |

Recommended order: see §6. **NBS track (§6.0) comes first** (user decision 2026-11-12).

---

## 1. Static analysis (`spillanalysis` and its callees)

### 1.1 What `TrapStructure` keeps resident

Constructed in `spillanalysis.jl:56-163`. Per-cell cost of each field:

| Field | Type | ≈ B/cell | Note |
|---|---|---|---|
| `topography` | `Matrix{T}` | 8 | the copy `gridcpy` (l.74) |
| `flowgraph` | `SimpleDiGraph{Int64}` | **130–150** | a `Vector{Int}` for forward and one for backward adjacency per vertex |
| `regions` | `Matrix{Int64}` | 8 | Int32 is enough while N < 2.1e9 |
| `footprints` | `Vector{Vector{Int}}` | 8 × depth per wet cell | see S3 |
| `clip_mask` | `Matrix{Bool}` | 1 | a `BitMatrix` input is converted to `Matrix{Bool}` |
| `waterbodies` | `Set{CartesianIndex{2}}` | 25–50 per water cell | only used for membership tests |
| `lowest_subtraps_for`, `supertraps_of` | `Vector{Vector{Int}}` | R × depth each | the same relation stored twice |

Total is about 160–180 B/cell plus footprints. A lean layout gets to about 13–20 B/cell: topography, Int32 regions, an Int32 downstream array, and compact footprints.

> **Constraint (user memory):** do not add *new* `TrapStructure` fields. Derived caches belong in `fill_sequence`. *Replacing* `flowgraph` with a more compact representation changes the layout of an existing field. It is breaking for users who access `tstruct.flowgraph`, so discuss it with the maintainer before doing it.

### S1. `flowgraph` as a full `SimpleDiGraph` [V]

- **Where:**
  - built in `spillregions.jl:104-109`;
  - field declared at `TrapStructure.jl:39`;
  - used in `watercourses.jl:60,128`, `utils.jl:911-937,1032,1280,1374,1765`, `fill_sequence/flow.jl:96,163`, and `dynamics/build_network.jl:85-109,398-447`.
- **Problem:**
  - Every cell has out-degree 0 or 1, because culverts are handled through region merging, not as extra out-edges. Confirm this in the `spillregions` culvert code before relying on it; `flow.jl:166` and `watercourses.jl:70` assert `length(ds)==1`.
  - The graph still stores forward and backward lists per vertex.
  - Building it also allocates a 16 B/cell edge vector plus degree arrays.
- **Type instability:** the field is declared as `Graphs.SimpleDiGraph` (a UnionAll), and `agglomerations` likewise (l.52). Every `outneighbors` / `inneighbors` call on it from a hot loop is dynamically dispatched.
- **Fix, quick:** declare `flowgraph::Graphs.SimpleDiGraph{Int}` and `agglomerations::Graphs.SimpleDiGraph{Int}`. This keeps the layout and only removes dispatch overhead.
- **Fix, full:** replace the graph with `downstream::Vector{Int32}`, where 0 means terminal (trap bottom, sink, or exit).
  - Callers needing upstream sets currently use `Graphs.dfs_parents(..., dir=:in)` (`utils.jl:1280,1765`) and `_map_inverse_flow`. Give them a CSR inverse built on demand: Int32 `offsets[N+1]` plus `sources[N]`, about 8 B/cell, freed after use.
  - The topological order `watercourses` needs comes from Kahn's algorithm on in-degree counts.
  - Provide accessor helpers (`_downstream(tstruct, c)`) so callers don't touch the field directly.
- **Risk:** public field change. Tests compare checksums and must stay bit-identical.
- **Prior work: branch `feature/compact-flowgraph`** (2026-06-06, pushed to origin, **not merged**; 4 commits ahead of `SUrbArea`).
  - `f2078ba` replaces the `SimpleDiGraph` flowgraph with a compact `FlowGraph` struct.
  - `350bd88` makes the `FlowGraph` API accept any Integer.
  - `1209c3e` fixes test destructuring for the `spillregions()` tuple return.
  - `c50f326` adds an Int32 migration analysis for `regions` and footprints, which is relevant to S3 and the Int32 item.
  - It predates all the July dynamics/NBS work (`dynamics/` uses `tstruct.flowgraph` in `build_network.jl` and `utils.jl:_downstream_cell`), so expect a non-trivial rebase.
  - **Start S1 from this branch, not from scratch.** Read its migration analysis before planning S1, S2 or S3.

### S2. Peak memory inside `spillregions` [R]

- **Where:**
  - `spillregions.jl:57` builds `edges::Set{Tuple{Int,Int}}` with ~1 entry per cell (25–51 B/cell at load factor ≤ 2/3).
  - `:327` does `union!(all_connections, edges)`, a second cell-sized Set.
  - `_determine_connected_components` (`:230-261`) allocates `findall(actives)`, `map = zeros(Int, max)`, an edge-vector comprehension, an undirected `SimpleGraph` (~70–90 B/cell), and `connected_components` output.
  - `:311` hashes all N edges even when there are no culverts.
  - `:370-371` creates three full-grid temporaries for the remap.
  - `findall(spilldomain .== -1)` is computed three times (`:374,403,697`).
- **Fix:** a single union-find over `parent::Vector{Int32}`:
  1. Union each cell with its downstream cell from `dir`.
  2. Union adjacent flat-bottom cells, checking cut edges with a direct `get(cut_edges, ...)` lookup instead of Set difference.
  3. Apply culvert unions.
  4. Compact the labels to 1..R in one pass.

  This gives `regions` (Int32), `downstream` (see S1) and `trap_bottoms` together, in O(N·α) time and about 8 B/cell of temporaries.
- **Related, quadratic leak-trap remap** (`spillregions.jl:352`): `remap[remap .== src_reg] .= ...` is O(R) per leak edge, so O(L·R) overall. Union-find removes it.
- **Related, renumbering** (`_renumerate_regions!`, `:452-472`):
  - `unique(regions)` hashes all N cells;
  - `in.(regnums, Ref(exitregions))` scans a Vector;
  - `_remap!` (`:435`) builds a full temporary from `mapping[regions]`.

  Use in-place loops and a mask.

### S3. Footprints and z/volume tables scale as O(N·depth) [V]

- **Where:**
  - `_compute_trap_footprints` (`spillanalysis.jl:215-249`) pushes a cell into every trap level whose spill elevation is ≥ the cell's z.
  - `_compute_z_vol_tables` (`fill_sequence.jl:676-701`) builds two Float64 vectors of length F per trap, i.e. 16 B × N × depth, plus sort temporaries.
  - `lowest_regions` (`sshierarchy.jl:82`) stores, for every trap, a full copy of its leaf regions: R × depth, O(R²) for chain-like hierarchies. `supertraps_of` (`spillanalysis.jl:179`) is the inverse of the same relation.
  - The same inverse is built three times: `_compute_supertraps_of`, `supertraps_for` (`spillanalysis.jl:223`), and `supertraps` (`trapvolumes.jl:42`).
- **Fix idea:**
  - Store per-region cell lists sorted by elevation (Int32).
  - The footprint of trap t is the concatenation over its leaf regions r of `cells_r[1:searchsortedlast(z_r, elev(t))]`. Expose it through an accessor or iterator; avoid materialising it.
  - The z/volume table for t can be derived by merging those sorted prefixes, or computed lazily per trap and cached only for traps actually used by `fill_sequence`.
  - Leaf lists can be represented as contiguous ranges, by ordering leaves in DFS order of the hierarchy (two Int32 per trap).
- **Caveat:** `footprints` is used in 50+ places (`grep -n footprints src`). Introduce an accessor first, migrate callers, then change the storage.

### S4. Waterbody flattening is O(N × #lakes) [V]

- **Where:** `_flatten_waterbody_regions!`, `spillanalysis.jl:166-175`. It runs `findall(wbody_labels .== label)` per label, and each call creates a full-grid BitArray.
- **Fix:** two passes. First compute `minz[label] = min(...)` over all cells, then assign `grid[c] = minz[label[c]]`.
- **Also:** `waterbodies` is stored as a `Set{CartesianIndex}` (`spillanalysis.jl:146`). Its only uses are membership tests (`watercourses.jl:134`, `utils.jl:942`); a `BitMatrix` is cheaper. This is a field-type change; see the constraint in §1.1.

### S5. `trapvolumes` allocates per cell [V]

- **Where:** `trapvolumes.jl:69-70`.
  ```julia
  dzvals = [max(elevations[r] - zval, 0) for r in allregs]
  resvec[allregs] .+= dzvals
  ```
  This runs for every cell in a multi-level region and makes two heap allocations each time.
- **Fix:**
  ```julia
  for r in allregs
      dz = elevations[r] - zval
      dz > 0 && (resvec[r] += dz)
  end
  ```
  If `allregs` is ordered top-down, `break` at the first `dz <= 0`, as `_compute_trap_footprints` does.
- **Also clean up:**
  - `res = Vector{Any}`, `Threads.@spawn` with `num_threads = 1` (`:52-88`). AGENTS.md says no threading without discussion, so drop this scaffolding.
  - `return Vector{Int64}()` vs Float64 elsewhere (`:37`) is type-unstable.
- **Note:** volumes are not scaled by cell area anywhere here. Check this is intended (cell-unit volumes).
- **Fusion opportunity:** `trapvolumes` and `_compute_trap_footprints` each do a full pass with the same per-cell logic and could be one pass.

### S6. `spillfield` temporaries and the unused slope grid [R]

- **Where:**
  - `spillfield.jl:79-80` allocates `slope::Matrix{Float64}`, which `spillanalysis.jl:100` discards.
  - `_find_downslopes` / `_spillfield!` (`:244-256,421-431`) allocate two `_diffgrid` results and three `_compare_slopes` pairs per call, plus logical-index copies. That is about 60–90 B/cell at peak.
  - With `clip_mask`, `:84-86` makes another full grid copy.
- **Fix:**
  - Use one fused loop over cells that checks the 8 neighbours and writes `dir::Int8` directly.
  - Add a `return_slope=false` keyword, defaulting to the current behaviour for API compatibility.
  - Set clipped cells to `Inf` once in `gridcpy`.

### S7. Spill-point hierarchy (`sshierarchy.jl`) [V]

- **Pairwise boundary merge:** `_identify_new_spillpoint!` (`:191-241`) runs a pairwise `setdiff` over the subtraps. That is O(k²·|B|), with a temporary reversed-tuple vector and Set per pair; k is large on flat or quantised DEMs.
  - **Fix:** mark the member leaf regions in a reusable `BitVector`. A boundary pair `(a,b)` is internal iff `mark[regions[b]]`. Filter and take the argmin z in one O(|B|) pass.
- **Boundaries never freed:** each higher-level trap pushes its whole outer boundary (`:223`), and the boundaries of merged subtraps are never released. Memory is O(perimeter × depth).
  - **Fix:** `empty!` the subtrap boundaries after the merge, if they are not needed later. Check `spillpoints` / `sshierarchy` callers first.
- **Stale comment** at `:201` ("by not copying…"): `setdiff` copies.
- **Cycle finding:** `Graphs.simplecycles` (`:62`) runs Johnson's algorithm on an out-degree ≤ 1 graph. A three-colour walk over a `next` array finds all cycles in O(R).
- **Climbing the hierarchy:** `_identify_highest_level_trap_of` (`:170`) climbs the hierarchy per call. Keep a `top_of` array, updated on merge.

### S8. Smaller static items [R unless marked]

- **Non-const global [V]:** `_SPILLFIELD_INACTIVE_FLAG = -10` at `spillregions.jl:8` is non-`const`, and it is compared per cell. Make it `const`.
- **Per-cell allocations:**
  - `dir in [-3, -1]` (`spillregions.jl:656`), and `tnodes = [...]` / `M[tn] ∉ [-1,-2]` (`:734-736`). Use scalar comparisons.
  - `trapvolumes.jl:69` (see S5).
- **Repeated `maximum(regions)`:** at `spillpoints.jl:123,129`, `sshierarchy.jl:49`, `trapvolumes.jl:33` and `spillanalysis.jl:131`. Pass the count down.
- **Double lookup:** `spillpoints.jl:146` checks `pos in keys(cut_edges) && ... in cut_edges[pos]`. Use one `get`.
- **`watercourses` temporaries (`watercourses.jl:50-112`):** 5–6 full-grid Float64 temporaries (`precipitation - infiltration`, `max.(-runoff,0)`, `infiltration - remaining`, and `fill`s for scalars). It also runs `topological_sort_by_dfs` over N vertices on every call. Fuse the loops and cache the order (see F5).
- **`flow_path_from`:** `x ∈ full_traps` on a Vector (`watercourses.jl:188`), and `setdiff(paths[end], footprints[t])` builds a Set of a whole footprint per step.

---

## 2. Dynamic sequencing: `fill_sequence/` (constant-rate path)

### F1. Full grid snapshots per weather period [V]

- **Where:** `fill_sequence.jl:102-103`.
  ```julia
  push!(seq, SpillEvent(cur_time, copy(cur_amounts), copy(filled_traps),
        copy(rateinfo.trap_inflow), copy(we.rain_rate), copy(rateinfo.runoff)))
  ```
- **Cost:** O(W·N). For example, 1440 one-minute rain steps on 10⁷ cells is about 115 GB for runoff alone, and double that with gridded rain. Events within a period are incremental and fine.
- **Invariant (AGENTS.md):** "The first event in each weather period is always a full snapshot." Changing this is a documented-invariant change. It affects `_reconstruct` (`spillevent.jl:188`), `trap_states_at_timepoints` / `network_states_at_timepoints` (`utils.jl`), and any user code reading `seq[i].runoff`.
- **Options:**
  1. Keep the full snapshot for the trap vectors (cheap, O(T)), but store runoff as a diff against the previous period's end state. `_reconstruct` already replays diffs; it would walk back to the *first* full runoff instead of the period start.
  2. Keep a full runoff grid only every k periods.
  3. Don't copy `we.rain_rate`; reference the `WeatherEvent`'s array, which the caller already owns.
- **Related [R]:** `getrunoffupdates` (`rateinfo.jl:215`) emits every cell touched since the savepoint, even when retract-then-add leaves it unchanged. Filter out entries equal to their stored value.
- **Related [R]:** `SpillEvent` fields are `Union`s (`spillevent.jl:78-92`), so every access branches at runtime. A flag plus concrete fields would avoid that. This is low priority.

### F2. Footprint rescans on every flow update [V]

- **Where:** `_track_flow!` (`flow.jl:137-204`), line 183, calls `_update_Smin_Smax!(rateinfo, tstruct, tstruct.supertraps_of[reg])` (`flow.jl:289`). That calls `_ponding_infiltration` (`flow.jl:266`) per supertrap, which loops over the full footprint.
- **Cost:** O(Σ F over the ancestor chain) ≈ O(N·depth) per call. `_update_flow!` calls it twice per graph update (once to retract, once to add), and a single event can trigger several.
- **Fix: incremental Smax.**
  - `Smax[t] = Σ_{c ∈ fp(t), z_c < sp_t} max(-runoff_c, 0)`.
  - When `setrunoff!` changes cell c from a to b, compute `δ = max(-b,0) - max(-a,0)`.
  - For each `t ∈ supertraps_of[regions[c]]` with `z_c < sp_t`, apply `Smax[t] += δ`. That is O(depth) per changed cell.
  - Then recompute `Smin` for the affected traps (sum over children).
  - Implement this in `_track_flow!`, which already knows `prev_value` and the new value, or in a variant of `setrunoff!` that has `tstruct` available.
  - The savepoint mechanism (`setsmax!`) must still record the original values.
- **Accuracy:** floating-point accumulation differs from a fresh sum, so checksum tests may drift in the last bits. Options:
  - periodic full recompute (e.g. per weather period, where `compute_flow` already does it);
  - compensated (Kahan) summation;
  - accept the change with a justification in the PR (AGENTS.md: don't change expected values without one).
- Run `Sequencing` and `Trapping structure`.

### F3. `fill_trap_until`: wasted allocations, and an unnecessary ODE [V]

- **Where:** `fill_sequence.jl:729-807`, with `_setup_dvdt` at `:820-842`.
- **Wasted allocations:**
  - `trap_bottom = tstruct.topography[footprint]` (`:733`) and `max.(trap_bottom, ...)` (`:742`) are computed *before* the `Smax == Smin` analytic early return (`:747`) and are unused on that path.
  - The function is called per inflow-updated trap per event, per branch-and-bound candidate, and per trap at period end.
  - **Quick fix:** move these lines below the early return.
- **ODE path costs:**
  - `_dvdt` allocates `trap_bottom .<= z` (a BitVector of length F) and `fprint_infil[ixs]` on every RHS call.
  - `condition` allocates `deriv=[0.0]` and calls `dvdt` again (another O(F)) on every callback evaluation.
  - `v2z` is either a closure or an interpolation object, so calls through it dispatch dynamically.
  - `solve` saves every step.
  - `condition_reached[1] ∈ [1,2]` allocates.
  - `fprint_infil` is built per call, using Dict `get`s when `use_saved`.
- **The ODE is unnecessary:**
  - `dV/dt = inflow − Σ_{c: bottom_c ≤ z(V)} infil_c` is piecewise constant in V, with breakpoints at the sorted bottom levels. Those are exactly the knots of `z_vol_tables[trap]`.
  - Precompute, per call or cached, the infiltration sorted by bottom level as a prefix sum `S_k`.
  - Then walk the segments from the current V: in segment k the rate is `r_k = inflow − S_k`, and the time to cross it is `ΔV_k / r_k`.
  - Stop at full, at empty, or where `r_k` changes sign. That last case is the "stagnation" condition, at a breakpoint.
  - This is exact (the current solve uses `abstol=1e-6, reltol=1e-4`), allocation-free apart from the prefix sum, and avoids `ODEProblem`, callback and integrator setup entirely.
- **Constraint:** AGENTS.md says `fill_trap_until` integrates with DifferentialEquations and that the three callback conditions (full, empty, sign change) must not be removed. The analytic version must reproduce all three outcomes and needs the maintainer's agreement to update AGENTS.md.
- **Expect** small numeric differences, since the result becomes exact instead of tolerance-limited. Justify any test updates.
- **Shared building block:** the same sorted prefix sum also serves D3 (network `wetted_infiltration`). Write one helper used by both.

### F4. Candidate scan over all traps every event [V]

- **Where:** `_identify_next_status_change!` (`fill_sequence.jl:511-590`).
  - `:534` loops `for ix in 1:num_traps` every event, O(E·T) in total.
  - `:557-560` do `findall(candidate_mintimes .> earliest)`, `sort!` and two `deleteat!` calls per branch-and-bound pass, O(C²) worst case.
- **Fix:** keep a min-heap keyed on `changetimeest[ix].min`, for example `DataStructures.MutableBinaryMinHeap` with one handle per trap. `DataStructures` is a new dependency, so check `Project.toml`.
  - Update the heap wherever `changetimeest` changes:
    - `_update_changetime_estimates!`;
    - pinning an exact value (`:549`);
    - fired traps set to Inf (`:572`);
    - children recomputed (`:583`);
    - network overrides (`_apply_network_changetimeest!`).
  - Pop while `top.min ≤ best_exact`, checking the `is_subsumed` and `all_subtraps_filled` gates at pop time. Reinsert candidates that survive.
- **Medium effort, moderate risk.** The tie handling (`==`) must be preserved exactly.

### F5. Per-weather-period rebuild in `compute_flow` [V/R]

- **Where:** `flow.jl:22-67` → `_compute_initial_rateinfo` → `watercourses`.
  - `topological_sort_by_dfs(tstruct.flowgraph)` (`watercourses.jl:61`) runs on every period: O(N) work plus allocations, although the graph never changes.
  - `_update_Smin_Smax!(…, 1:T)` costs O(Σ F).
- **Fix:** compute the topological order once in `fill_sequence` and pass it down. Per memory: no new `TrapStructure` fields, so add an optional argument to `watercourses`.
- **Oddity [V]:** `_compute_initial_rateinfo` calls `watercourses(tstruct, [false], …)`, a length-1 `full_traps`. It works only because `any([false])` is false. Pass `falses(numtraps(tstruct))` or add a dedicated entry point.

### F6. `trap_states_at_timepoints` rebuilds state per timepoint [R]

- **Where:** `utils.jl:139-169`. For each timepoint it calls `filled_at`, `runoff_at`, `inflow_at` and `amount_at`. Each walks back to the last full snapshot, copies it (O(N) for runoff), and replays all diffs since. Then each transitory trap is re-projected with `fill_trap_until`, which is O(F) because of F3.
- **Cost:** O(P·(N + events since boundary)).
- **Fix:** keep running state and apply diffs forward as the timepoints advance, as `network_states_at_timepoints` (`utils.jl:253-308`) already does. Also type `result = []` (`:139`).

### F7. Per-event small overheads [R]

- `setsavepoint!` → `_clear_saved_dicts!` (`rateinfo.jl:169-174`) allocates 4 new Dicts per event. Use `empty!`; this is safe because the update getters return fresh vectors.
- `getinflowupdates(rateinfo)` is rebuilt 3–4 times per event (`fill_sequence.jl:181`, `:215`, `:320`, `:426`). Compute it once per event and pass it along.
- `copy(net_covered_set)` (`fill_sequence.jl:168`) runs every event, even with no networks. Guard it with `isempty(driver.contexts)`.
- `filter(u -> u.index ∉ net_covered_set, …)` (`:162`) allocates per event.
- `_compute_exact_fill` (`:609`) uses the allocating `!all(filled_traps[subtrapsof(…)])`. Use `all_subtraps_filled` (`:477`).
- `@assert(all([ca.time == cur_time for ca ∈ cur_amounts]))` (`:82`) builds an O(T) array per period. Use a generator.
- `flow.jl:296`: `sum(getsmax(rateinfo, subtrapsof(…)))` indexes with a vector, so it allocates.
- `_is_sink` (`flow.jl:114-119`) builds `CartesianIndices` and does a linear `in` over the `sinks` Vector per call. The `sinks !== nothing` check is always true given the field type. Precompute a `BitVector` / `Set{Int}` of linear indices in `fill_sequence`.
- `SpillGraph.edges::Dict{Int,Int}` (`spillgraph.jl:17`): a `Vector{Int}` with 0 meaning "no edge" would make `_propagate_amount!`'s lookups plain array loads.
- `spillgraph.jl:95`: `full_traps[children_of(x)]` allocates per trap in `compute_complete_spillgraph`.
- `:180`: `all([…])` allocates.
- `:219`: `_siblings_of` returns `[]` (`Vector{Any}`).

---

## 3. Dynamic networks: `dynamics/`

### D1. Full dense ODE trajectory kept [V], the top quick win

*Fixed 2026-11-13 (`85ff9db`), in `fill_trap_until` too.*

- **Where:** `solveDynNetwork!`, `networksolver.jl:1672-1675`.
  ```julia
  sol = solve(ODEProblem(dynNetworkRateFunction!, V0, (0.0, tmax_ode), p),
              AutoTsit5(Rodas5P(autodiff = AutoFiniteDiff()));
              callback = callback, abstol = abstol, reltol = reltol)
  ```
- **Problem:** the defaults `save_everystep=true, dense=true` store `u` plus all stage vectors for every step. The comment at `:1666-1670` notes solves of about 10⁶ steps. Only `sol.u[end]` and `sol.t[end]` are read.
- **Fix:** add `save_everystep=false, save_start=false, dense=false`. Callbacks still interpolate within the current step. Alternatively use `init` + `solve!` and read `integrator.u` / `integrator.t`.
- Also check any other `solve(` calls in `src/` for the same issue.

### D2. Whole-driver recommit and rebuild on every touched event [V]

- **Re-checked 2026-11-12.** An earlier per-context gating (`76f9181` on `strip-zeno-machinery`, merged) predates the gate Phase C driver. The current code gates only **all-or-nothing**:
  - `_touch_affected_networks!` (`fill_sequence.jl:213-222`) returns early only if *no* context is touched.
  - If any one context is touched, `_touch_networks_driver!` (`network_driver.jl:283-323`) commits **every** context, one solve each (`:296-298`).
  - It then calls `_rebuild_contexts!` (`:127-139`), which discards all contexts and builds and re-predicts **every** component (`_make_context` + `_predict_network!`), another solve each.
  - The finding below therefore stands for runs with more than one network component. With a single network component there is no waste.
- **Where:**
  - `network_driver.jl:296-298` commits every context, one ODE solve each.
  - `:317` → `_rebuild_contexts!` (`:127-139`) rebuilds and re-predicts every context, another solve each.
  - The trigger in `fill_sequence.jl:217-221` is "any context touched".
  - `step_network_driver!` (`:213-243`) duplicates this logic.
- **Fix:**
  - Pass the indices of the touched components through.
  - Contexts that are neither touched nor fused keep `(state, last_solve_time, next_event)`, which stay valid because each solve starts from the context's own `last_solve_time`.
  - Only fused, regrown or detached components get new contexts.
  - Components whose `inflow_regions` overlap the changed inflows must also count as touched.
- **Related costs:**
  - `_build_rate_params` is called inside every `solveDynNetwork!` (`networksolver.jl:1617`). It rebuilds `_build_trap_geometry` (three footprint-sized arrays and an interpolation object per trap, `:72-108`), `_build_nbs_plan` (a `trap_of_cell` Dict over all footprint cells, `:857-860`), the order, the event templates and the culvert plan. Predict and commit build identical sets.
    - **Fix:** build `p` once in `_make_context` and store it on the context. Only `extern_inflow` changes per touch.
  - `Matrix{Float64}(precipitation)` (`network_context.jl:247-248`) copies the grid per context per rebuild, even though `driver.precip` is already a `Matrix{Float64}`. Store the reference.
  - `_reconcile_spillgraph!` (`network_context.jl:298-299`) rebuilds `compute_complete_spillgraph` from scratch, O(T). `_net_covered_set` / `_net_trap_set` walk all descendants of every node per event. `Set(findall(filled_traps))` (`network_driver.jl:244,318`) is O(T).

### D3. Network ODE right-hand side [R]

- **Untyped field:** `TrapGeometry.v2z` (`networksolver.jl:43`) holds either an interpolation object or a closure, so every `water_level` call (`:127`) dispatches dynamically. Store the knot vectors and interpolate with `searchsortedlast`.
- **`wetted_infiltration`** (`:147-156`) is O(F) per trap per RHS call, probably the largest CPU term. Use a sorted bottom array with a prefix sum of infiltration and look up `cumsum[searchsortedlast(zsorted, z)]`. Share the helper with F3.
- **`_attenuate_range`** (`:797-803`) walks every path cell per call. Each cell's map `d ↦ max(V+d,0) − max(V,0)` has the form `max(d + a, b)`, and that form is closed under composition: `max(max(d+a1,b1)+a2, b2) = max(d+(a1+a2), max(b1+a2,b2))`.
  - Precompute one `(a,b)` pair per segment between path events in `_path_event_templates`. Routing becomes O(events + nodes).
  - This must stay exact (mass conservation); the composition is exact apart from rounding.
- **About 8 allocations per RHS call:**
  - `:1123` `spilling = Bool[...]`
  - `:1126` `trap_level`
  - `:958-959,978` `NBSRouting` (three vectors)
  - `:545` `Float64.(external_inflow)`
  - `:547-549` `path_flow`, `trib_output`, `culvert_actual`

  Preallocate a scratch struct in `p` and `fill!` it. If the buffers are typed by `eltype(V)` (or use `PreallocationTools.DiffCache`), the RHS becomes Dual-safe, which removes the stated reason for `AutoFiniteDiff`.
- **Solver choice:** `Rodas5P(autodiff=AutoFiniteDiff())` with no `jac_prototype` means a dense Jacobian: n+1 RHS calls each, O(n²) memory, O(n³) LU. The true Jacobian is very sparse: the step functions have zero derivative almost everywhere, leaving only culvert heads and NBS layer ↔ carrier couplings. Options:
  - plain `Tsit5` when the net has no NBS (the stiffness comes from NBS overflow, per the comment);
  - a sparse `jac_prototype` built from the routing graph;
  - `linsolve = KrylovJL_GMRES()`.
- **Callbacks:**
  - `VectorContinuousCallback` (`:1438`) uses the default `interp_points=10`, and its condition runs a full `_routed_inflow` each time.
  - `:1418` routes even when there are no `:unspill` conditions.
  - The NBS steady-state callback (`:1375`) re-evaluates the RHS on each step; `get_du!` gives the same value for free.
  - Fixes: route only when `:unspill` conditions exist, and reduce `interp_points`.
- **Minor:**
  - `@assert` in the RHS (`:1164`, `:539-540`).
  - Symbol-keyed event kinds branch inside the condition (`:1420-1426`).
  - `DynNetworkContext.next_event::NamedTuple` is an abstract type (`network_context.jl:42`).
  - `NBSLayer.EVCoeff` is untyped (`nbs_elements.jl:34`).

### D4. Tracer, indexing and setup [V/R]

- **`_map_inverse_flow`** [V] (`build_network.jl:85-95`) builds a `Dict{Int,Vector{Int}}` over the whole grid. It is called unconditionally from `setup_network` through `_compute_nbs_inflow_outflow_cells!` (`:46`, `:109`), once per weather period and once per `network_watercourses` timepoint.
  - **Fix:** return early when `isempty(nbs)`. Otherwise find inflow cells by checking the 8 neighbours of each footprint cell, which is O(8·|footprint|).
- **Linear membership tests:**
  - `build_network.jl:161`: `filter(c -> !(c in tstruct.footprints[target_trap]), path)` is O(|path|·F). Strip only the trailing suffix inside the footprint, or use a mask.
  - `t in full_traps`, where `full_traps = findall(filled_traps)` is a Vector (`build_network.jl:146,156,238,245`; `network_updating.jl:328`). Pass the Bool vector instead.
  - `findfirst(dtrap -> dtrap.trap_ix == trap_ix, network.traps)` (`build_network.jl:199`) is O(nt²) during growth.
  - `findfirst(cells .== CI[isect_cell])` (`:225`) allocates a BitVector.
  - Culvert inlet/outlet `in footprint` (`:304-305`) scans `CI[footprint]`, a fresh copy.
  - `ds in footprint` (`:402`) on a Vector is O(F²).
  - `u.index ∈ ctx.global_ix` (`fill_sequence.jl:218`) is a linear scan.
- **`apply_fill!` indexing:** it calls `_index_components` (`network_updating.jl:349`, `:291-294`), which indexes every footprint cell of every component on each fill, although only path and NBS cells are needed. `grow_spill!` (`:244-247`) rebuilds the path-cell map on every fill. Make both persistent on the driver.
- **`network_watercourses`:** `:73-74` fills full-grid arrays even for scalar inputs, `:84` copies the runoff grid, and `:91` hits the `_map_inverse_flow` cost above.

---

## 4. Correctness issues

Fix these independently of the performance work. Each one needs a regression test that fails before the fix.

### Dynamic path

- **C1 [V]: wrong time origin in the submerged-child changetime.** *Fixed 2026-11-13, with C2: the estimate is timed from the parent's checkpoint; an already-empty parent exposes the trap at once.*
  - `_compute_changetime_estimate`, `fill_sequence.jl:383-391`, computes:
    ```julia
    parent_volume = cur_amounts[parent].amount
    ... ChangeTimeEstimate(false, cur_time + min_time, cur_time + max_time)
    ```
    `parent_volume` was recorded at `cur_amounts[parent].time`, which can be earlier than `cur_time`. The unfilled branch at `:399-405` correctly uses `starttime = cur_amounts[trap].time`. The children re-estimate at `:584` passes `earliest_changetime` and hits the same issue.
  - **Effect:** the lower bound can come out too late. Branch-and-bound then drops the candidate (`:557`), and the event is missed or misordered.
  - **Fix:** use `cur_amounts[parent].time` as the origin, or project the parent's volume to `cur_time` first. Compare with the earlier stale-time fix (memory `staletime-bug-fix.md`), which addressed the same pattern elsewhere.
- **C2 [R]: rate-only changes are not followed up.** *Fixed 2026-11-13. `b2bc8d7`: estimate refresh and amount checkpoints cover every trap whose inflow, Smax or Smin changed (`_rate_updated_traps`). Follow-up commit: a refreshed trap that is not full also refreshes its full subtraps, whose exposure time is derived from it. Tests in `test/fill_sequence_rates_test.jl`; the invariant "a trap holding water of its own has all its subtraps full" is now asserted per event (`_assert_subtraps_full_under_water`).*
  - `_track_flow!` updates Smax/Smin for all `supertraps_of[reg]`, but sets inflow only if the flow reaches a trap bottom (`flow.jl:188-190`).
  - When infiltration absorbs the flow en route, the trap's rates change but `_update_changetime_estimates!` (keyed on inflow updates only, `:426`) doesn't refresh its estimate.
  - Cached exact estimates (`min==max`, `:638`) are returned unchanged.
  - No amount checkpoint is taken, so a later `fill_trap_until` projects the whole interval since the last checkpoint at the new rate.
  - Submerged children of a parent whose rates changed are not refreshed either.
  - **Fix:** track Smin/Smax keys in `RateInfo`'s stored dicts (they already are: `stored_Smin_values` / `stored_Smax_values`). Treat those keys like inflow keys for estimate refresh and amount checkpointing, and refresh the children of affected parents.
- **C3 [V]: slow network flow silently dropped (mass loss).** *Fixed 2026-11-13: commits and timepoint samples run to `tmax` (`stop_at_steady = false`); traps settle only at a zero rate or a sign change, NBS layers below a relative rate tolerance and not by a sign change. Remaining: (a) in a prediction, a trap in a culvert/NBS network whose rate crosses zero smoothly (a real turn, not a stall) is still taken as settled, so a later event of that network can go unpredicted until the next touch; (b) new finding: the solver's fill-event time is off by ~4e-4 to 1e-3 relative on a single constant-rate trap (mini.txt trap 57), more than its tolerances suggest — not investigated.*
  - `solveDynNetwork!` returns early with `time=Inf`, without advancing `state`, when every `|du0| ≤ abstol` (`networksolver.jl:1648-1649`). The steady-state callbacks (`:1339-1341`, `:1378-1379`) can also stop a *commit* solve before `tmax`.
  - `_commit_network!` (`network_context.jl:146-153`) still sets `last_solve_time = T_commit`, so sub-`abstol` net flow over the window is lost.
  - `abstol` is a volume tolerance reused as a rate threshold, so the units don't match.
  - The "sign flipped" settled test can accept large `|rate|`, e.g. an NBS lower layer that drains and then refills.
  - **Fix:** in commit mode (finite `tmax`), disable the early return and `cb_ss`, or advance analytically with `state .+= du * remaining`. Use a dedicated, scale-aware rate tolerance.
- **C4 [R]: commit result ignored.** `_commit_network!` returns the solver result, but no caller checks it (`network_driver.jl:214,297`; `network_context.jl:393`). A `:fill` firing before `T_commit` leaves state clamped at an earlier time but labelled `T_commit`. Handle the event, or assert that none occurred.
- **C5 [P]: re-rooting drops co-located tributaries (mass loss).**
  - `_promote_tributary!` (`network_updating.jl:206`):
    ```julia
    mg = vcat(Q.merges, [(m, reb(p)) for (m, p) in P.merges if p != j])
    ```
    This drops *every* tributary merging at junction `j`, not just the promoted `q`. Two tributaries meeting at the same D8 cell is common.
  - **Effect:** the other tributary keeps `target_trap = 0`, so its flow is routed out of the domain.
  - **Fix:** filter on `m != q`, but first check that `Q.merges` doesn't already cover it. Write a test with two tributaries joining at one cell, then a re-root.
- **C6 [R]: culvert draw not capped by available water (mass creation).**
  - `_route_trap_node!` (`networksolver.jl:474-490`) subtracts `q = _culvert_flow(...)` from `trap_inflow[i]` with no cap.
  - The RHS floor (`:1179`) zeroes `dV` at V ≤ 0 while the outlet still receives `q`.
  - It is triggered for parent traps, where `zmin` is raised to the children's spill level (`:82`) but `inlet_invert` is raw terrain, giving a positive head at V=0. It is also triggered by overshoot to V<0.
  - **Fix:** mirror the path rule `min(capacity, available)` (AGENTS.md "Mass conservation"), and zero the draw at V ≤ 0.
- **C7 [R]: NBS leaks.**
  - `compute_outflow` uses `max(S − Smax, eps())^n` (`nbs_elements.jl:263`), so an empty layer emits `K·eps^n`. That is about 1.5e-8·K for n=0.5.
  - `:1194`: `max(O_0_total + nbs_draw, 0)` discards a negative remainder.
  - `_promote_tributary!` silently drops a culvert inlet in the cut stub (already flagged with `@@@` at `:196`).
  - `_update_pathmap!` (`build_network.jl:272-279`) splices the vector while iterating over it. Add a `break`.
- **C14 [V] — FIXED 2026-11-12 in `d0e296a`: a network commit read the runoff grid as updated for the event that triggered it (mass loss).**
  - **Found** while reproducing a SUrbArea webapp NBS anomaly, on the synthetic surface (`SUrbArea/debug/synthetic_dtm*`) with three green roofs (`synthetic_userdefined_2.geojson`, 120 mm/h for 1000 min).
  - **Symptoms:**
    - The global balance was 7.2% short (about 97 l/s) between 260 and 480 min.
    - Trap 3 (the downstream pond) was replayed as full from about 257 min but only marked filled at 483 min, so in the viewer its inflow vanished and its spill path appeared late.
    - Splitting the same rain into two weather periods gave the right answer.
  - **Cause:**
    - `DynNetworkContext.runoff` was a reference to the live `rateinfo.runoff`.
    - `_fill_sequence_for_weather_event!` runs `update_spillgraph!` / `_update_flow!` for an event *before* `_touch_affected_networks!` commits the interval leading up to it.
    - The commit therefore combined the cached pre-event `extern_inflow` with the post-event runoff. Here that was trap 4's new 11.09 m³/min spill crossing roof 1, whose downpipe drains elsewhere.
    - The roof's signed diff removed water that never arrived. Trap 3's rate fell to exactly 0 (5.76 with the correct grid), the solver stopped as "steady" (see C3), and ~1305 m³ already in trap 3 was lost.
  - **Fix:**
    - `DynNetworkContext` now holds `params::DynNetworkRateParams`, built once in `_make_context`. It freezes the inflow and everything read from the runoff and precipitation grids, and replaces the old fields `extern_inflow`, `runoff` and `precip`.
    - A new `solveDynNetwork!(state, tstruct, p; tmax, abstol, reltol)` method solves from prebuilt params; predict, commit, finalize and the replay sampling all use it.
    - The now-unused `infiltration` / `z_vol_tables` arguments were dropped from `_predict_network!`, `_commit_network!`, `_finalize_networks!`, `finalize_network_driver!`, `_carry_nbs!` and `_sample_network`.
    - This also delivers part of D2: the params are no longer rebuilt per solve, and the precipitation grid is no longer copied per context.
  - **Test:** "NBS: a network commit is not altered by the event that triggers it" (`test/nbs_routing_test.jl`).
    - The setup is a two-pit plane with a green roof on the upstream pit's spill path.
    - The downstream pit's fill time must be the same with and without a period split just before the upstream pit fills.
    - Before the fix: 6170.6 vs 6057.2 min. After: identical.
  - **Result on the synthetic case:** trap 3 fills at 256.6 min, and the balance closes. Repro scripts: `scratch/nbs_debug/` (`run2.jl` for the balance, `run9.jl` for the proof of the stale-runoff mismatch).
  - **Not a touch-gating problem.** The change comes from the same event, which does touch the network.

### Static path and API

- **C8 [V]: `lengths` cannot be used.** *Fixed 2026-11-13 (`085addb`).* `lengths::Union{Tuple{<:Real}, Nothing}` (`spillanalysis.jl:61`, `spillfield.jl:68,296,369`) is a *1-tuple* type, so `(Lx, Ly)` raises a MethodError. Change it to `NTuple{2,<:Real}`.
- **C9 [V]: integer DEMs crash.** *Fixed 2026-11-13 (`abdcda8`): grids are converted to float.* `_diffgrid` allocates `Array{eltype(grid)}` and then does `target[:] .= NaN` (`spillfield.jl:495-496`). That throws `InexactError` for an `Int` grid, yet `Matrix{<:Real}` is accepted. `TrapStructure{eltype(grid)}` would also fail on Float64 volumes. Convert the input to float in `spillanalysis`, or restrict the signature.
- **C10 [V]: `_setsamesize!` does nothing** *Fixed 2026-11-13 (`8b53b23`): removed.* (`spillfield.jl:560-564`); it rebinds a local variable. Find the callers and either return the array or delete the function.
- **C11 [R]: `domain=` and tiling are broken.** *2026-11-13: the options now raise "not fully implemented yet" (part of a suspended parallelisation effort, to be revisited).*
  - `dir` is sized to the domain but written with global indices (`spillfield.jl:79-80,255`).
  - `tiles()` is annotated `::Vector{Domain2D}` but returns a Matrix (`domain.jl:90-97`).
  - `spillpoints.jl:100` uses an undefined `return_boundaries`.
  - `update_spillregions!` starts with `@assert false` (`spillregions.jl:145`).
  - `_process_extended_domain!` calls a function with the wrong number of arguments (`:191`).
  - The tiling tests are commented out.
  - The docs already say domain is "not fully supported". Either fix it or make `domain`/`tiling` raise a clear error.
- **C12 [R]: traps with no outer boundary.** *Fixed 2026-11-13 (`d1c1ad4`): such a basin keeps spill elevation and volume `Inf` (it never fills, consistent with infinitely high clipped cells); its z/volume table, which ended in `NaN` and so ignored its spare infiltration, now extends one unit above its highest cell with the full footprint area.* `_identify_new_spillpoint!` returns `Spillpoint()` with elevation `typemax(T)` (`sshierarchy.jl:224-228`), and `trapvolumes` then computes `typemax − z`. Treat such traps as closed basins explicitly.
- **C13 [R]: out-of-domain spill targets dropped silently.** *Fixed 2026-11-13 (`bc00a90`).* `SimpleDiGraph(Edge.((i, downstream_region)))` with `downstream_region ≤ 0` (`sshierarchy.jl:57`) is skipped by Graphs without warning. It works only by accident; filter explicitly.

---

## 5. Cleanup and cruft

- **Unconditional prints [V]:** `flow.jl:31` (`println("compute initial rateinfo…")`, once per weather period) and `spillregions.jl:111,112,280,338`. Route them through `verbose`.
- **Dead code:**
  - `toptraps` [V] (`spillanalysis.jl:119-124`)
  - `regioncells` (`TrapStructure.jl:111`)
  - `_fix_boundary_seams!`, `_compute_region_correspondences!`, `_ensure_unique_region_numbering!`, `_update_correspondences!` (tiling only)
  - `compute_spillfield_graph`
  - `tull` (`utils.jl:1773-1774`)
  - commented blocks in `flow.jl:101-105,169-177,192-202`, and an unused `icount` with `@@@` at `flow.jl:41`
  - `num_traps`/`num_regions` unused at `flow.jl:27-28`
- **Unused parameter:** `time_slack` is accepted but not used (`fill_sequence.jl:29,129`), and its doc text "Support for this currently verified" is garbled. Implement it or remove it; removal is a keyword-interface change.
- **Unreachable branches:** the `Real` branch at `fill_sequence.jl:60`, and `sinks == nothing` at `spillanalysis.jl:79`.
- **`vcat()` on empty `nbs`:** `fill_sequence.jl:63` indexes with `vcat()` when `nbs` is empty, while `utils.jl:248` guards the same call. Make them consistent.
- **`_remove_outedge!`** (`spillgraph.jl:46-54`) has both an `@assert` and an `if` for the same condition. With asserts disabled it returns `nothing`, which is then pushed.
- **Exports and imports:**
  - `rateinfo.jl:3` exports a `copy` method that doesn't exist.
  - Unused imports: `SparseArrays`, `center` (`spillfield.jl:1-2`).
  - Commented-out Infiltrator imports (`utils.jl:3`, `spillregions.jl:3-4`).
- **Stale docs and comments:**
  - "Spillpoint.elevation is ::Real" (`fill_sequence.jl:768`, `flow.jl:267`); `Spillpoint{T}` is parametric now.
  - The `stored_inflow_values` docstring says "Smax" (`rateinfo.jl:38-41`).
  - `spillgraph.jl:14` writes `edges[i]=k` but means `j`.
  - The `TrapStructure` docstring omits `trap_bottoms` and `waterbodies`.
  - `networksolver.jl:43` "perf is a later concern" (that is D3).
  - `networksolver.jl:1070` "(tributaries only when there are no culverts)" is stale.
  - `build_network.jl:206` "@@ handle spill_path index later" is stale.
  - `sshierarchy.jl:201` "by not copying" is wrong.
- **Dead work:** `build_network.jl:203` computes `nbs_out` for traps and discards it.
- **Duplicate logic:** `step_network_driver!` is used only by `test/network_driver_test.jl` and duplicates `_touch_networks_driver!`. Retire it or make the test use the real path.
- **`reverse!.(supertraps_for)`** (`spillanalysis.jl:231`) allocates a throwaway result vector. Use `foreach(reverse!, …)`.

---

## 6. Recommended work plan

> **Priority decision (user, 2026-11-12):** get the **NBS code right first**, meaning correct and reasonably efficient. Only then turn to SWIM performance overall. Do the NBS track below before §6.1 onward.

### 6.0 NBS track (do first)

NBS state as of 2026-11-12:
- The last NBS commit is `e0157db` (2026-07-18); nothing NBS-related has been committed since.
- The current design is in `agent/NBS_ROUTING_NODE_PLAN.md`. The router-element model predates the 2026-07-18 fixes (area factor, exit weights, `runoff_orig`, replay carry, stiff solver).
- `agent/NBS_REMAINING.md`, `NBS_INTEGRATION_PLAN.md`, `NBS_OPTION1_OVERLAY_PLAN.md` and `NBS_ROUTING_REDESIGN.md` describe superseded designs. `NBS_REMAINING.md` also references a test file that no longer exists, `test/nbs_overlay_test.jl`.

Steps:
1. **Baseline.** Run the NBS-relevant suites in isolation (`nbs_routing_test.jl`, `dynamics_test.jl`, `network_driver_test.jl`, `network_states_test.jl`, `network_watercourses_test.jl`, `dynamic_membership_test.jl`) and confirm they are green.
2. **Correctness, test-first.**
   - C7: empty-layer `eps()^n` leak, and the discarded negative remainder at `networksolver.jl:1194`.
   - C3: the steady-state / early-return mass loss, including the NBS settled test that wrongly accepts a draining-then-refilling layer.
   - C4: unchecked commit result.
   - Also check C5 and C6, since NBS paths share the router.
   - Add a whole-network mass-balance test with NBS on a multi-cell footprint.
3. **Reasonable efficiency.**
   - D1: solver save flags.
   - D4: skip `_map_inverse_flow` when there is no NBS; otherwise do a local neighbour scan.
   - D2: cache `p`, including `_build_nbs_plan`, on the context.
   - D3: scratch buffers for `NBSRouting`, and a ForwardDiff-compatible RHS so the dense `AutoFiniteDiff` Jacobian can go. Use plain Tsit5 when a network has no NBS.
   - D2: per-component touch gating.
4. **Doc cleanup.** Update `NBS_ROUTING_NODE_PLAN.md` to the as-built state, and retire or delete the stale NBS documents.
5. **Then decide** whether the deferred NBS features (ET, real cell area, separate per-layer outlets, NBS→NBS in series) come before the general performance work.

### 6.1 General SWIM performance (after the NBS track)

Each step should be its own branch and PR. Run the relevant test sets after each (§7).

1. **Quick, local, behaviour-preserving wins.** Checksums should be unchanged.
   - D1 save flags.
   - F3 part 1: move the O(F) allocation below the analytic early return.
   - S1 quick fix: concrete `SimpleDiGraph{Int}` field types.
   - S5 `trapvolumes` scalar loop.
   - S4 two-pass waterbody flattening.
   - S8 `const` global and per-cell array literals.
   - D4: `_map_inverse_flow` early return when there is no NBS.
   - Remove the unconditional prints.
   - F7 items.
2. **Correctness fixes** C1, C3, C5, C6, then C2, C4, C7, each test-first. C8–C10 are small API fixes.
3. **Shared sorted-prefix-sum helper,** then:
   - F3 part 2: the analytic `fill_trap_until`. Needs maintainer approval to change the AGENTS.md ODE statement.
   - D3 `wetted_infiltration`.
4. **F2 incremental Smax/Smin,** with a strategy for floating-point drift.
5. **Network driver:** D2 (lazy untouched contexts, `p` cached on the context, precipitation by reference), then the rest of D3 (scratch buffers, Dual-safe RHS, solver and Jacobian choice, callback cost), then D4 tracer membership.
6. **F1 diff-based runoff and rain storage** in `SpillEvent`. Changes a documented invariant; discuss first.
7. **F4 heap,** F5 cached topological order, F6 running state.
8. **Large refactor, needs discussion:**
   - S2 union-find construction;
   - S1 `downstream::Vector{Int32}` replacing `flowgraph`;
   - Int32 `regions`;
   - S3 compact footprints and lazy z/volume tables;
   - S6 fused spillfield kernel;
   - S7 hierarchy improvements.

   This changes the `TrapStructure` layout and public fields: bump the minor version (AGENTS.md) and respect the "no new fields" rule.

---

## 7. Verification protocol

- **Baseline first:** before each step, record `@time` / `@allocated` and peak RSS (`Sys.maxrss()`) on the largest artifact grid, for `spillanalysis` and for `fill_sequence` with and without networks.
  - Scripts already exist under `agent/reports/perf_*.jl` (`perf_overhead_two_grids.jl`, `perf_full_coverage.jl`, `perf_localized_vs_spread.jl`, `perf_alloc_with_without.jl`).
  - Run with `exec taskset -c 0,2,4,6,8,10,12,14 julia --threads=8 .` (memory note).
- **Tests:** `Pkg.test()`. Tests are deterministic checksums.
  - Pure performance steps must leave them bit-identical.
  - Steps that change numerics (F2, F3 part 2) or fix bugs (C*) need a written justification for any expected-value change (AGENTS.md).
- **Minimum test sets:**
  - after `spillfield`/`spillregions` changes: `SpillfieldTests`, `SpillregionTests`;
  - after `fill_sequence/` changes: `Sequencing`, `Trapping structure`;
  - dynamics: run the dynamics tests in isolation.
- **Ignore as known failures:** `MaskingTests` (deprecated) and the pre-existing bay-grid drift / `raise_buildings!` failures (memory notes). They are not regressions.
- **Mass conservation:** for C3, C5, C6 and C7, and any D3 routing change, add a test that integrates total inflow − infiltration − outflow and checks the stored-volume change.
