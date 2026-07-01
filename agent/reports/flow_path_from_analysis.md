# `flow_path_from` changes & `_subnetwork` fix — analysis

## 1. Is the new `flow_path_from` sound?

Yes — both changes are individually correct, and each maps to a clean
termination mode.

### Change A — wrap-around pop (watercourses.jl lines 194-199)

When the trace re-enters a full trap already in `downstream_filled_traps`, it
pops the returning segment and breaks. This is exactly the "water pools among
the full children of a not-yet-full common parent" situation.

It's sound: every footprint cell of the unfilled parent `P` belongs to a *leaf
descendant's* region (your own B1 note — e.g. cell 2366 / region 18 of parent
414), and all those descendants are full, so tracing out of the last child
always lands back on an already-visited sibling. There is no terrain cell where
flow can "rest" inside `P`'s exclusive area, so wrap-around is the *only* way
this case can terminate. Popping the segment that loops back is right — it
carries no external flow.

### Change B — return the final out-of-domain segment (lines 213-221)

Setting `cur_node = downstream_region_cell` and continuing means a full trap
that spills onto terrain that *runs off the edge* now gets that final segment
traced and returned (terminating via the `cur_reg <= 0` break, mode 1). The
`current_region_cell == downstream_region_cell` guard correctly singles out the
*direct* boundary-spill case: `_add_outer_bounderies!` (spillpoints.jl:199) sets
`Spillpoint(0, cell, cell, …)` precisely there, so the equality test is a
faithful detector. Sound.

The `@assert abs(length(paths) - length(ftraps)) <= 1` and the strict
path/trap alternation still hold under both new modes (parity checked for
start-inside and start-outside in each).

### Caveat to flag (not introduced by these changes, but exposed by B)

The direct-out case (mode 4) leaves the terminal *full* trap with
`spill_path = 0`. `_build_network`'s own tests bless that ("trap exits domain"
→ `spill_path == 0`, dynamics_test.jl lines 54/72), but the solver's
`_validate_network` *requires* `spill_path > 0` for any FULL trap. So a
domain-exiting trap that is simultaneously full *and* a network node would fail
validation. Whether that combination can actually arise in the integration is a
design question — likely not for an interior network, but worth a conscious
decision rather than a validator error later.

---

## 2. What's broken in `_subnetwork`, and why

`_subsume_terminal_parent` was wired to fire only when `ends_with_path` is true
(chain ends on a path discharging into the unfilled parent). After Change A, the
terminal-unfilled-parent scenario **no longer ends on a path** — the wrap-around
pop makes it end on a *trap* (the last full sibling). For a trap-terminated
chain `ends_with_path` is `false`, so the subsume branch is skipped entirely.

Confirmed empirically on `mini.txt` (parent 414, children 9 & 18 full, parent
unfilled) — both B1 test cases regress:

| case | `setup_network` now returns | B1 test expects |
|------|------------------------------|-----------------|
| (a) start inside child 9 | `traps=[9,18]`, spill_paths `[1,0]` | `traps=[414]` |
| (b) start upstream on slope | `traps=[9,18]` | `traps=[414]` |

Not cosmetic: the old B1 bug is fully back, and trap 18 comes out FULL with
`spill_path=0` (solver-invalid).

### Termination-mode map (reference)

| Mode | flow_path_from break | chain ends on | `ends_with_path` | needed action |
|------|----------------------|---------------|------------------|---------------|
| 1 | `cur_reg <= 0` (line 178) | path (exits domain) | true | none (`_unfilled_trap_at`→0) |
| 2 | `isempty(full_supertraps)` (187) | path (into unfilled trap) | true | subsume if parent |
| 3 | wrap-around (194) **NEW** | trap (full sibling) | **false** | **subsume into unfilled parent (currently missing)** |
| 4 | direct-out (217) **NEW** | trap (full, exits domain) | false | none — but `spill_path=0` caveat |

---

## 3. Proposed fix to `_subnetwork`

Extend the terminal detection to also handle the trap-terminated (wrap-around)
case. The discriminator is the same one `flow_path_from` itself uses — did the
terminal trap spill *directly out of domain*? If not, it pooled into its
unfilled parent, which is what we subsume into.

```julia
traps = collect(ftraps)
ends_with_path = length(paths) > length(ftraps) ||
                 (starts_with_trap && length(paths) == length(ftraps))

# The unfilled trap the chain ultimately pools into (0 = exits the domain).
tix = if ends_with_path
    # chain ends on a path discharging into an unfilled trap (or out of domain)
    _unfilled_trap_at(tstruct, paths[end][end], full_traps)
elseif !isempty(ftraps) && !_spills_out_of_domain(tstruct, ftraps[end])
    # chain ends on a full trap that wrapped back among full siblings of an
    # unfilled parent (flow_path_from's wrap-around termination, Change A):
    # subsume that suffix into the parent.
    _unfilled_parent_of(tstruct, ftraps[end], full_traps)
else
    0
end

if tix > 0
    paths, traps, starts_with_trap =
        _subsume_terminal_parent(tstruct, paths, ftraps, tix, starts_with_trap)
end
```

with two small helpers:

```julia
# A trap whose spillpoint is the domain boundary cell (see _add_outer_bounderies!).
_spills_out_of_domain(tstruct, t) =
    (sp = tstruct.spillpoints[t]; sp.current_region_cell == sp.downstream_region_cell)

# Lowest ancestor of `t` that is not full (0 if none) — the basin water pools in.
function _unfilled_parent_of(tstruct, t, full_traps)
    p = parentof(tstruct, t)
    while p !== nothing && p ∈ full_traps
        p = parentof(tstruct, p)
    end
    return p === nothing ? 0 : p
end
```

`_subsume_terminal_parent` itself needs no change — for case (a) it finds
`d=1`, collapses to `[414]` with no paths; for case (b) it keeps the external
approach path and trims the footprint cells. Both traced by hand and reproduce
the expected `[414]`.

### Design notes

- Keep the fix in `_subnetwork` rather than have `flow_path_from` return its
  termination reason: `flow_path_from` is exported and has another caller
  (`utils.jl:729`); re-deriving the discriminator locally is cheap and avoids
  touching the public signature.
- The mode-4 direct-out / `spill_path=0` tension (§1 caveat) is deliberately
  left untouched here — separate decision.

---

## Open items for you

1. Apply this fix and re-run the B1 testset to confirm (a)/(b) go green?
2. How to handle the mode-4 direct-out / `spill_path=0` tension (full trap that
   spills straight out of the domain vs. the solver validator requiring
   `spill_path > 0`)?
