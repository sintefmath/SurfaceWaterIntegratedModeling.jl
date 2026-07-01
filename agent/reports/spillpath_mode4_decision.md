# `_subnetwork` fix result + mode-4 `spill_path` decision

## 1. Fix applied & tested

The `_subnetwork` terminal-detection fix (handle the wrap-around / trap-terminated
case) is in `src/dynamics/elements.jl`, with two new helpers `_spills_out_of_domain`
and `_unfilled_parent_of`.

Verification:

- Both B1 probe cases now collapse to `[414]` (were `[9,18]` before):
  - (a) start inside child 9 → `traps=[414]`, no paths
  - (b) start upstream on slope → `traps=[414]`
- Full dynamics testset: **21 testsets, all pass, 0 failures** (B1 testset 161/161).

---

## 2. Mode-4 `spill_path = 0` tension — decision

### Key finding: the consuming code already handles `spill_path <= 0`

`_route_flow` (networksolver.jl):

```julia
if spilling[i]
    spill = max(trap_inflow[i] - footprint_infil[i], 0.0)
    sp = net.traps[i].spill_path
    sp > 0 && (path_flow[sp] += spill)            # sp == 0: spills out of domain
end
```

`_network_order` (networksolver.jl):

```julia
trap.spill_path > 0 && Graphs.add_edge!(g, np + t, trap.spill_path)
```

Both already treat `spill_path <= 0` as "terminal, drop the spill, it exits the
domain" — the routing comment literally says so, and a node with no out-edge is
fine for the topological sort. **The only place that rejects a FULL trap with
`spill_path == 0` is `_validate_network`.** So neither option is disruptive — it
is essentially a one-function decision.

### What you'd be giving up with bare `0`

`spill_path == 0` currently does double duty as an **integrity tripwire**: it is
how the validator catches "a trap filled but the network wasn't rebuilt to give
it a spill path." That bug and the legitimate domain-exit case become
indistinguishable if you blanket-allow `FULL + spill_path == 0`.

### Why `spill_path = -1` is the better option

| state | spill_path | meaning |
|---|---|---|
| TRANSITORY | `0` | not full yet, no spill path |
| FULL, interior | `> 0` | spills into that path |
| FULL, exits domain | `-1` | spills off the edge (drop it) — **explicit** |

With `-1` the validator becomes a clean 3-way check
(`FULL ⇒ spill_path ≠ 0`, `TRANSITORY ⇒ spill_path == 0`), so the
"forgot to rebuild" tripwire is **preserved** — a FULL trap with `0` is still an
error. And because every consumer already guards with `> 0`, `-1` flows through
routing and ordering with **no change** (no edge added, spill dropped) —
identical behavior to `0` there, just self-documenting.

### Cost of `-1` (small, localized)

1. `_validate_network`: relax FULL to allow `-1`, keep rejecting `0`.
2. `_build_network` / `_subnetwork`: emit `-1` (not `0`) for the mode-4 terminal
   trap — i.e. when the terminal full trap `_spills_out_of_domain`. The `link(...)`
   helper currently yields `0`; special-case the domain-exit terminal.
3. Update the `valid_network` test helper and the two `_build_network` tests
   (dynamics_test.jl lines 54/72) that currently assert `spill_path == 0` for
   "trap exits domain".

(Third option — synthesize a degenerate exit `DynFlowPath` with `target_trap = 0`
so the `>0` invariant holds literally — works too, but a fake one-cell path is
less honest than an explicit sentinel.)

### Recommendation

Go with **`-1`**. It keeps the safety check, needs zero routing/ordering changes,
and makes the intent legible.

---

## Next step (awaiting your go-ahead)

Implement `-1`:
- validator 3-way check
- `-1` emission in `_build_network` / `_subnetwork` for the domain-exit terminal
- corresponding test updates (`valid_network`, the two `_build_network` asserts)
