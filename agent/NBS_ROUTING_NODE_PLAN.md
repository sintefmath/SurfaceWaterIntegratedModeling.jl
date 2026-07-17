# NBS as a signed-diff router element — as-built design

> **Status: IMPLEMENTED** on branch `nbs-dynamic-core` (commits `2ff9f6a` — add
> `DynFlowPath.nbs_inlets`; `fa21568` — the router rewrite; `c4bb889` — this doc).
> This file now documents what was *built*, not a forward plan. Where the built code
> diverged from the original plan, the divergence is called out inline (**BUILT:**).

## Core idea

`watercourses` stays **NBS-oblivious**: NBS installations cannot own watersheds, they
intersect existing ones. So the oblivious overland flow still passes through the footprint
and is already routed downstream by the ordinary machinery. The NBS is a **router element**
layered on top of that baseline:

1. It captures its **live total input** `I_1` (oblivious background throughput + any network
   flow crossing the footprint), runs the layered-storage ODE (`NBSSystem`), and
2. emits only the **diff** `O_1 − O_0` of its live output over what the oblivious routing
   already delivered. The diff rides the ordinary flow router (extended to signed values),
   attenuated against the read-only oblivious `runoff` grid, and lands in the downstream trap
   inflow — where it cascades onward via the normal trap spill.

Net downstream `= O_0 (baseline) + diff = O_1`. `watercourses` untouched; no double-count.

## The oblivious runoff grid (`rateinfo.runoff`) is the reference

`compute_flow` runs on the `FULL ∧ ¬COVERED` spillgraph, so `runoff` is
**network/NBS-oblivious**: per cell, **net overland flow (positive)** or **remaining
infiltration capacity (negative)**, static within a solve (see `RateInfo` docstring). The
grid and the network router are therefore *complementary* — the grid excludes network traps,
the router carries only network flow — so quantities read from each can be summed without
double-counting. `runoff` is threaded into the solver and read **read-only** (never mutated,
no per-step copy): as `path_runoff` (per-path cell values, for diff attenuation) and to
precompute each NBS element's endpoint sums.

## Input `I_1` (drives the layer ODE)

```
I_1  =  O_0_total            # oblivious background throughput (static; = inflow + on-footprint rain)
      + Σ :nbsin draws       # network flow crossing the footprint (dynamic, from the router)
```

Layer-1 (top) inflow in the rate function is `el.O_0_total + nbs_draw[k]`.

**Why `O_0_total`, and not a separate `Σ runoff[inflow_cell] + precip` sum** — this is the
key subtlety. The footprint has **zero infiltration** (the contract), so by mass balance
every drop entering (boundary inflow + rain falling on the footprint) must leave at a
boundary exit or pond at an internal sink:

```
inflow + on-footprint rain  =  Σ(boundary-exit flow)  +  Σ(internal-sink ponding)  =  O_0_total
```

`_build_nbs_plan` sums `runoff` over **every** footprint drainage endpoint — boundary exits
(`ds[1] ∉ footprint`) *and* internal depressions (flowgraph sinks) — skipping interior cells
that drain to another footprint cell. So each drop is counted once, at its exit/sink, with
its full accumulated flow, and `O_0_total` equals the **total live background input**. It only
*looks* like an output name.

- **Precipitation is included** — folded into `O_0_total` via the `runoff` grid (`compute_flow`
  accumulates rain), *not* as a separate `+precip` term. A separate term would double-count
  the on-footprint rain that is already in the endpoint runoff.
- `:nbsin` draws are the *network* flow (upstream network spills, upstream outlet discharge)
  arriving on `DynFlowPath`s. The oblivious grid excludes these, so `O_0_total + draws` does
  not double-count.

**BUILT vs plan:** the plan proposed `I_bg = Σ runoff[footprint_inflow_cells]` (boundary
inflow only) plus a separate precip term. That omits on-footprint rain and broke pass-through
(see *History* below). Replaced by `O_0_total`, which is mass-consistent with the output
baseline.

### `:nbsin` capture

`DynFlowPath.nbs_inlets` marks the footprint-inflow boundary cells a path crosses. In
`_path_delivered!`, a `:nbsin` event **consumes** the passing router flow:
`nbs_draw[idx] += current; current = 0`. The path then continues through the footprint
carrying 0 (interior cells have zero infiltration, so they do nothing). `nbs_draw` comes back
from `_route_flow` and feeds the layer ODE.

## Output diff

`O_1 = O_terrain(V) + E(V)` from the live layer state (`O_terrain` = overflow of the top
`n_terrain` layers; `E` = piped-layer discharge). The NBS emits only the diff from the
oblivious baseline `O_0`, distributed across the same endpoints the oblivious flow used:

- **terrain endpoint** `e`: `diff_e = (O_terrain(V) − O_0_total) · ratio_e`,
  `ratio_e = O_0[e] / O_0_total` (guard `O_0_total ≈ 0` → even split over endpoints). New total
  at `e` `= O_0[e] + diff_e = O_terrain · ratio_e`.
  - *boundary exit*: head-injected on the carrier path departing from the exit cell.
  - *internal depression*: deposited straight into the accumulation trap (no path).
- **piped outlet** `l`: `+E_l(V)`, head-injected on the carrier path departing from the
  outlet cell.

The diffs total `O_1 − O_0_total`. The terrain split math equals the old model's `−X·Q_c`
(`Q_c = O_0[e]`); what changed is the dynamic, mass-consistent input `I_1` and that the diff
now rides the router.

**BUILT vs plan — head injection, not `:nbsout` events.** The plan emitted the output diff as
a mid-path `:nbsout` event (a path crossing the output cell, culvert-outlet style). But every
NBS output cell — terrain exit cells (`footprint_outflow_cells`) and piped outlet cells — is
**already a seed** in `build_network`, so a `DynFlowPath` departs *from* it (`departure_point ==
that cell`) and already runs to the target trap. The diff is therefore injected at the **head**
of that seed path (`NBSRouting.path_diff[carrier]`), which is equivalent and needs no `:nbsout`
event kind and no output-side tracking. `_build_nbs_plan` maps each endpoint to its carrier via
a `departure_point → path` dict. The `DynFlowPath.nbs_outlets` field is retained only for
component coupling in `_nbs_coupled_nodes`, not for routing.

## Signed-diff attenuation

**One rule for all router flow (commit `f22ce5e`).** Every flow the router carries — network
trap spills, tributary/culvert additions, and the NBS output diff — is a **signed** value
attenuated per cell against that cell's oblivious residual `runoff` value `V` by the single rule
(`_attenuate_range` / `_attenuate_diff`, read-only):

```
out = max(V + d, 0) − max(V, 0)
```

| `d`  | cell `V`        | effect                                                            |
|------|-----------------|-------------------------------------------------------------------|
| `+`  | `+` (flow)      | passes unchanged (background already saturated the cell — no re-charge) |
| `+`  | `−` (capacity)  | remaining capacity `\|V\|` re-fills first, absorbing it            |
| `−`  | `+` (flow)      | the flow absorbs it, capped at `−V` (can't remove more than is there) |
| `−`  | `−` (no flow)   | moot → 0                                                          |

This is the **same** rule `_track_flow!` uses to build the oblivious grid (`flow.jl`); the
router now mirrors it instead of the old infiltration-prefix (which charged network flow the
full per-cell infiltration even on cells the background had saturated). `path_runoff` is the
real grid where available, else `−infiltration` (cells at their capacity floor), which
reproduces the prefix behaviour exactly for hand-built nets with no background.

The NBS output diff is **head-injected** into the path's signed `current` (`head = path_flow[p]
+ path_diff[p]`) and rides this one accumulator to `target_trap`; the normal trap spill carries
it onward, closing the "doesn't cascade past the first trap" gap. Because it shares `current`, a
same-path downstream `:nbsin` intercepts the whole signed flow (`nbs_draw += current; current =
0`) — an upstream NBS's positive release *or* negative storing deficit both fold into the
downstream NBS's input `I_1`. (All four sign cases verified against
`agent/prompts/flowtracking_adjust.org`.)

## No double-count (option ii)

`external_inflow[trap]` keeps the oblivious baseline — `O_0`'s contribution as it reached the
trap. The router delivers only the **attenuated diff** (the change), so the trap ends at
`O_0 + diff = O_1`. Track the diff, not the total; `runoff` stays a read-only reference.

Because the oblivious total-in equals total-out-or-ponded, the number used as the ODE input
baseline and as the diff's output baseline is the **same** `O_0_total` — the mass-balance
identity is what makes a true pass-through NBS emit `diff ≈ 0`.

## Code map (as built)

`src/dynamics/networksolver.jl`:
- `NBSElement` — `system`, `A`, `state_base`, `n_terrain`, `O_0_total`, and the output routes:
  `terrain_paths::[(carrier_path, ratio_e)]`, `terrain_traps::[(acc_trap, ratio_e)]`,
  `piped_paths::[(carrier_path, layer)]`.
- `NBSPlan` — `elems`, `nlayer_total`.
- `NBSRouting` — per-step scratch: `path_diff`, `trap_extra` (internal-depression deposits),
  `nbs_draw` (router output = `:nbsin` captures).
- `_build_nbs_plan` — endpoint walk → `O_0`, `O_0_total`, `ratio_e`, carrier paths.
- `_nbs_routing` — per-step output diffs from the live layer state → a fresh `NBSRouting`.
- `_attenuate_range` (the one rule) / `_attenuate_diff` (whole-path), `_path_cell_values`,
  `DynNetworkRateParams.path_runoff` (always present — real grid or `−infiltration`).
- `_path_event_templates` emits `:nbsin`; `_path_delivered!` walks cells against `path_runoff`
  and handles the `:nbsin` draw; `_route_flow` seeds `trap_extra` and head-injects `path_diff`
  into the signed `current`; `_routed_inflow` builds the `NBSRouting` and returns `nbs_draw`;
  the rate function drives layer-1 with `O_0_total + nbs_draw`.
- **Deleted:** `NBSDelivery`, `_walk_to_trap`, `_apply_nbs_corrections!`,
  `_propagate_correction`, the output-side static `capture`, `NBSLayerParams`.

`src/dynamics/elements.jl` — `DynFlowPath.nbs_inlets::Vector{Tuple{Int,Int}}` (additive; no
`target_trap` change, no path termination). `build_network.jl` — populated via
`_intersecting_on_path` on `footprint_inflow_cells`; threaded through the membership/localize
reconstructions (`network_updating.jl`, `network_utils.jl`). `nbs_elements.jl` — `NBSLayer`
numeric fields typed `Float64`.

Untouched: `watercourses.jl`, `_network_order`, the membership/trace machinery.

## History — the spurious negative diff (fixed)

The first cut followed the plan literally: `I_bg = Σ runoff[footprint_inflow_cells]`, i.e. only
flow crossing the **upstream boundary**, omitting rain on the footprint. But the output
baseline `O_0_total` includes that rain, so `O_0_total > I_bg`. A fast pass-through NBS
(`O_1 = O_terrain ≈ I_1 ≈ I_bg`) then produced `diff ≈ I_bg − O_0_total = −(on-footprint rain)
< 0` — the NBS spuriously subtracted the footprint's own rain, filling the downstream pit ~20%
late. Root cause: input and output baseline were measured over *different water*. Fix: feed the
layer the full throughput `O_0_total`, restoring the `input = output` identity that pass-through
relies on.

## Deferred work

- **Submergence** — a containing trap flooding the footprint from below is **not** modeled; the
  footprint is treated like any terrain. To be revisited later (leaving it lie for now).
- **Internal-depression re-emit** — the output diff is distributed to internal sinks too
  (`terrain_traps`), so in pass-through an internal sink still ponds its oblivious share
  (matches the old model). Physically a storage NBS arguably captures internal ponding and
  re-emits it at the boundary; reconsider together with submergence.
- ~~NBS → NBS on the same path~~ **resolved** (router unification `f22ce5e` + signed `:nbsin`):
  the NBS output diff rides the single signed `current`, and a same-path downstream `:nbsin`
  intercepts the *whole* signed flow into its input — a positive release *and* a negative deficit
  (upstream NBS storing) both land in the downstream NBS's `I_1 = O_0_total + nbs_draw`, correct
  and non-negative (attenuation caps a deficit at the flow actually present, which is part of the
  downstream `O_0_total`). Unit-tested (`:nbsin intercepts the whole signed flow`).
- **Per-layer distinct outlets**, **evapotranspiration** (currently a `0.0` placeholder in the
  layer ODE), and **real cell area** (`A` = footprint cell count; `@@@ 1 m²/cell`) — the last is
  the biggest quantitative approximation (needs a resolution value threaded from the grid).
- **Output split by static oblivious ratio** (`ratio_e = O_0[e]/O_0_total`) rather than the live
  drainage split — an approximation while the distribution is stable.
- **Integration tests** — the terrain **cascade / network-inflow / upstream-outlet** tests are
  not yet written: hard to isolate on rain-flooded terrain (every pit fills from its own
  catchment before the cascade arrives), and the cascade path is generic router behavior
  already exercised by the driver suite.

## Hardening (done)

- **Silent-drop fallbacks → asserts.** A boundary-exit / piped-outlet carrier path and a
  trap-bottom internal-depression trap are all guaranteed by construction (output cells are
  seeds; accumulation regions map to exactly one net trap). `_build_nbs_plan` now asserts these
  instead of defaulting a missing lookup to `0` (which would silently leak mass). A domain-exit
  sink (`regions[f] <= 0`) is still a legitimate `0`-drop (its share leaves the domain).
- **`I_1` floor.** Layer-1 inflow is `max(O_0_total + nbs_draw, 0)` — physically `≥ 0` (a deficit
  is capped by attenuation at the flow present, part of the footprint's own `O_0_total`); the
  `max` only guards float noise.

## Verification

Green: `test/nbs_routing_test.jl` (16 — the `_attenuate_diff` V-rule unit test, signed `:nbsin`
interception (both signs), retention + pass-through integration, **layer-cascade mass
conservation** `Σ dV_layer = I_1 − O_surface − ground`), `dynamics_test.jl` (729),
`dynamic_membership_test.jl` (280), `network_driver_test.jl` (31), the `Sequencing` clean cases
(5). Package loads clean. Only the pre-existing bay-grid drift / `raise_buildings!` cases fail
(unrelated).

The mass-conservation test drives the real rate function post-solve (off the ODE hot path).
Still not written: cascade / network-inflow / upstream-outlet terrain integration tests, and a
routing-side ledger (delivered diff + en-route infiltration = emitted diff).
