# Resume note — culvert drought fix (done) + runtime work (next)

_Working doc / note-to-self. Branch `SUrbArea`. Last updated 2026-07-02._

## UPDATE 2026-07-02 (session 2): levers #1 + #2 done — bottleneck is per-solve stiffness

Both quick/structural runtime levers are implemented and the full dynamics suite is
green (Zeno 14/14, slice-3 parity 112/112, culverts end-to-end 6/6).  **Changes are
NOT yet committed** (working tree has the three `src/` edits below).

- **Lever #1 — decouple ODE tolerance (`solveDynNetwork!`).** Added `ode_abstol`
  (1e-8) / `ode_reltol` (**1e-7**) kwargs passed to `solve(...)`; the topology
  `abstol` stays pinned at 1e-8 (noise-floor / `-abstol` dead-bands).  Measured the
  safe headroom: `ode_reltol=1e-4` and `1e-6` both **break the three-state contract**
  on the pinned culvert outlet (trap 13) and blow full-coverage parity past 1e-6;
  `1e-7` is the loosest safe value (full-cov parity 1.3e-8, culvert scenario runs
  clean).  The RESUME "20% faster" figure below came from the *unsafe* 1e-4.  Net win
  from #1 alone: 40.6 s → 37.3 s (~8%).  Tight ODE tol is load-bearing for the
  contract on balanced culvert traps — do not loosen further.

- **Lever #2 — touch-gating (`fill_sequence` main loop, plan D4/§8).** The call site
  now gates entry to `_touch_networks!`: it runs only when some network fired a member
  (`u.index ∈ ctx.global_ix`) or its inflow changed (`inflow_sources ∩
  getinflowupdates`).  Verified safe wholesale-skip: network growth/merge boundaries
  are *members* (the terminal unfilled trap is a network node in `_subnetwork`), so
  `fired` captures every topology change.  Untouched networks keep their cached state /
  prediction / changetimeest; `_network_amount_updates` is now also gated on
  `network_touched` (untouched state is not committed to `cur_time`; amounts interpolate
  per §9).  The per-context (partial-rebuild) gate is still deferred — `@@@` note moved
  into `_touch_networks!`.

- **The finding:** on the culvert-drought mini scenario, gating **skips 417/512
  events** and cuts `solveDynNetwork!` calls **1154 → 318** (3.6×) — but runtime only
  40.6 → 36.7 s.  The skipped solves were the *cheap* steady-drain ones; the ~318
  remaining are the **stiff culvert-active** solves (~116 ms each now vs 36 ms avg
  before).  **Bottleneck = per-solve cost, not solve count.**  So the real remaining
  win is **lever #3 (stiff / auto-switching solver)**, which needs an AD-clean rate
  function or a finite-difference Jacobian (blocked today by explicit `Float64(...)`
  casts in `water_level`/`_surface_level`).  Levers #1/#2 are still worth keeping:
  #2 is a genuine 3.6× solve-count cut that will dominate in multi-network / less-stiff
  scenarios and bounds worst case.

**Next step:** lever #3.  Make `dynNetworkRateFunction!` (and the `water_level` /
`_surface_level` / interpolation path it calls) AD-clean, then try
`AutoTsit5(Rosenbrock23())` or `Rodas5P()` with an `ADTypes` Jacobian; re-measure the
318 stiff solves.  See "Runtime levers, ranked" §3 below.

## TL;DR

- The **culvert-drought Zeno-chatter hang is FIXED and committed** (`945430a`). Full
  dynamics suite green (slice-3 parity 112/112, leg-2 culverts 6/6).
- The open item is **RUNTIME**: a culvert network draining during drought is slow
  (~42–54 s warm for 514 events on the mini grid; the *same network without a
  culvert* is 0.5 s). Correctness is done; this is pure performance.
- **Next step I recommended:** decouple the ODE integration tolerance from the
  topology `abstol`, loosen the ODE side (`1e-8 → 1e-4`), re-run the suite to confirm
  parity holds, measure speedup. Then touch-gating (the big structural win).

## Where things stand (commit history on `SUrbArea`)

```
945430a Fix Zeno chatter of a culvert outlet pinned at its spillpoint (drought)  <-- latest
a24b4a7 Add end-to-end culvert tests through fill_sequence (validation leg 2)
cbf8efc Add §10 weather-boundary finalization for network traps
4393aaf Fix drought drain-rate formula and stale network-absorption state
85faa7e Add agent reports, source prompt, and verification image ...
```

Working tree: clean except two emacs temp files under `agent/prompts/` (`#...#`,
`.#...`) — untracked, leave them. `examples/Project.toml` is skip-worktree, never
commit it.

Integration plan status (see `agent/reports/integrate_networks_plan.md`, 5-point
order): points 1–4 done; point 5 (validation) legs 1+2 done, leg 3 (animation
script) not done. Also see agent memory `dynnetwork-integration-status.md`.

## The chatter fix (already landed — context so I don't re-diagnose)

**Symptom:** `fill_sequence(ts, w; culverts=[cv])` hung when the culvert network had
to DRAIN during drought. Repro scenario (mini grid):
```julia
w   = [WeatherEvent(0.0, 1.0), WeatherEvent(0.1, 0.0)]
inf = fill(0.02, size(ts.topography))
cv  = DynCulvert(CartesianIndex(7,119), CartesianIndex(199,4), 0.5,0.6,0.5,0.02,1.7)  # inlet trap 233, outlet trap 13
fill_sequence(ts, w; infiltration=inf, culverts=[cv])
```

**Root cause:** trap 13 (culvert OUTLET) sits at its spillpoint with culvert supply
≈ its infiltration, so its net rate is at the numerical noise floor (~1e-8). The ODE
stepper's ULP-level rounding walked its volume back and forth across the capacity
boundary, firing `:fill`/`:unspill` at picosecond spacing (Zeno). No solver-level
threshold fixes it — the event root-finder locates the crossing regardless of step
size or rate clamp (five failed attempts documented in the memory file).

**Fix (two coordinated pieces, keyed on "ns resolution is meaningless; ms is plenty"):**
1. *Noise-floor steady state* (`networksolver.jl`): the t=0 `:unspill`/`:empty` fast
   paths and the in-integration `:unspill` root fire only past `-abstol` (not `<0`),
   matching the steady-state callback. A trap balanced within `abstol` stays full.
2. *Localized event-time slack* (`network_context.jl`, `LOCAL_SLACK = 10 ms`): when a
   rebuilt network predicts an event that REVERSES the transition its trap just fired,
   within `LOCAL_SLACK`, march the network forward one `LOCAL_SLACK` with topology
   detection off (new `solveDynNetwork!(...; topology_events=false)` blind-advance
   mode), letting the pinned trap drift physically off the boundary, then re-predict.
   Keyed on *reversal-within-slack* (not raw spacing) so genuine sub-ms fills are
   untouched → parity preserved. `_reverses`, `_resolve_zeno_chatter!` are the helpers.
   `LOCAL_SLACK` is flagged `@@@` to fold into the global `time_slack` later.

Tests: `@testset "Zeno guard: noise-floor steady-state and blind-advance"` in
`test/dynamics_test.jl` (`_reverses`, the `-abstol` dead-band, the blind-advance mode).

## RUNTIME — evidence already gathered (don't re-measure)

Culvert drought, mini grid, warm second run:

| Metric | Value |
|---|---|
| Events | 514 |
| `solveDynNetwork!` calls | **1154** (~2.2 / event) |
| Runtime (warm) | ~42–54 s |
| Same network, NO culvert | 0.5 s |
| Blind-advances (Zeno marches) | ~60 (NOT the bottleneck) |
| Per-solve cost | ~36 ms for a **2-trap** ODE (pathological → stiffness) |

Two multiplicative problems: **too many solves** (most redundant) and **each solve
too expensive** (culvert makes the 233↔13 pair stiff).

Quick experiments run (all reverted):
- `reltol 1e-8 → 1e-4`: ~20% faster. Modest — tolerance is not the main cost.
- Stiff solver `Rosenbrock23()`: FAILS — rate function isn't ForwardDiff-clean
  (explicit `Float64(...)` in `water_level`/`_surface_level`; the interpolation).
  Would need AD-clean rate fn or `AutoFiniteDiff()` Jacobian (note: `autodiff=false`
  is deprecated; use an `ADTypes` specifier).
- `dtmin` / minimum ODE step: NOT useful (root-finder circumvents it; risks accuracy).

## Runtime levers, ranked (my recommendation)

1. **Touch-gating (biggest win, ~10× fewer solves).** The culvert net is solved 1154×
   but changes on only a handful of events; during drought its *external* inflow is
   constant so its cached prediction stays valid until a member fires or inflow
   changes. This is plan Decision **D4**, deferred (`@@@` in `_touch_networks!`) due
   to context-merge complexity for *growing* networks — simpler for draining. Cut
   1154 → tens of solves.
2. **Decouple + loosen ODE tolerance (easiest quick win, ~20%).** `abstol` now doubles
   as the topology noise-floor threshold, so DON'T just loosen `abstol`. Add separate
   `ode_abstol`/`ode_reltol` (default ~1e-4) passed to `solve(...)`, keeping the
   topology `abstol` (=noise floor / `-abstol` guards) at 1e-8. Re-run full suite —
   parity is 1e-6, so loosening the ODE side must be validated not to shift fill times.
3. **Stiff / auto-switching solver (medium, big per-solve win, compounds with #1).**
   Requires making the rate fn AD-clean OR a finite-difference Jacobian. 36 ms/solve
   for 2 traps says the explicit default is taking thousands of tiny stiff steps.
4. **Cache `DynNetworkRateParams` across repeated solves of an unchanged net** (commit
   + predict + blind-advances all rebuild params via `_build_rate_params`; invalidate
   only on rebuild).
5. **Cheapen `culvert_rate`** — only if a profile shows per-*eval* cost dominates over
   step count (probably it doesn't). Defer.

Start order: **#2 (validate) → #1 (structural) → #3 (compounds)**.

## Key files / functions

- `src/dynamics/networksolver.jl`: `solveDynNetwork!` (fast paths, `topology_events`
  blind-advance branch, the `solve(...)` call ~line 1102, tolerances), `dynNetworkRateFunction!`,
  `_build_event_callback` (`:unspill` root shift), `_build_rate_params`, `culvert_rate`.
- `src/dynamics/network_context.jl`: `LOCAL_SLACK`/`ZENO_MAX_STEPS`, `_reverses`,
  `_resolve_zeno_chatter!`, `_touch_networks!` (the `@@@` touch-gating deferral note is
  here — that's where #1 lands), `_predict_network!`, `_commit_network!`.
- `test/dynamics_test.jl`: "Zeno guard" testset; "fill_sequence dyn_traps parity
  (slice 3)"; "culverts through fill_sequence" (leg 2).

## How to run (P-core pinned)

```
exec taskset -c 0,2,4,6,8,10,12,14 julia --threads=8 --project=. -e '...'
```
Dynamics tests in isolation (ignore MaskingTests in runtests.jl):
```julia
using Test, SurfaceWaterIntegratedModeling, LazyArtifacts
const SWIM = SurfaceWaterIntegratedModeling
grid = loadgrid(joinpath(artifact"swim_testdata","data","small","mini.txt"))
ts   = spillanalysis(grid, usediags=true)
include("test/dynamics_test.jl")
```
