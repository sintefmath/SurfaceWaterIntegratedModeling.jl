# Steady-state callback design — ongoing discussion

---

## Round 1

### (1) Trap 13 and the test change

**Claude (initial):** Trap 13 has two cells at elevations 34.34751 and 34.41824; the
spillpoint is 34.41824.  Cell 2's bottom equals the spillpoint, so cell 2 only becomes
wet at exactly V = C.  With infil = 0.3/cell and inflow = 0.5 there is no sub-capacity
equilibrium: dV = +0.2 for all 0 < V < C, then :fill fires at V = C.  The old test
expected :none because the sign-flip stagnation caught the rate flip AT capacity.
Switching to trap 57 corrects a wrong test, not a real problem.

**User:** Thanks for the clarification.  I still don't understand why switching to trap 57
was necessary though, if the correct result for trap 13 (:fill) was returned.  Or was the
test looking for something else, and trap 13 didn't fit what the test was needing?

**Claude:** The test had two jobs:

1. Verify that :fill is NOT returned when there is a genuine sub-capacity equilibrium.
2. Verify that :none IS returned when the trap settles below capacity.

Trap 13 with infil = 0.3/cell correctly returns :fill (no sub-capacity equilibrium), so
it can serve job 1 but not job 2.  The stagnation test was specifically written to cover
job 2 — a trap that SETTLES below capacity rather than filling.  For that we need a trap
where the equilibrium actually is below capacity.  Trap 57 (11-cell footprint, inflow 0.3,
infil 0.1/cell → equilibrium at 3 cells wet ≈ 27 % of capacity) provides this.  The
switch was necessary because trap 13 simply does not exercise the code path the test was
meant to cover.

**User (round 2):** 1 sounds fine — make sure the solution to 2 is properly tested.

---

### (2) The `:unspill`-condition fix

**Claude (initial):** My fix: only register `:unspill` VCB conditions for full traps with
strictly positive initial net inflow.  Traps at zero prevent pathological rootfinding.

**User:** Wouldn't your solution mean that traps with initial net inflow zero from the
beginning will never be able to unspill, even if water level starts shrinking at some
point during integration?

**Claude:** You are right, that is a real bug.  Consider a full trap A that starts with
zero net inflow because its upstream feeder (trap B) is still filling.  Once B fills and
starts spilling into A, A's net inflow becomes positive — and if it later drops back
below zero, an `:unspill` event for A must fire.  My filter would have dropped A's
condition entirely, causing the event to be missed.

**Correct fix:** revert the filter; instead add `rootfind = LeftRootFind` to the
`VectorContinuousCallback` constructor.

With `LeftRootFind` the solver only looks for *positive → negative* (downward) crossings.
A condition permanently at 0 is not a downward crossing (0 → 0 does not change sign from
positive to negative), so the root-finder leaves it alone.  This eliminates the
pathological step-shrinking without discarding any conditions.  All three condition kinds
want downward crossings anyway:

| Condition | Expression | Starts at | Fires when |
|-----------|-----------|-----------|------------|
| `:fill`   | C − V     | C (positive) | V reaches C (crosses 0 downward) |
| `:empty`  | V         | V₀ (positive) | V reaches 0 (crosses 0 downward) |
| `:unspill`| inflow − infil | ≥ 0  | net inflow goes negative (crosses 0 downward) |

**User (round 2):** 2 sounds fine — make sure the solution is properly tested.

---

### (3) Stagnation criterion

**Claude (initial):** I proposed a `DiscreteCallback` that fires after `n_consecutive = 3`
successive steps with `|Δu| < tol`.  This number is ad-hoc.

**User (round 1):** Imposing a magic number of tiny steps sounds questionable.  Didn't we
agree that trying to detect stagnation was a red herring, since a trap could switch from
filling to emptying (or vice versa) without actually stagnating?  The only reasonable
stagnation measure I can think of is a global test of current rates: if the maximum (in
absolute value) rate across the net inflows among *all* traps is smaller than some chosen
lower threshold, the system as a whole has effectively stagnated and we can stop
integration based on some kind of "system stagnation" event.

**Claude (round 1):** Agreed.  This is the VCB max-rate condition.  Put `:steadystate` as
a single entry in `conds`, computed as `max(|dV[e]|) − abstol` over all evolving traps:

```julia
elseif ec.kind == :steadystate
    max_rate = 0.0
    for e in evolving           # ← aggregation loop, NOT a per-trap condition
        spilling[e] && (max_rate = abstol + 1.0; break)
        r = abs(dubuf[e])
        r > max_rate && (max_rate = r)
    end
    out[k] = max_rate - abstol  # ← ONE scalar output slot
```

The `for e in evolving` inside the branch is computing the global maximum; there is still
only ONE `:steadystate` entry in `out`.  Add `rootfind = LeftRootFind` to the VCB
constructor (this also fixes issue 2).  Keep the spilling-veto (forces `out[k]` positive
when any trap is spilling so `:fill` always wins the race).  Add a `sol.retcode` check
for the rare Zeno/maxiters case.

**User (round 2):** I am still confused about your proposed strategy.  Stagnation is a
*global* phenomenon, so why attach a condition to each specific evolving trap and test it
inside a loop in the function 'condition'?  Shouldn't stagnation rather be identified as
a collective measure, comparing all rates in the system and determining if they are *all*
negligible?

**Claude (round 2):** The architecture in my round-1 pseudocode does exactly that — the
single `out[k]` entry is the collective measure — but the way the code is structured makes
it look per-trap.  Let me separate the two concerns more clearly.

There are two different callbacks in play:

**Per-trap VCB** (the existing one) — one output slot per topological event:

| Slot | Kind | Expression |
|------|------|-----------|
| 1 | `:fill` for trap e₁ | C[e₁] − V[e₁] |
| 2 | `:fill` for trap e₂ | C[e₂] − V[e₂] |
| … | … | … |
| k | `:unspill` for full trap f₁ | inflow[f₁] − infil[f₁] |
| … | … | … |

**Separate scalar `ContinuousCallback`** — one global slot for stagnation:

```julia
cb_ss = ContinuousCallback(
    (u, t, integrator) -> begin
        compute dubuf for all traps (one routing pass)
        any(spilling[e] for e in evolving) && return abstol + 1.0  # veto
        maximum(abs(dubuf[e]) for e in evolving) - abstol           # global max
    end,
    integrator -> (event.kind = :steadystate; terminate!(integrator));
    rootfind = LeftRootFind
)

sol = solve(...; callback = CallbackSet(cb_topo, cb_ss))
```

This is architecturally cleaner: the two concerns live in two separate callbacks.  The
scalar `ContinuousCallback` is explicitly global — it has no trap argument, its condition
returns one number, and its affect sets `event.kind = :steadystate`.

The `CallbackSet` combining a VCB and a scalar CC is fine; it is only mixing *continuous*
callbacks (unlike the previous attempt which mixed a VCB and a `DiscreteCallback`).

**One remaining concern:** a transient near-zero crossing of the collective rate (not a
true equilibrium) would also fire `:steadystate`.  For the current implementation (simple
networks, no NBS dynamics that could make rates transiently zero), this seems unlikely in
practice.  If it becomes an issue we can add hysteresis (require the condition to stay
negative for a minimum solver-time window) but that adds complexity; deferring for now.

**User (round 3):** Regarding 3, it's not just that the way your code was structured
"makes it look per trap".  You are nesting a loop inside a loop.  The inner loop loops
over all evolving traps and computes a max.  But you recompute a max every time you
encounter an ec.kind=:steadystate in the outer loop, and ec is tied to a specific trap.
Or am I misunderstanding something?  Are you now suggesting to split this up in two
callbacks instead?  If so, I think that is fine.  But I might need a bit more
clarification on that veto case.  How can a trap be registered as spilling and at the
same time part of the evolving set?

**Claude (round 3):** You are right to push back.  To be precise about the round-1 code:
there would have been exactly ONE `:steadystate` entry in `conds` (added once at the end
of `_event_conditions`, with `ec.trap = 0`), so the inner loop would have run only once
per ODE evaluation.  But your objection stands on a different level: mixing a global
condition into a list whose other entries are all per-trap is inconsistent and confusing.
The two-callback split (agreed) eliminates that entirely.

**On the veto — how can an evolving trap be spilling?**

`evolving` records which traps started *below* capacity at the beginning of this
integration period.  During integration the ODE drives their volume upward.  The routing
function marks a trap as `spilling` whenever `V[e] ≥ geom[e].capacity` at the *current
evaluation point*.

The ODE solver evaluates the condition functions at multiple candidate states *within* a
step before accepting it.  If a candidate state happens to have `V[e]` slightly above C
(an overshoot before the rootfinder locates the exact `:fill` time), the routing function
will momentarily see `spilling[e] = true` for that evolving trap.  At that same candidate
state, `dubuf[e] = 0` (the rate formula zeroes out a spilling trap's change), which makes
`max_rate = 0` and the global condition `= −abstol < 0`.

Without the veto, the steadystate callback would see a downward crossing at this candidate
state and attempt to fire *before* `:fill` does.  The veto short-circuits that: if any
evolving trap shows `V ≥ C`, return `abstol + 1` (positive), keeping the steadystate
condition above zero so `:fill` wins.

In summary: the veto is not about a trap being both "in the network's permanent full set"
and "evolving" at the same time — that cannot happen by construction.  It is purely about
transient ODE overshoot during rootfinding within a step.

---

## Proposed revised design

1. **`_event_conditions`**: restore unconditional `:unspill` registration (revert the
   `inflow0` filter).  Remove the `inflow0` parameter from the signature.

2. **`_build_event_callback`**: add `rootfind = LeftRootFind` to the VCB constructor.
   Remove the `:steadystate` entry from `conds` (it moves to a separate callback).

3. **New `_build_steadystate_callback`** (replaces the `DiscreteCallback` version):
   returns a scalar `ContinuousCallback` with the global max-rate condition and the
   spilling-veto.

4. **`solveDynNetwork!`**: `CallbackSet(cb_topo, cb_ss)` as before, but now both
   members are continuous callbacks.  Add `sol.retcode` check for Zeno/maxiters.

5. **Tests**: add a test that specifically exercises a full trap with initial zero net
   inflow that later receives positive inflow and then must unspill (verifying the
   `LeftRootFind` fix for issue 2).

Shall I implement this?

---

## Implementation note (divergence from proposed design)

Item 3 above was revised during implementation.  The proposed `ContinuousCallback`
for `cb_ss` turned out to be wrong for this system.

**The problem:** `ContinuousCallback` evaluates its condition function at *dense-output*
(polynomial interpolant) states within each ODE step, not only at accepted endpoints.
`wetted_infiltration` is a step function of V — it jumps at cell-volume boundaries.
When the ODE polynomial interpolant evaluates V at a point that straddles a boundary
mid-step, the condition function sees a discontinuous jump and interprets it as a
downward crossing.  In the two-trap steady-state test, the CC fired when trap 57 was
still filling at rate 0.1 (far from equilibrium), because an interpolated state briefly
had V[57] above the 7-cell equilibrium boundary while V[313] was already settled.

**The fix:** `cb_ss` remains a `DiscreteCallback`, but the `n_consecutive` counter and
`reltol` magic parameters are removed.  The new condition checks the global maximum
`|dV/dt|` at each accepted step endpoint (where the ODE state is genuine, not
interpolated), and fires on the first step where that maximum is below `abstol`.
The spilling-veto (return false if any evolving trap has V ≥ capacity) replaces the
old consecutive-counter approach for handling the race with the VCB `:fill` event.

**Final `_build_steadystate_callback` signature:**
```julia
function _build_steadystate_callback(p::DynNetworkRateParams,
                                     evolving::AbstractVector{<:Integer},
                                     abstol::Real)
    function condition(u, t, integrator)
        inflow, spilling = _routed_inflow(u, p)
        for e in evolving
            spilling[e] && return false   # veto: :fill wins
        end
        for e in evolving
            abs(inflow[e] - wetted_infiltration(p.geom[e], u[e])) >= abstol && return false
        end
        return true
    end
    return DiscreteCallback(condition, terminate!; save_positions = (false, false))
end
```

Items 1, 2, 4, and 5 were implemented as proposed (item 4 without the `sol.retcode`
check — not needed since `event.kind == :none` already handles the maxiters case).
