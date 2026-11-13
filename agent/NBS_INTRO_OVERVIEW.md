# The NBS module: an overview for SWIM developers

This note is for readers who know SWIM's static analysis and `fill_sequence`, but not the Nature-Based Solutions (NBS) addition. It explains the main idea and follows one simulation from the entry point to the results.

Companion documents:
- `NBS_INTRO_DETAILED.md`: the same tour, function by function.
- `NBS_INVARIANTS.md`: the authoritative terms, model and invariants.

This overview uses the same terms as the glossary there.

---

## 1. What an NBS is in SWIM

An **NBS placement** is a green roof, a rain bed, a permeable patch or similar. It is described by three things:

- **A footprint:** the grid cells it covers.
- **A stack of storage layers** (`NBSSystem`). Water enters the top layer.
  - A layer holding more than its threshold `Smax` **overflows**.
  - Above its threshold `Smin`, it **percolates**: it passes water down to the layer below. The bottom layer's percolation leaves the model to the ground.
  - Both rates grow as a power of the excess storage (in mm).
  - Ready-made stacks:
    - `puddle`: one layer that overflows;
    - `elhadiGreenRoof`: a soil layer percolating into a drainage layer that overflows;
    - `mantillaRRmodel`: surface over soil over drainage.
- **Where each layer's overflow goes.**
  - The top `n_terrain` layers (normally 1) overflow **onto the terrain**. The water leaves across the footprint's natural downhill edge, at the **exits**: the outside cells the footprint drains to.
  - Every lower layer that overflows is **piped** to an **outlet** cell you give, for example a green roof's downpipe.

The layers form an ODE in time whose rates depend on what happens elsewhere. Base SWIM's `fill_sequence` works with constant inflow rates and solves each trap on its own (`fill_trap_until`). The only place where several elements are solved together, as one coupled ODE, is the **dynamic network**: a set of traps and flow paths originally built for **culverts** (pipes under roads, whose rate depends on the water levels at both ends). NBS placements live in the dynamic network for the same reason.

## 2. The central idea: base SWIM never sees the NBS

The design keeps base SWIM unchanged. `spillanalysis`, `watercourses` and the constant-rate parts of `fill_sequence` treat a footprint as ordinary terrain, and route water straight through it. This NBS-free result is the **baseline** (in the code: *oblivious*).

As in base SWIM's runoff grid, a cell's **runoff value** is positive surface flow, or, if negative, unused infiltration capacity.

The NBS module corrects the baseline, inside the dynamic network solver only:

1. **Throughput.** It computes what the baseline sends through the footprint: `O_0` = rain on the footprint + baseline flow entering it.
2. **Input.** The top layer receives `O_0`, plus the **draw**: any change the network makes to the flow entering the footprint. An example is the spill of a full trap upstream that belongs to the network (base SWIM does not route such spills, see below). The draw can be negative if an NBS further upstream is holding water back.
3. **Integration.** The layer ODE gives the NBS's real output:
   - `O_t`, the overflow returning onto the terrain;
   - the **pipe flow** at each outlet.
4. **Correction.** In the baseline, the throughput leaves the footprint through its exits, and, if the footprint contains trap bottoms (**ponds**), partly ends in those. The **share** of an exit or pond is the fraction of `O_0` that went that way; the shares of all exits and ponds together sum to 1. The network adds `(O_t − O_0)·share` at each exit, the same into each pond's trap, and the pipe flow at each outlet.

**Why only the difference is added.**
- In the network solver, each trap starts from its **`external_inflow`**: base SWIM's inflow to that trap, read from the baseline. It belongs to a trap, not to the footprint.
- Because the baseline routed water through the footprint, the downstream trap's `external_inflow` already contains its share of `O_0`, reduced by infiltration on the way.
- The difference `(O_t − O_0)·share` travels the same way, so the trap ends up with `O_t·share` reduced on the way instead: what the NBS really releases.

**How a change travels downstream (the per-cell rule).** Everything the network carries along its paths is a change relative to the baseline: an NBS correction, a culvert's flow, or the spill of a full network trap. (Base SWIM does not route the spill of traps that belong to the network; the network does.) Such a change `d` moves downhill cell by cell, following the **per-cell rule**. At a cell with runoff value `R` it becomes `max(R+d, 0) − max(R, 0)`. In words:
- added water first fills a cell's unused infiltration capacity;
- a removal cannot take more than the flow that is there.

For added water, this is the rule base SWIM uses when it routes a trap's spill.

## 3. Workflow, from the user's call to the results

```julia
# elhadiGreenRoof(Smin_soil, Smin_drain, kINF, kDRAIN, nINF, nDRAIN, EVC, EVS11):
#   soil layer Smin (field capacity) 10 mm, drainage layer Smax 5 mm, percolation coefficient
#   kINF and drainage-overflow coefficient kDRAIN, their exponents nINF and nDRAIN, then the
#   (unused) evapotranspiration parameters
roof = DynNBSPlacement(elhadiGreenRoof(10, 5, 0.5, 0.5, 1, 1, 0, 1),
                       footprint_cells,                 # Vector{Int}: linear grid indices
                       [CartesianIndex(80, 121)])       # downpipe outlet (one per piped layer)

seq    = fill_sequence(tstruct, weather; nbs = [roof], infiltration = infil)
states = all_states_at_timepoints(tstruct, seq, times;   # times: Vector of sample times
                                  nbs = [roof], infiltration = infil)
s      = states[k]                                       # state at times[k]
field  = network_watercourses(tstruct, s.filled,
                              Dict{Int,Float64}(),       # trap volumes; only culverts need them
                              s.nbs_layers;              # layer volumes per placement id
                              nbs = [roof], precipitation = rain,   # rain rate in effect at times[k]
                              infiltration = infil)
```

### Step 1: `fill_sequence` (entry point)

Before anything else, infiltration is set to zero on every footprint, because the layers own that water. Then, **for each weather period:**

1. **Baseline, as usual.** `compute_flow` builds the rates exactly as base SWIM does.
2. **Build the dynamic network.** This is done by a *driver*, the object that keeps the network's state through the period.
   1. Each placement gets an **`id`** (its position in the `nbs` vector) and its link to the terrain:
      - its **inflow cells** (outside cells draining into it);
      - its **exits**;
      - its **ponds**.
   2. From every **seed**, SWIM traces the flow downhill to the next trap, and on through full traps along their spill paths. The seeds are the NBS exits, outlets and ponds, plus culvert ends and any traps you list in `dyn_traps`. The traps reached are the **network traps**: their volumes become ODE states, and base SWIM stops estimating them.
   3. The network is split into independent **components**. Each NBS ends up in one component with everything connected to it: the paths leaving its exits and outlets, the paths crossing its inflow cells, and the traps under its ponds. That way all the water reaching it is computed in the same ODE.
3. **Run the event loop.** This is base SWIM's loop: find the next trap that fills or empties, update the spill graph and rates, and record a `SpillEvent`.
   - A network trap's next event is found by integrating its component's ODE ahead.
   - When an event concerns a component (one of its traps changed status, or the baseline inflow to one of its traps changed), the network is **touched**. For all components, not only the affected one:
     1. their ODEs are advanced to the event time;
     2. the affected component is updated: it may gain traps, lose traps, or merge with another component;
     3. base SWIM's spill graph is told which traps the network now routes;
     4. every component is rebuilt with rates from the current baseline;
     5. each one's next event is predicted again.
4. **End of the period.** Each component is advanced to the end. The layer volumes are stored under each placement's `id`, and the next weather period continues from them.

**What `seq` contains.** It holds the trap events and volumes, including those of network traps. It does not hold the NBS layer volumes.

### Step 2: `all_states_at_timepoints` / `network_states_at_timepoints`

These replay only the dynamic network through `seq` and sample it at the requested times. Per time you get the trap states and, per placement `id`, the water in each layer.

### Step 3: `network_watercourses`

This builds the flow-intensity field for one instant, which is what the web viewer of SUrbArea (an application built on SWIM) draws.
1. Start from the baseline `watercourses` field.
2. Apply each culvert's flow and each NBS's correction, with the same functions and per-cell rule as the solver. One difference: this baseline field also routes the spills of full network traps, which the solver receives as draw instead. So at each exit the field removes what its own baseline delivers there and adds the solver's share of the NBS's outflow, with the shares taken from the solver's baseline.
3. Mark footprint cells `NaN`: surface flow is undefined there, since the layers hold the water.

`interpolate_timeseries(...; nbs = ...)` paints the footprints as their own class on the trap/river map.

## 4. What the solver computes, at every evaluation of the ODE

At each evaluation the solver does a **routing pass**: it walks the component's paths and traps from upstream to downstream and works out what arrives where.

```
layer volumes ──► O_t and pipe flow ──► (O_t − O_0)·share at exits/ponds, pipe flow at outlets
                                                     │
        routing pass over the component ◄────────────┘   (trap spills, culverts, NBS corrections)
                                                     │
        draw = change reaching the inflow cells ─────┴──► top-layer input = O_0 + draw
```

The NBS output depends only on what the layers hold. So the solver computes it before the routing pass, and the draw collected during the pass only affects how fast the layers fill. Nothing is circular.

**Events.** Apart from the stop when everything has settled (a steady state), only traps end an integration:
- `:fill`: a trap becomes full;
- `:empty`: a trap drains enough to expose its sub-traps;
- `:unspill`: a full trap stops spilling.

## 5. Where to go next

- **Terms, model and invariants:** `NBS_INVARIANTS.md`.
- **Function-level walkthrough:** `NBS_INTRO_DETAILED.md`.
- **Code, in reading order:** `nbs_elements.jl` (layers), `elements.jl` (placement), `nbs_terrain.jl` (link to the terrain), `build_network.jl` + `network_utils.jl` (network building), `networksolver.jl` (the NBS section, `_route_flow`, `dynNetworkRateFunction!`).
