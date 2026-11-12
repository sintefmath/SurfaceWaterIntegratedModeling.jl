# The NBS module: an overview for SWIM developers

This note is for readers who know SWIM's static analysis and `fill_sequence`, but not the Nature-Based Solutions (NBS) addition. It explains the main idea and follows one simulation from the entry point to the results.

Companion documents:
- `NBS_INTRO_DETAILED.md`: the same tour, function by function.
- `NBS_INVARIANTS.md`: the authoritative model, code map and invariants.

---

## 1. What an NBS is in SWIM

An NBS placement is a green roof, a rain bed, a permeable patch or similar. It is described by three things:

- **A footprint:** a set of grid cells.
- **A stack of storage layers** (`NBSSystem`). Water enters the top layer and seeps (*percolates*) downward. Each layer can also *overflow*. The bottom layer's percolation leaves the model to the ground. Three ready-made stacks exist: `puddle`, `elhadiGreenRoof` and `mantillaRRmodel`.
- **Where each layer's overflow goes:**
  - The top `n_terrain` layers (normally 1) overflow **onto the terrain**. The water leaves across the footprint's natural downhill edge and continues as surface flow.
  - Every lower layer that overflows does so through a **pipe** to an outlet cell you specify, for example a green roof's downpipe.

The layers are an ODE in time. That is why an NBS can only live in the part of SWIM that integrates ODEs: the **dynamic networks**, the same machinery that handles culverts.

## 2. The central idea: base SWIM never sees the NBS

The design keeps base SWIM unchanged. `spillanalysis`, `watercourses` and the constant-rate parts of `fill_sequence` treat a footprint as ordinary terrain. Water is routed straight through it. This NBS-free flow field is called the **baseline** (in the code: *oblivious*).

The NBS module then corrects the baseline, but only inside the dynamic network solver:

1. **Throughput.** It computes how much water the baseline sends through the footprint. This is `O_0`: the rain on the footprint plus the flow entering it.
2. **Input.** The NBS receives that water as input, plus any extra water the network currently delivers to it, such as a spilling trap upstream.
3. **Integration.** The layer ODE runs and gives the NBS's actual outflow.
4. **Correction.** Downstream, the network replaces the baseline throughput by the actual outflow. It adds `(actual − O_0)` at the footprint's exits, and adds the pipe flow at the outlets.

Since a downstream trap's inflow already contains `O_0` from the baseline, adding the difference gives it exactly what the NBS really releases.

This difference, `actual − O_0`, can be negative (the NBS is holding water back) or positive (it is releasing stored water). It travels downstream like any other network flow. Cells with spare infiltration absorb added water, and a removal can never push a cell below zero flow. That is the same rule `_track_flow!` uses in base SWIM.

## 3. Workflow, from the user's call to the results

```julia
roof = DynNBSPlacement(elhadiGreenRoof(10, 5, 0.5, 0.5, 1, 1, 0, 1),
                       footprint_cells,                 # Vector{Int}, linear indices
                       [CartesianIndex(80, 121)])       # downpipe outlet

seq    = fill_sequence(tstruct, weather; nbs = [roof], infiltration = infil)
states = all_states_at_timepoints(tstruct, seq, times; nbs = [roof], infiltration = infil)
field  = network_watercourses(tstruct, states[k].filled, Dict{Int,Float64}(),
                              states[k].nbs_layers; nbs = [roof], precipitation = rain,
                              infiltration = infil)
```

### Step 1: `fill_sequence` (entry point)

Before anything else, infiltration is set to zero on every footprint, because the NBS layers own that water. Then, **for each weather period:**

1. **Baseline flow, as usual.** `compute_flow` builds the rates exactly as base SWIM does. It knows nothing about the NBS.
2. **Build the dynamic networks** (`build_network_driver` → `setup_network`).
   - Each placement is first resolved against the terrain. This finds its *inflow cells* (outside cells draining into it), its *exits* (outside cells it drains to), and its *ponds* (trap bottoms inside it).
   - From every exit, outlet and pond, SWIM traces the flow downhill, trap to trap. That builds a network of flow paths and traps.
   - Paths that cross an inflow cell are tied to the NBS too.
   - Everything an NBS touches ends up in one **network component**, solved as one coupled ODE.
3. **Run the event loop.** This is base SWIM's loop: find the next trap that fills or empties, update the spill graph and the rates, record a `SpillEvent`.
   - Traps inside a network are not estimated by base SWIM. Their next event is *predicted* by integrating the network's ODE.
   - When an event concerns a network (one of its traps filled, or its inflow changed), the network is advanced to that time, its membership is updated (it may grow, shrink or merge), and it is rebuilt with the new rates. Then its next event is predicted.
4. **End of the period.** Each network is advanced to the end. The NBS layer storage is saved under the placement's `id`, and the next weather period continues from it.

**What `seq` contains.** It holds the trap events and volumes, including those of network traps. It does not hold the NBS layer volumes; those are recovered in the next step.

### Step 2: `all_states_at_timepoints` / `network_states_at_timepoints`

These replay only the networks forward through `seq` and sample them at the requested times. You get the trap volumes and, per placement `id`, the water in each layer.

### Step 3: `network_watercourses`

This builds the flow-intensity field for one instant, which is what the webapp draws:
1. Start from the baseline `watercourses` field.
2. Apply the culvert diversions and each NBS's correction, using the same throughput, exit shares and downstream rule as the solver.
3. Mark footprint cells `NaN`: surface flow is undefined there, since the layers hold the water.

`interpolate_timeseries(...; nbs = ...)` paints the footprints as their own class on the trap/river map.

## 4. What happens inside the solver, per time step

```
layer storage ──► NBS outflow ──► (outflow − O_0) at exits, pipe flow at outlets
                                          │
          routing pass over the network ◄─┘   (trap spills, culverts, NBS corrections)
                                          │
          flow reaching the inflow cells ─┴──► top-layer input = O_0 + that flow
```

The NBS outflow depends only on what the layers currently hold. So the solver can compute it before routing and feed the routed inflow back into the layers' rates, with no circular dependency.

NBS layers never trigger events of their own. Only traps do (`:fill`, `:empty`, `:unspill`).

## 5. Where to go next

- **Model, code map and invariants:** `NBS_INVARIANTS.md`.
- **Function-level walkthrough:** `NBS_INTRO_DETAILED.md`.
- **Code, in reading order:** `nbs_elements.jl` (layers), `elements.jl` (placement), `nbs_terrain.jl` (coupling to the terrain), `build_network.jl` + `network_utils.jl` (network membership), `networksolver.jl` (the NBS section, `_route_flow`, `dynNetworkRateFunction!`).
