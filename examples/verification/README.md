# Verification scripts

Throwaway, interactive scripts for *visually* checking SWIM internals — the kind
of thing that complements, but does not belong in, the headless test suite.

These are deliberately **not** wired into `docs/make.jl`, so nothing here is
swept into the documentation. They require a display (GLMakie) and are meant to
be run by hand from the REPL while developing.

They share the parent `examples/` project environment:

```julia
using Pkg
Pkg.activate("examples")   # GLMakie, Graphs, SWIM path-dep, ...
include("examples/verification/<script>.jl")
```

| Script | What it checks |
|--------|----------------|
| `dynamic_network.jl` | Overlays the `DynNetwork` produced by `setup_network` on the terrain so the flow-path / trap / merge wiring can be eyeballed. |
| `solve_dynamic_network.jl` | Runs `solveDynNetwork` in a loop (uniform inflow, optional infiltration, all traps start empty) and plots each trap's fill fraction against cumulative time, so the cascade order and fill-rate acceleration can be eyeballed. |
| `culvert_rate.jl` | Sweeps the boundary conditions of two contrasting culverts and plots `culvert_rate` vs inlet head: the top row contrasts free-outfall and submerged-outlet cases (weir→orifice switch, inlet/outlet-control crossover); the bottom row shows the negative reverse-flow branch in the drowned region when `allow_reverse = true`. |
| `culvert_network.jl` | Runs `solveDynNetwork` on the `:long` scenario with and without a culvert that short-circuits ~80 % of the chain, showing the inlet trap draining slower and the outlet trap filling faster. Side-by-side fill-fraction panels with `:unspill` events marked. |
