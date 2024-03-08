# Development over time

## Constructing a sequence of events

The central function for analysing development over time is
[`fill_sequence`](@ref), which computes a sequence of [`SpillEvent`](@ref)s that
describe how water accumulates, flows and drains on a terrain over time,
considering infiltration and potentially changing weather conditions.

```@autodocs
Modules = [SurfaceWaterIntegratedModeling]
Pages = ["fill_sequence.jl"]
```

## Associated structs and functions

```@autodocs
Modules = [SurfaceWaterIntegratedModeling]
Pages = ["weatherevent.jl", "spillevent.jl", "spillgraph.jl", "rateinfo.jl", "flow.jl"]
```
