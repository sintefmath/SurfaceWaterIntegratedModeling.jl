# Visual verification of `solveDynNetwork`.
#
# Runs a fill-cascade scenario on the mini.txt terrain:
#   * constant uniform inflow (one unit per trap per time unit)
#   * zero infiltration
#   * all traps start empty
#
# `solveDynNetwork` is called in a loop, each call advancing to the next
# topology-changing event (:fill or :empty), until the network reaches
# steady state (:none).
#
# The result is plotted as a fill-fraction time series: one line per trap,
# coloured upstream (blue) to downstream (red), with faint vertical markers
# at each fill event.  Correct cascade behaviour looks like a fan of lines
# that rise steeply and then flatten at 1, with downstream traps (red) filling
# noticeably faster than upstream ones once the spill chain reaches them.
#
# Run from the repo root:
#   julia --project=examples examples/verification/solve_dynamic_network.jl
# or from the REPL (with the examples environment active):
#   include("examples/verification/solve_dynamic_network.jl")
#   fig = verify_solver()                      # default: inflow=1.0, no infiltration
#   fig = verify_solver(infiltration=0.3)      # partial-fill / stagnation scenario

using SurfaceWaterIntegratedModeling, LazyArtifacts
import GLMakie

const SWIM = SurfaceWaterIntegratedModeling

function _mini_trapstructure()
    grid = loadgrid(joinpath(datapath_testdata(), "data", "small", "mini.txt"))
    return spillanalysis(grid, usediags=true)
end

# ---------------------------------------------------------------------------

"""
    run_cascade(ts, net; inflow=1.0, infiltration=0.0, max_events=2000)

Evolve `net` forward under constant uniform `inflow` (per trap) and uniform
`infiltration` (per grid cell), starting with all traps empty.  Calls
`solveDynNetwork` in a loop until steady state or `max_events` is reached.

Returns `(times, fracs, events)`:
- `times`  – cumulative event times, length `n+1` (includes t=0)
- `fracs`  – `(ntraps × (n+1))` matrix of fill fractions at each snapshot
- `events` – vector of `(t, kind, trap_ix)` named tuples for non-`:none` events
"""
function run_cascade(ts, net; inflow=1.0, infiltration=0.0, max_events=2000)
    nt   = length(net.traps)
    ifil = fill(infiltration, size(ts.topography))
    qin  = fill(Float64(inflow), nt)

    geom = SWIM._build_trap_geometry(ts, net, ifil)
    caps = [g.capacity for g in geom]

    state = zeros(Float64, nt)
    t_abs = 0.0
    fracs_norm(s) = s ./ max.(caps, 1e-12)

    times     = [0.0]
    snapshots = [fracs_norm(state)]
    events    = NamedTuple{(:t, :kind, :trap), Tuple{Float64, Symbol, Int}}[]

    for _ in 1:max_events
        res = solveDynNetwork(ts, net, ifil, qin, state)

        t_abs += isfinite(res.time) ? res.time : 0.0
        state  = copy(res.state)

        push!(times, t_abs)
        push!(snapshots, fracs_norm(state))
        res.kind != :none && push!(events, (t=t_abs, kind=res.kind, trap=res.trap))

        res.kind == :none && break
    end

    return times, hcat(snapshots...), events
end

# ---------------------------------------------------------------------------

"""
    verify_solver(; inflow=1.0, infiltration=0.0)

Build the long single-chain network on mini.txt, run the fill cascade, and
plot the fill-fraction time series.  Returns the GLMakie figure.

With `infiltration=0` every trap fills completely and the plot shows a fan of
lines converging to 1 (downstream traps, red, fill faster once the cascade
reaches them).  Set `infiltration` high enough that some traps stagnate below
full to see lines plateau below 1.
"""
function verify_solver(; inflow=1.0, infiltration=0.0)
    ts   = _mini_trapstructure()
    nets = setup_network(ts, [CartesianIndex(7, 119)], collect(1:numtraps(ts)))
    net  = only(nets)
    nt   = length(net.traps)

    @info "Network: $nt traps, $(length(net.flow_paths)) flow paths"

    times, fracs, events = run_cascade(ts, net; inflow=inflow, infiltration=infiltration)

    nfill  = count(e -> e.kind == :fill,  events)
    nempty = count(e -> e.kind == :empty, events)
    @info "Simulation complete: $nfill fill, $nempty empty events" *
          "  t_final = $(round(times[end]; digits=3))"

    fig = GLMakie.Figure(size=(960, 500))
    ax  = GLMakie.Axis(fig[1, 1];
        title   = "Fill cascade — $nt traps   " *
                  "inflow=$(inflow)/trap   infiltration=$(infiltration)/cell",
        xlabel  = "cumulative time",
        ylabel  = "fill fraction  V / capacity",
        limits  = (nothing, (-0.04, 1.12)))

    # one line per trap, upstream (blue) -> downstream (red)
    cmap = GLMakie.cgrad(:RdYlBu, nt; categorical=false, rev=true)
    for i in 1:nt
        GLMakie.lines!(ax, times, fracs[i, :];
                       color=cmap[i], linewidth=1.2)
    end

    # faint vertical ticks at fill events
    fill_ts = [e.t for e in events if e.kind == :fill]
    isempty(fill_ts) ||
        GLMakie.vlines!(ax, fill_ts; color=(:black, 0.12), linewidth=0.7)

    GLMakie.Colorbar(fig[1, 2]; colormap=GLMakie.cgrad(:RdYlBu; rev=true),
                     limits=(1, nt),
                     label="trap index  (upstream → downstream)",
                     flipaxis=false)

    return fig
end

# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    fig = verify_solver()
    GLMakie.display(fig)
    GLMakie.wait(GLMakie.Screen())
end
