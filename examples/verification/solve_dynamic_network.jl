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
#   julia --project=examples examples/verification/solve_dynamic_network.jl pair
# or from the REPL (with the examples environment active):
#   include("examples/verification/solve_dynamic_network.jl")
#   fig = verify_solver()                      # default :long scenario
#   fig = verify_solver(:pair)                 # branched (two tributaries)
#   fig = verify_solver(:long, infiltration=0.3)  # partial-fill / stagnation

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

# Named scenarios (start cells on mini.txt, all traps assumed full for topology).
# Branched scenarios (:pair, :triple) produce a single merged network with tributaries.
const SCENARIOS = Dict(
    :long     => [CartesianIndex(7, 119)],                                          # single chain
    :pair     => [CartesianIndex(6, 66),  CartesianIndex(6, 162)],                  # two tributaries
    :triple   => [CartesianIndex(6, 6),   CartesianIndex(6, 102), CartesianIndex(6, 162)],  # three
    :mixed    => [CartesianIndex(38, 182), CartesianIndex(6, 38)],                  # lake + slope start
)

"""
    verify_solver(scenario=:long; inflow=1.0, infiltration=0.0)

Build the named `scenario` network on mini.txt, run the fill cascade, and plot
the fill-fraction time series.  Returns the GLMakie figure.

Branched scenarios (`:pair`, `:triple`) produce a merged network with tributaries;
the plot colour-codes traps by their index in the topologically-sorted network
(upstream=blue, downstream=red) so tributary branches show as interleaved colours.

With `infiltration=0` every trap fills completely.  Set `infiltration` high enough
relative to `inflow` for some traps to stagnate below full.
"""
function verify_solver(scenario::Symbol=:long; inflow=1.0, infiltration=0.0)
    ts    = _mini_trapstructure()
    starts = SCENARIOS[scenario]
    nets  = setup_network(ts, starts, collect(1:numtraps(ts)))
    net   = only(nets)   # all built-in scenarios produce a single connected network
    nt    = length(net.traps)
    np    = length(net.flow_paths)

    @info "Scenario :$scenario — $nt traps, $np flow paths"

    geom = SWIM._build_trap_geometry(ts, net, fill(infiltration, size(ts.topography)))
    caps = [g.capacity for g in geom]

    times, fracs, events = run_cascade(ts, net; inflow=inflow, infiltration=infiltration)

    nfill  = count(e -> e.kind == :fill,  events)
    nempty = count(e -> e.kind == :empty, events)
    @info "Simulation complete: $nfill fill, $nempty empty events" *
          "  t_final = $(round(times[end]; digits=3))"

    fig = GLMakie.Figure(size=(960, 850))
    cmap = GLMakie.cgrad(:RdYlBu, nt; categorical=false, rev=true)

    # --- top panel: fill-fraction time series ---
    ax_ts = GLMakie.Axis(fig[1, 1];
        title   = "Fill cascade — :$scenario   $nt traps   " *
                  "inflow=$inflow/trap   infiltration=$infiltration/cell",
        xlabel  = "cumulative time",
        ylabel  = "fill fraction  V / capacity",
        limits  = (nothing, (-0.04, 1.12)))

    for i in 1:nt
        GLMakie.lines!(ax_ts, times, fracs[i, :]; color=cmap[i], linewidth=1.2)
    end

    fill_ts = [e.t for e in events if e.kind == :fill]
    isempty(fill_ts) ||
        GLMakie.vlines!(ax_ts, fill_ts; color=(:black, 0.12), linewidth=0.7)

    GLMakie.Colorbar(fig[1, 2]; colormap=GLMakie.cgrad(:RdYlBu; rev=true),
                     limits=(1, nt),
                     label="trap index  (upstream → downstream)",
                     flipaxis=false)

    # --- bottom panel: capacity vs fill time, connected in fill order ---
    fill_events = filter(e -> e.kind == :fill, events)
    sc_li   = [findfirst(t -> t.trap_ix == e.trap, net.traps) for e in fill_events]
    sc_cap  = [caps[li] for li in sc_li]
    sc_time = [e.t for e in fill_events]

    ax_sc = GLMakie.Axis(fig[2, 1];
        title   = "Capacity vs fill time  (connected in fill order)",
        xlabel  = "trap capacity  (volume)",
        ylabel  = "fill time  (cumulative)")

    # connecting line in fill order (grey), dots coloured by topological position
    for k in 1:length(fill_events)-1
        GLMakie.lines!(ax_sc, sc_cap[k:k+1], sc_time[k:k+1];
                       color=(:grey, 0.4), linewidth=1.0)
    end
    GLMakie.scatter!(ax_sc, sc_cap, sc_time;
                     color = sc_li, colormap = GLMakie.cgrad(:RdYlBu; rev=true),
                     colorrange = (1, nt),
                     markersize = 8, strokewidth = 0.5,
                     strokecolor = :black)

    GLMakie.Colorbar(fig[2, 2]; colormap=GLMakie.cgrad(:RdYlBu; rev=true),
                     limits=(1, nt),
                     label="trap index  (upstream → downstream)",
                     flipaxis=false)

    return fig
end

# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    name = isempty(ARGS) ? :long : Symbol(ARGS[1])
    haskey(SCENARIOS, name) ||
        error("unknown scenario \"$name\"; choose one of: " *
              join(sort(string.(keys(SCENARIOS))), ", "))
    fig = verify_solver(name)
    GLMakie.display(fig)
    GLMakie.wait(GLMakie.Screen())
end
