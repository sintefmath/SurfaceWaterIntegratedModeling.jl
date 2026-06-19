# Visual verification of `setup_network`.
#
# Loads the `mini.txt` test grid, runs `setup_network` from one or more start
# cells, and overlays the resulting DynNetwork on the terrain so the wiring can
# be checked by eye:
#
#   * grey heatmap            : terrain elevation
#   * pale coloured blobs     : trap footprints (one colour per network)
#   * coloured poly-lines     : flow paths
#   * white dots              : trap bottoms
#   * black arrows            : path -> target_trap  and  trap -> spill_path
#   * dashed arrows           : tributary merges (path -> path)
#   * green stars             : the start cells handed to setup_network
#
# Run from the repo root:
#   using Pkg; Pkg.activate("examples")
#   include("examples/verification/dynamic_network.jl")
#   fig = verify_dynamic_network([CartesianIndex(7, 119)])

using SurfaceWaterIntegratedModeling
using Pkg.Artifacts
import GLMakie, ColorSchemes

const SWIM = SurfaceWaterIntegratedModeling

# load the mini.txt artifact grid and its trap structure
function _mini_trapstructure()
    grid = loadgrid(joinpath(datapath_testdata(), "data", "small", "mini.txt"))
    return spillanalysis(grid, usediags=true)
end

# pick a visually distinct colour per network
_net_color(i, n) = GLMakie.get(ColorSchemes.tab10, (i - 1) / max(n - 1, 1))

# centroid of a trap footprint, as a plottable Point2 (x = row, y = col)
function _trap_point(ts, trap_ix)
    CI = CartesianIndices(ts.topography)
    cells = CI[ts.footprints[trap_ix]]
    GLMakie.Point2f(sum(c[1] for c in cells) / length(cells),
                    sum(c[2] for c in cells) / length(cells))
end

_cell_point(c::CartesianIndex) = GLMakie.Point2f(c[1], c[2])

# draw an arrow from p to q on the given axis
function _arrow!(ax, p, q; color=:black, style=:solid)
    GLMakie.lines!(ax, [p, q]; color, linestyle=style, linewidth=1.5)
    GLMakie.scatter!(ax, [q]; color, marker=:utriangle, markersize=10)
end

# overlay one network's footprints, paths, traps and links
function _draw_network!(ax, ts, net, col)
    CI = CartesianIndices(ts.topography)

    # trap footprints (pale) and trap bottoms (white dots)
    for t in net.traps
        cells = _cell_point.(CI[ts.footprints[t.trap_ix]])
        GLMakie.scatter!(ax, cells; color=(col, 0.20), markersize=4)
        GLMakie.scatter!(ax, [_trap_point(ts, t.trap_ix)];
                         color=:white, strokecolor=col, strokewidth=1, markersize=9)
    end

    # flow paths as poly-lines
    for p in net.flow_paths
        isempty(p.cells) && continue
        GLMakie.lines!(ax, _cell_point.(p.cells); color=col, linewidth=2)
    end

    # path -> target_trap links
    for p in net.flow_paths
        (p.target_trap == 0 || isempty(p.cells)) && continue
        _arrow!(ax, _cell_point(p.cells[end]), _trap_point(ts, net.traps[p.target_trap].trap_ix))
    end

    # trap -> spill_path links
    for t in net.traps
        (t.spill_path == 0 || isempty(net.flow_paths[t.spill_path].cells)) && continue
        _arrow!(ax, _trap_point(ts, t.trap_ix), _cell_point(net.flow_paths[t.spill_path].cells[1]))
    end

    # tributary merges (path -> path), drawn dashed
    for p in net.flow_paths, m in p.merges
        trib = net.flow_paths[m]
        (isempty(trib.cells) || isempty(p.cells)) && continue
        _arrow!(ax, _cell_point(trib.cells[end]), _cell_point(p.cells[1]); style=:dash)
    end
end

"""
    verify_dynamic_network(starts; ts=_mini_trapstructure(), full_traps=:all)

Run `setup_network(ts, starts, full_traps)` and plot the result. `starts` is a
vector of `CartesianIndex` start cells. Returns the GLMakie figure. Pass
`full_traps=:all` (default) to fill every trap (a valid downward-closed fill),
or supply your own vector of filled trap indices.
"""
function verify_dynamic_network(starts; ts=_mini_trapstructure(), full_traps=:all)
    ftraps = full_traps === :all ? collect(1:numtraps(ts)) : full_traps
    nets = setup_network(ts, starts, ftraps)
    @info "setup_network produced $(length(nets)) network(s)" *
          " with $(sum(n -> length(n.traps), nets; init=0)) trap(s)" *
          " and $(sum(n -> length(n.flow_paths), nets; init=0)) flow path(s)"

    fig = GLMakie.Figure(size=(900, 900))
    ax = GLMakie.Axis(fig[1, 1]; aspect=GLMakie.DataAspect(),
                      title="setup_network: $(length(nets)) network(s) from $(length(starts)) start(s)")
    GLMakie.heatmap!(ax, ts.topography; colormap=:greys)

    for (i, net) in enumerate(nets)
        _draw_network!(ax, ts, net, _net_color(i, length(nets)))
    end
    GLMakie.scatter!(ax, _cell_point.(starts);
                     color=:green, marker=:star5, markersize=16, strokecolor=:black, strokewidth=1)
    return fig
end

# A couple of ready-made scenarios mirroring the Tier-3 integration tests, so a
# quick `include` shows something immediately.
if abspath(PROGRAM_FILE) == @__FILE__
    fig = verify_dynamic_network([CartesianIndex(7, 119)])
    GLMakie.display(fig)
    GLMakie.wait(GLMakie.Screen())  # keep the window open when run as a script
end
