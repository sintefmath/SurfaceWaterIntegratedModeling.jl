# Visual verification of `setup_network`.
#
# Builds a colour texture from the DynNetwork(s) returned by `setup_network` and
# drapes it on the 3D terrain surface using SWIM's own `plotgrid`/`drape_surface`
# utilities, so the flow-path / trap / merge wiring can be checked against the
# relief:
#
#   * grey                : terrain (shading conveys relief)
#   * coloured blobs      : trap footprints, one hue per network
#   * red lines           : main flow paths (water routes)
#   * orange lines        : tributary paths that merge into another path
#   * green discs         : the start cells handed to setup_network
#
# Run from the repo root:
#   using Pkg; Pkg.activate("examples")
#   include("examples/verification/dynamic_network.jl")
#   fig = verify_dynamic_network([CartesianIndex(7, 119)])

using SurfaceWaterIntegratedModeling
using Pkg.Artifacts
import GLMakie, Images, ColorSchemes

const SWIM = SurfaceWaterIntegratedModeling
const RGBf = Images.RGB{Float64}

# colour constants
const BG_COLOR    = RGBf(0.82, 0.82, 0.82)
const PATH_COLOR  = RGBf(0.85, 0.10, 0.10)   # main flow paths
const MERGE_COLOR = RGBf(1.00, 0.55, 0.00)   # tributary (merging) paths
const START_COLOR = RGBf(0.00, 0.70, 0.00)   # start cells

# load the mini.txt artifact grid and its trap structure
function _mini_trapstructure()
    grid = loadgrid(joinpath(datapath_testdata(), "data", "small", "mini.txt"))
    return spillanalysis(grid, usediags=true)
end

# a vivid, distinct hue per network for the trap footprints
_net_color(i) = _blend(RGBf(ColorSchemes.tab10[mod1(i, 10)]), RGBf(1, 1, 1), 0.10)
_blend(a, b, t) = RGBf((1 - t) * a.r + t * b.r,
                       (1 - t) * a.g + t * b.g,
                       (1 - t) * a.b + t * b.b)

# paint a cell, optionally widening to a (2*width+1)-square so paths read as lines
function _paint!(tex, cell::CartesianIndex, color; width=0)
    for di in -width:width, dj in -width:width
        c = cell + CartesianIndex(di, dj)
        checkbounds(Bool, tex, c) && (tex[c] = color)
    end
end

# indices of paths that are tributaries (registered in some other path's merges)
function _tributary_paths(net)
    trib = Set{Int}()
    for p in net.flow_paths, m in p.merges
        push!(trib, m)
    end
    return trib
end

# build the drape texture (grid-shaped Matrix{RGB}) from the networks.
# Paths are drawn first and the (vivid) trap footprints painted on top, so the
# lakes read as solid blobs and the thin path lines show in the connectors
# between them rather than drowning the lakes.
function _network_texture(ts, nets, starts; path_width=0)
    CI = CartesianIndices(ts.topography)
    tex = fill(BG_COLOR, size(ts.topography))

    # flow paths first, tributaries highlighted
    for net in nets
        trib = _tributary_paths(net)
        for (k, p) in enumerate(net.flow_paths)
            color = k in trib ? MERGE_COLOR : PATH_COLOR
            for cell in p.cells
                _paint!(tex, cell, color; width=path_width)
            end
        end
    end

    # trap footprints on top, one vivid hue per network
    for (i, net) in enumerate(nets)
        col = _net_color(i)
        for t in net.traps
            tex[CI[ts.footprints[t.trap_ix]]] .= col
        end
    end

    # start cells last, so they sit on top of everything (single cell, like paths)
    for s in starts
        _paint!(tex, s, START_COLOR; width=0)
    end
    return tex
end

"""
    verify_dynamic_network(starts; ts=_mini_trapstructure(), full_traps=:all,
                           heightfac=3.0, downsamplefac=1.0, path_width=0)

Run `setup_network(ts, starts, full_traps)` and drape the resulting network(s) on
the 3D terrain. `starts` is a vector of `CartesianIndex` start cells. Returns the
GLMakie figure (the surface and scene are also returned, as `(fig, surf, scene)`).

`full_traps=:all` (default) fills every trap — a valid downward-closed fill — or
pass your own vector of filled trap indices. `heightfac` exaggerates the relief.
`path_width` widens flow-path lines to a (2*path_width+1)-cell band (0 = 1 cell).
"""
function verify_dynamic_network(starts; ts=_mini_trapstructure(), full_traps=:all,
                                heightfac=3.0, downsamplefac=1.0, path_width=0)
    ftraps = full_traps === :all ? collect(1:numtraps(ts)) : full_traps
    nets = setup_network(ts, starts, ftraps)
    @info "setup_network produced $(length(nets)) network(s)" *
          " with $(sum(n -> length(n.traps), nets; init=0)) trap(s)" *
          " and $(sum(n -> length(n.flow_paths), nets; init=0)) flow path(s)"

    tex = _network_texture(ts, nets, starts; path_width=path_width)
    surf, fig, scene = plotgrid(ts.topography; texture=tex,
                                heightfac=heightfac, downsamplefac=downsamplefac)
    return fig, surf, scene
end

# Ready-made scenario when run as a script (mirrors the Tier-3 "long" case).
if abspath(PROGRAM_FILE) == @__FILE__
    fig, = verify_dynamic_network([CartesianIndex(7, 119)], heightfac=0.2)
    #fig, = verify_dynamic_network([CartesianIndex(6,66), CartesianIndex(6, 162)], heightfac=0.2)
    #fig, = verify_dynamic_network([CartesianIndex(6, 6), CartesianIndex(6, 102), CartesianIndex(6,162)], heightfac=0.2)
    GLMakie.display(fig)
    GLMakie.wait(GLMakie.Screen())  # keep the window open
end
