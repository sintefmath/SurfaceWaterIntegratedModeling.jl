# Visual verification of `solveDynNetwork` with culverts.
#
# Scenario: the :long chain on mini.txt, starting from cell (7,119).  A single
# concrete culvert (r=0.3 m) is placed from trap 233 — the first trap the
# flow path reaches after the start cell — to trap 13, which lies roughly
# 80 % of the way through the same chain.  The culvert short-circuits a large
# part of the normal cascade: trap 233 is drained while filling (via the
# culvert inlet) and trap 13 is fed ahead of schedule (via the culvert outlet).
#
# Three panels on a single figure:
#   left   (no culvert)  — reference fill-fraction cascade; traps coloured
#                          upstream (blue) → downstream (red).
#   centre (with culvert) — same cascade with the culvert active; the inlet
#                           trap is magenta, the outlet trap is cyan; thin
#                           magenta vlines mark :unspill events.
#   right  — shared colourbar for the topological index.
#
# What correct behaviour looks like:
#   * The inlet trap (magenta) should plateau below 1 or fill noticeably
#     more slowly than its same-coloured counterpart on the left.
#   * The outlet trap (cyan) should fill earlier/faster than in the reference
#     because it receives extra water via the culvert.
#   * Intermediate traps between 233 and 13 may fill a little more slowly
#     (less spill from 233 feeding them via the normal path).
#
# Run from the REPL (needs a display; shares the examples/ project):
#
#   using Pkg
#   Pkg.activate("examples")
#   include("examples/verification/culvert_network.jl")
#   fig = verify_culvert_network()
#
# Or from the repo root:
#   julia --project=examples examples/verification/culvert_network.jl

using SurfaceWaterIntegratedModeling, LazyArtifacts
import GLMakie

const SWIM = SurfaceWaterIntegratedModeling

const _START_CELL = CartesianIndex(7, 119)   # :long scenario start
const _CV_INLET   = CartesianIndex(7, 119)   # culvert inlet  (trap 233)
const _CV_OUTLET  = CartesianIndex(199, 4)   # culvert outlet (trap 13)
const _CV_RADIUS  = 0.5                      # barrel radius, metres
# Trap 233 is only ~0.1 m deep, so its head at full gives Q_culvert ≈ 0.055 m³/s.
# inflow=0.02 per trap gives Q_culvert/inflow ≈ 2.5 at equilibrium, producing a
# clear plateau of the inlet trap at ~50 % of its capacity.

function _mini_trapstructure()
    grid = loadgrid(joinpath(datapath_testdata(), "data", "small", "mini.txt"))
    return spillanalysis(grid, usediags=true)
end

# Run the fill cascade from all-empty to steady state, rebuilding the network
# after every topology-changing event (fill / unspill) as required by the
# solveDynNetwork! caller contract.
#
# Returns:
#   times      — cumulative event times (length n+1, includes t=0)
#   fracs      — (ntraps × (n+1)) fill-fraction matrix, traps in first-network order
#   events     — vector of (t, kind, trap) named tuples for non-:none events
#   gix_order  — global trap_ix in the order used for rows of `fracs`
function _run_cascade(ts, start_cells, cvlts=DynCulvert[];
                      inflow=1.0, infiltration=0.0, max_events=5000)
    ifil       = fill(infiltration, size(ts.topography))
    full_traps = Int[]

    # First network (no full traps yet) establishes the full trap list and capacities.
    first_nets = setup_network(ts, start_cells, full_traps; culverts=cvlts)
    isempty(first_nets) && return [0.0], zeros(0, 1), NamedTuple{(:t,:kind,:trap),Tuple{Float64,Symbol,Int}}[], Int[]
    first_net  = first_nets[1]
    gix_order  = [t.trap_ix for t in first_net.traps]  # stable row order for `fracs`
    ntotal     = length(gix_order)

    # Per-trap capacity (by global trap_ix) and current volume.
    cap_by_gix = Dict{Int,Float64}()
    vol_by_gix = Dict(gix => 0.0 for gix in gix_order)

    function _caps_and_state!(net)
        geom = SWIM._build_trap_geometry(ts, net, ifil)
        state = Vector{Float64}(undef, length(net.traps))
        for (i, trap) in enumerate(net.traps)
            cap_by_gix[trap.trap_ix] = geom[i].capacity
            state[i] = trap.trap_ix in full_traps ? geom[i].capacity :
                                                    get(vol_by_gix, trap.trap_ix, 0.0)
        end
        return state
    end

    snap_fn() = [get(vol_by_gix, gix, 0.0) / max(get(cap_by_gix, gix, 1.0), 1e-12)
                 for gix in gix_order]

    # Populate initial capacities.
    _caps_and_state!(first_net)

    t_abs     = 0.0
    times     = [0.0]
    snapshots = [snap_fn()]
    events    = NamedTuple{(:t, :kind, :trap), Tuple{Float64, Symbol, Int}}[]

    for _ in 1:max_events
        nets = setup_network(ts, start_cells, full_traps; culverts=cvlts)
        isempty(nets) && break
        net   = nets[1]
        qin   = fill(Float64(inflow), length(net.traps))
        state = _caps_and_state!(net)

        res = solveDynNetwork!(state, ts, net, ifil, qin)

        # Update global volumes from the (mutated) state.
        for (i, trap) in enumerate(net.traps)
            vol_by_gix[trap.trap_ix] = state[i]
        end

        t_abs += isfinite(res.time) ? res.time : 0.0
        push!(times, t_abs)
        push!(snapshots, snap_fn())
        res.kind != :none && push!(events, (t=t_abs, kind=res.kind, trap=res.trap))

        if res.kind == :fill
            push!(full_traps, res.trap)
        elseif res.kind == :unspill
            # Caller must set the unspilling trap to just below capacity.
            li = findfirst(t -> t.trap_ix == res.trap, net.traps)
            li !== nothing && (vol_by_gix[res.trap] = prevfloat(cap_by_gix[res.trap]))
            filter!(!=(res.trap), full_traps)
        elseif res.kind == :none
            break
        end
    end

    return times, hcat(snapshots...), events, gix_order
end

"""
    verify_culvert_network(; inflow=1.0, infiltration=0.0) -> GLMakie.Figure

Build the :long scenario on mini.txt with and without a $(2*_CV_RADIUS) m-diameter
culvert from the head-of-chain trap to a trap ~80 % of the way downstream, run
both fill cascades from empty, and plot the fill fractions side by side.
"""
function verify_culvert_network(; inflow=0.02, infiltration=0.0)
    ts = _mini_trapstructure()
    cv = DynCulvert(ts, _CV_INLET, _CV_OUTLET; r=_CV_RADIUS)

    times_bare, fracs_bare, evts_bare, gix_bare = _run_cascade(ts, [_START_CELL];
                                                                inflow, infiltration)
    times_cv,   fracs_cv,   evts_cv,   gix_cv   = _run_cascade(ts, [_START_CELL], [cv];
                                                                inflow, infiltration)

    nt_bare = length(gix_bare)
    nt_cv   = length(gix_cv)

    @info "Without culvert : $nt_bare traps"
    @info "With culvert     : $nt_cv traps, 1 culvert  r=$(_CV_RADIUS) m"

    # For the culvert network, find the local row indices of the inlet / outlet traps.
    # (We know them by their global trap_ix via the culvert endpoints.)
    cv_nets = setup_network(ts, [_START_CELL], Int[]; culverts=[cv])
    cv_net  = isempty(cv_nets) ? nothing : cv_nets[1]
    ti_in   = cv_net === nothing ? nothing :
              findfirst(i -> !isempty(cv_net.traps[i].culvert_inlets),  eachindex(gix_cv))
    ti_out  = cv_net === nothing ? nothing :
              findfirst(i -> !isempty(cv_net.traps[i].culvert_outlets), eachindex(gix_cv))

    nfill_bare   = count(e -> e.kind == :fill,    evts_bare)
    nfill_cv     = count(e -> e.kind == :fill,    evts_cv)
    nunspill_cv  = count(e -> e.kind == :unspill, evts_cv)
    @info "Without culvert : $nfill_bare fill events  t_final=$(round(times_bare[end]; digits=2))"
    @info "With culvert     : $nfill_cv fill + $nunspill_cv unspill events  " *
          "t_final=$(round(times_cv[end]; digits=2))"

    # shared colourmap: upstream (blue) -> downstream (red)
    nt_max = max(nt_bare, nt_cv)
    cmap   = GLMakie.cgrad(:RdYlBu, nt_max; categorical=false, rev=true)

    fig = GLMakie.Figure(size=(1200, 560))

    # --- left panel: reference (no culvert) ---
    ax_bare = GLMakie.Axis(fig[1, 1];
        title   = "Without culvert — $nt_bare traps",
        xlabel  = "cumulative time",
        ylabel  = "fill fraction  V / capacity",
        limits  = (nothing, (-0.04, 1.12)))

    for i in 1:nt_bare
        GLMakie.lines!(ax_bare, times_bare, fracs_bare[i, :];
                       color=cmap[i], linewidth=1.2)
    end
    fill_ts = [e.t for e in evts_bare if e.kind == :fill]
    isempty(fill_ts) ||
        GLMakie.vlines!(ax_bare, fill_ts; color=(:black, 0.12), linewidth=0.7)

    # --- right panel: with culvert ---
    ax_cv = GLMakie.Axis(fig[1, 2];
        title   = "With culvert (r=$(_CV_RADIUS) m) — $nt_cv traps",
        xlabel  = "cumulative time",
        ylabel  = "fill fraction  V / capacity",
        limits  = (nothing, (-0.04, 1.12)))

    for i in 1:nt_cv
        if i == ti_in
            GLMakie.lines!(ax_cv, times_cv, fracs_cv[i, :];
                           color=:magenta, linewidth=2.5, label="inlet trap (drained)")
        elseif i == ti_out
            GLMakie.lines!(ax_cv, times_cv, fracs_cv[i, :];
                           color=:darkcyan, linewidth=2.5, label="outlet trap (boosted)")
        else
            GLMakie.lines!(ax_cv, times_cv, fracs_cv[i, :];
                           color=(cmap[i], 0.45), linewidth=1.0)
        end
    end

    fill_ts_cv    = [e.t for e in evts_cv if e.kind == :fill]
    unspill_ts_cv = [e.t for e in evts_cv if e.kind == :unspill]
    isempty(fill_ts_cv)    ||
        GLMakie.vlines!(ax_cv, fill_ts_cv;    color=(:black,   0.12), linewidth=0.7)
    isempty(unspill_ts_cv) ||
        GLMakie.vlines!(ax_cv, unspill_ts_cv; color=(:magenta, 0.30), linewidth=0.9,
                        label=":unspill event")

    (ti_in !== nothing || ti_out !== nothing) &&
        GLMakie.axislegend(ax_cv; position=:lt, framevisible=true)

    # --- shared colourbar ---
    GLMakie.Colorbar(fig[1, 3];
                     colormap   = GLMakie.cgrad(:RdYlBu; rev=true),
                     limits     = (1, nt_max),
                     label      = "trap index  (upstream → downstream)",
                     flipaxis   = false)

    return fig
end

# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    fig = verify_culvert_network()
    GLMakie.display(fig)
    GLMakie.wait(GLMakie.Screen())
end
