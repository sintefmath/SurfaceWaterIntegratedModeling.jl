# Second pass: sample both runs over time; check global mass balance, flow-field
# bounds, and gauge series (with vs without NBS).
include("setup.jl")
using Printf
JLD2.@load "run1.jld2" r1 r2 tstruct dtm
const G = SWIM.Graphs
gs = size(dtm)
const RAINR = 120.0 / (1000 * 60)       # m/min per cell
const NCELL = prod(gs)

# terminal cells that drain out of the domain
exitcells = [c for c in 1:NCELL if isempty(G.outneighbors(tstruct.flowgraph, c)) && tstruct.regions[c] <= 0]
println("exit cells: ", length(exitcells), "; regions ", unique(tstruct.regions) |> sort)
for t in 1:numtraps(tstruct)
    println("trap $t vol=", round(tstruct.trapvolumes[t], digits=1), " sp=", tstruct.spillpoints[t],
            " subtraps=", subtrapsof(tstruct, t), " fp=", length(tstruct.footprints[t]))
end

function sample(r, tps; withnbs=true)
    seq = withnbs ? r.seq : r.seq0
    nbs = withnbs ? r.nbs : DynNBSPlacement[]
    st = all_states_at_timepoints(tstruct, seq, tps; nbs = nbs)
    out = []
    for (k, s) in enumerate(st)
        f = network_watercourses(tstruct, s.filled, Dict{Int,Float64}(), s.nbs_layers;
                                 nbs = nbs, precipitation = RAINR)
        trapvol = sum(s.amounts)
        nbsvol  = isempty(s.nbs_layers) ? 0.0 : sum(sum(v) for v in values(s.nbs_layers))
        push!(out, (; t = tps[k], f, trapvol, nbsvol, s))
    end
    return out
end

flowls(x) = x * 1000 / 60   # m3/min -> l/s

function report(name, r)
    println("\n================ ", name)
    tps = collect(0.0:10.0:1000.0)
    A = sample(r, tps; withnbs=true)
    B = sample(r, tps; withnbs=false)
    nbscells = isempty(r.nbs) ? Int[] : reduce(vcat, [n.footprint for n in r.nbs])
    # bound: static all-full no-NBS field (max conceivable flow per cell)
    fmax, = watercourses(tstruct, fill(true, numtraps(tstruct)); precipitation = RAINR)
    println(" t    | S_trap     S_nbs     | dS/dt   out     rain   resid(%) | negflow  overmax(cells,max excess l/s)")
    for k in 2:length(A)-1
        a = A[k]
        dS = ((A[k+1].trapvol + A[k+1].nbsvol) - (A[k-1].trapvol + A[k-1].nbsvol)) / (A[k+1].t - A[k-1].t)
        out = sum(max(a.f[c], 0.0) for c in exitcells)
        rain = RAINR * NCELL
        resid = (dS + out - rain) / rain * 100
        fin = [c for c in 1:NCELL if !(c in nbscells) && isfinite(a.f[c])]
        neg = count(c -> a.f[c] < -1e-9 && false, fin)  # negatives are legit (infil capacity); none here
        over = [c for c in fin if a.f[c] > fmax[c] * (1 + 1e-6) + 1e-9]
        mx = isempty(over) ? 0.0 : maximum(flowls(a.f[c] - fmax[c]) for c in over)
        if mod(a.t, 50) == 0 || abs(resid) > 1 || !isempty(over)
            @printf("%5.0f | %9.2f %9.2f | %7.4f %7.4f %7.4f %8.2f | %d  %d %.3f\n",
                    a.t, a.trapvol, a.nbsvol, dS, out, rain, resid, neg, length(over), mx)
        end
    end
    println("--- gauges (l/s for flow cells; m3 for trap cells), NBS vs no-NBS")
    for g in r.gp
        c = LinearIndices(gs)[g]
        reg = tstruct.regions[c]
        intrap = any(c in tstruct.footprints[t] for t in 1:numtraps(tstruct))
        innbs = c in nbscells
        print(rpad(string(g), 22), intrap ? "trap " : innbs ? "nbs  " : "flow ")
        for k in 1:10:length(A)
            v  = A[k].f[c]; v0 = B[k].f[c]
            @printf(" %6.0f:%7.3f/%7.3f", A[k].t, flowls(v), flowls(v0))
        end
        println()
    end
    return A, B
end

A1, B1 = report("setup 1", r1)
A2, B2 = report("setup 2", r2)
JLD2.@save "run2.jld2" A1 B1 A2 B2
