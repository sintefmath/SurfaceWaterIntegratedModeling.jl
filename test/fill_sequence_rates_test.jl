# fill_sequence: refreshing a trap whose rates change without its inflow changing.
#
using Test, SurfaceWaterIntegratedModeling

# A stream absorbed inside a non-full trap's footprint, before reaching its bottom, changes only
# that trap's Smax.  Its volume and fill-time estimate must still be refreshed at that event.
@testset "fill_sequence: Smax change without inflow change" begin
    # channel along row 3: trap Y (bottom at j=3) spills down into trap X (bottom at j=9)
    prof = [20.0, 10.0, 5.0, 5.5, 6.0, 4.0, 3.0, 2.0, 1.0, 3.5, 0.5, 0.0]
    z = [10.0 + 0.1j + 0.01i for i in 1:5, j in 1:length(prof)]  # sloped rim, no flats
    z[3, :] = prof
    z[1, :] .+= 10; z[end, :] .+= 10; z[:, 1] .= 30; z[:, end] .= 0.0
    ts = spillanalysis(z)
    @test numtraps(ts) == 2
    Y, X = 1, 2
    @test ts.trap_bottoms[X] == CartesianIndex(3, 9)

    # one strongly infiltrating cell in X's footprint, on the stream's path, absorbing all of
    # Y's spill: X stagnates before Y is full and fills only thanks to Y's spill
    infil = zeros(size(z))
    infil[3, 8] = 42.0
    fill_time(seq) = seq[findfirst(k -> filled_at(seq, k)[X], eachindex(seq))].timestamp
    seq = fill_sequence(ts, [WeatherEvent(0.0, 1.0)]; infiltration = infil)
    tY = seq[findfirst(k -> filled_at(seq, k)[Y], eachindex(seq))].timestamp

    # reference: a new weather period when Y fills, which rebuilds all rates from scratch
    ref = fill_sequence(ts, [WeatherEvent(0.0, 1.0), WeatherEvent(tY, 1.0)];
                        infiltration = infil)
    @test any(k -> filled_at(seq, k)[X], eachindex(seq))
    @test fill_time(seq) ≈ fill_time(ref)
end

# True iff no trap holds water of its own (amount > 0) while one of its subtraps is not full, at
# any event of `seq`.
function _subtraps_full_under_water(ts, seq)
    for k in eachindex(seq)
        filled, amounts = filled_at(seq, k), SurfaceWaterIntegratedModeling.amount_at(seq, k)
        for t in 1:numtraps(ts)
            amounts[t].amount > 0.0 &&
                !all(filled[c] for c in SurfaceWaterIntegratedModeling.subtrapsof(ts, t)) &&
                return false
        end
    end
    return true
end

# A full subtrap under a draining parent is exposed when the parent has drained; that time is
# derived from the parent's rates but stored on the subtrap.  When the parent's rates change, its
# subtraps' estimates must change with them.
@testset "fill_sequence: a parent's rate change reaches its submerged subtraps" begin
    # channel along row 3: U (bottom j=2:3) spills into A (j=7); A and B (j=9) are the subtraps of P
    prof = [30.0, 0.0, 0.0, 9.5, 7.0, 5.0, 1.0, 3.0, 1.0, 5.0, 6.0, 0.5, 0.0]
    z = [10.0 + 0.1j + 0.01i for i in 1:5, j in 1:length(prof)]  # sloped rim, no flats
    z[3, :] = prof
    z[1, :] .+= 10; z[end, :] .+= 10; z[:, 1] .= 30; z[:, end] .= 0.0
    ts = spillanalysis(z)
    @test numtraps(ts) == 4
    U, A, B, P = 1, 2, 3, 4
    @test SurfaceWaterIntegratedModeling.subtrapsof(ts, P) == [A, B]

    infil = zeros(size(z))
    infil[3, 7] = 3.0; infil[3, 9] = 2.0; infil[3, 6] = 1.0; infil[3, 10] = 1.0
    # heavy rain fills A and B and part of P; then light rain lets P drain, until U fills and
    # its spill (through A) makes P fill instead
    rain2 = fill(0.1, size(z)); rain2[:, 1:4] .= 0.4
    w = [WeatherEvent(0.0, 1.0), WeatherEvent(0.5, rain2)]
    seq = fill_sequence(ts, w; infiltration = infil)
    @test _subtraps_full_under_water(ts, seq)

    # reference: a new weather period when U fills, which rebuilds all estimates from scratch
    tU = seq[findfirst(k -> filled_at(seq, k)[U], eachindex(seq))].timestamp
    ref = fill_sequence(ts, [w; WeatherEvent(tU, rain2)]; infiltration = infil)
    fill_time(s, t) = s[findfirst(k -> filled_at(s, k)[t], eachindex(s))].timestamp
    @test fill_time(seq, P) ≈ fill_time(ref, P)
    @test all(filled_at(seq, length(seq))) && all(filled_at(ref, length(ref)))
end

# A basin enclosed by clipped cells has no outer boundary: spill elevation and volume are Inf.
# Once the water stands above its highest cell, every footprint cell is wet and infiltrates.
@testset "fill_sequence: closed basin above its highest cell" begin
    z = [10.0 + 0.1j + 0.01i for i in 1:9, j in 1:12]
    z[5, 5] = 2.0; z[5, 6] = 3.0; z[5, 7] = 1.5        # two pits and a ridge
    mask = falses(size(z)); mask[3:7, 3:9] .= true; mask[4:6, 4:8] .= false  # ring of buildings
    ts = spillanalysis(z; clip_mask = mask)
    @test numtraps(ts) == 3
    P = 3                                               # the merged, enclosed basin
    @test isinf(ts.spillpoints[P].elevation) && isinf(ts.trapvolumes[P])
    @test !any(isnan, SurfaceWaterIntegratedModeling._compute_z_vol_tables(ts)[P][2])

    # rain 1 on the 15 cells; the 10 cells off the middle row infiltrate 1.2, so they hand on
    # nothing and keep 0.2 of spare capacity each: the basin gains 5 - 10 * 0.2 = 3 when all wet
    infil = zeros(size(z)); infil[4, 4:8] .= 1.2; infil[6, 4:8] .= 1.2
    seq = fill_sequence(ts, [WeatherEvent(0.0, 1.0)]; infiltration = infil)
    st = trap_states_at_timepoints(ts, seq, [40.0, 60.0]; verbose = false)
    @test (st[2][2][P] - st[1][2][P]) / 20 ≈ 3.0
end
