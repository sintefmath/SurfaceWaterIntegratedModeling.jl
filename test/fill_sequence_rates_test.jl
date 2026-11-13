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
