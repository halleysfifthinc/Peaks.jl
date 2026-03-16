@testset "Peak height" begin
    x = [0,5,2,3,3,1,4,0]
    pks, heights = findmaxima(x)

    @test_throws DimensionMismatch peakheights(pks, heights[1:end-1])

    filtpks, Y = peakheights(pks, heights)
    @test pks == filtpks
    @test Y == heights

    _, Y = peakheights(pks, heights; max=4)
    @test all(≤(4), Y)

    _, Y = peakheights(pks, heights; min=4)
    @test all(≥(4), Y)

    _, Y = peakheights(pks, heights; max=4.5, min=3.5)
    @test all(x -> 3.5 ≤ x ≤ 4.5, Y)

    @testset "Filtering empty peaks (issue #84)" begin
        ipks, Y = peakheights(Int[], Float64[])
        @test isempty(ipks) && isempty(Y)

        # test both the raw/original API and the NamedTuple API
        # (they are filtered by different codepaths)
        _pks, _ = peakheights(pks, heights; max=0)
        @test isempty(_pks)

        _pks = peakheights(findmaxima(x); max=0)
        @test isempty(_pks.indices)
    end
end
