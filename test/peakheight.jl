@testset "Peak height" begin
    x = [0,5,2,3,3,1,4,0]
    pks, heights = findmaxima(x)

    @test_throws DimensionMismatch peakheights(pks, heights[1:end-1])

    # TODO: Remove after next breaking release (v0.5)
    @test_logs (:warn, r"renamed") peakheights(pks, heights; maxheight=1)
    @test_logs (:warn, r"renamed") peakheights(pks, heights; minheight=1)
    @test_logs (:warn, r"renamed") peakheights!(copy(pks), copy(heights); maxheight=1)
    @test_logs (:warn, r"renamed") peakheights!(copy(pks), copy(heights); minheight=1)

    filtpks, Y = peakheights(pks, heights)
    @test pks == filtpks
    @test Y == heights

    _, Y = peakheights(pks, heights; max=4)
    @test all(≤(4), Y)

    _, Y = peakheights(pks, heights; min=4)
    @test all(≥(4), Y)

    _, Y = peakheights(pks, heights; max=4.5, min=3.5)
    @test all(x -> 3.5 ≤ x ≤ 4.5, Y)
end
