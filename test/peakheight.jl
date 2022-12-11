@testset "Peak height" begin
    x = [0,5,2,3,3,1,4,0]
    pks, heights = findmaxima(x)


    filtpks, Y = peakheights(pks, heights)
    @test pks == filtpks
    @test Y == heights

    _, Y = peakheights(pks, heights; maxheight=4)
    @test all(≤(4), Y)

    _, Y = peakheights(pks, heights; minheight=4)
    @test all(≥(4), Y)

    _, Y = peakheights(pks, heights; maxheight=4.5, minheight=3.5)
    @test all(x -> 3.5 ≤ x ≤ 4.5, Y)
end
