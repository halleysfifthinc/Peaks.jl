@testset "Peak height" begin
    x = [0,5,2,3,3,1,4,0]
    pks = argmaxima(x)

    filtpks, Y = peakheights(pks, x)
    @test pks == filtpks
    @test Y == x[pks]

    _, Y = peakheights(pks, x; maxheight=4)
    @test all(≤(4), Y)

    _, Y = peakheights(pks, x; minheight=4)
    @test all(≥(4), Y)

    _, Y = peakheights(pks, x; maxheight=4.5, minheight=3.5)
    @test all(x -> 3.5 ≤ x ≤ 4.5, Y)
end
