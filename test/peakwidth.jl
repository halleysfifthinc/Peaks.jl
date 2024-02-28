@testset "Peak width" begin
    fs = 100
    T = 1/fs
    t = 0:T:6pi+T
    sint = sin.(t)

    @testset "minima/maxima" begin
        sinpks = argmaxima(sint)
        _, widths, _, _ = peakwidths(sinpks, sint, sint[sinpks]; strict=false, relheight=1)
        @test widths ≈ fill(pi*100, length(sinpks)) atol=.01

        sinpks = argminima(sint)
        _, widths, _, _ = peakwidths(sinpks, sint, abs.(sint[sinpks]); strict=false, relheight=1)
        @test widths ≈ fill(pi*100, length(sinpks)) atol=.01
    end

    _, widths, _, _ = peakwidths([2], [0.,1.,0.], [1.])
    @test widths == [1.]
    _, widths, _, _ = peakwidths([2], [0.,1.,0.], [NaN])
    @test widths[1] === NaN
    _, widths, _, _ = peakwidths([2], [0.,1.,0.], [missing])
    @test widths[1] === missing
    _, widths, _, _ = peakwidths([2], [0.,1.,NaN], [1.]; strict=true)
    @test widths[1] === NaN
    _, widths, _, _ = peakwidths([2], [0.,1.,0.,-1.], [1.]; strict=false)
    _, widthsnan, _, _ = peakwidths([2], [0.,1.,NaN,-1.], [1.]; strict=false)
    @test widths == widthsnan

    @testset "strict" begin
        sinpks = argmaxima(sint)
        _, widths, _, _ = peakwidths(sinpks, sint, ones(length(sinpks)); strict=true, relheight=1)
        @test first(widths) === NaN
        _, widths, _, _ = peakwidths(sinpks, sint, ones(length(sinpks)); strict=false, relheight=1)
        @test first(widths) !== NaN
    end

    @testset "Min/max width" begin
        sinpks = argmaxima(sint)
        _, proms = peakproms(sinpks, sint)

        @test length(first(peakwidths(sinpks, sint, proms; min=pi*75))) == 2
        @test length(first(peakwidths(sinpks, sint, proms; max=pi*75))) == 1

        @test_throws ArgumentError peakwidths(1:3, ones(3), ones(3); max=0.1, min=1)

        # TODO: Remove after next breaking release (v0.5)
        @test_logs (:warn, r"renamed") peakwidths(sinpks, sint, proms; maxwidth=1)
        @test_logs (:warn, r"renamed") peakwidths(sinpks, sint, proms; minwidth=1)
        @test_logs (:warn, r"renamed") peakwidths!(copy(sinpks), copy(sint), copy(proms); maxwidth=1)
        @test_logs (:warn, r"renamed") peakwidths!(copy(sinpks), copy(sint), copy(proms); minwidth=1)
    end

end
