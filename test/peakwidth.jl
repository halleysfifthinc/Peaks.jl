@testset "Peak width" begin
    fs = 100
    T = 1/fs
    t = 0:T:6pi+T
    sint = sin.(t)

    @testset "minima/maxima" begin
        sinpks = argmaxima(sint)
        _, widths, _, _ = peakwidth(sinpks, sint, sint[sinpks]; strictbounds=false, relheight=1)
        @test widths ≈ fill(pi*100, length(sinpks)) atol=.01

        sinpks = argminima(sint)
        _, widths, _, _ = peakwidth(sinpks, sint, abs.(sint[sinpks]); strictbounds=false, relheight=1)
        @test widths ≈ fill(pi*100, length(sinpks)) atol=.01
    end

    _, widths, _, _ = peakwidth([2], [0.,1.,0.], [1.])
    @test widths == [1.]
    _, widths, _, _ = peakwidth([2], [0.,1.,0.], [NaN])
    @test widths[1] === NaN
    _, widths, _, _ = peakwidth([2], [0.,1.,0.], [missing])
    @test widths[1] === missing
    _, widths, _, _ = peakwidth([2], [0.,1.,NaN], [1.]; strictbounds=true)
    @test widths[1] === NaN
    _, widths, _, _ = peakwidth([2], [0.,1.,0.,-1.], [1.]; strictbounds=false)
    _, widthsnan, _, _ = peakwidth([2], [0.,1.,NaN,-1.], [1.]; strictbounds=false)
    @test widths == widthsnan

    @testset "strictbounds" begin
        sinpks = argmaxima(sint)
        _, widths, _, _ = peakwidth(sinpks, sint, ones(length(sinpks)); strictbounds=true, relheight=1)
        @test first(widths) === NaN
        _, widths, _, _ = peakwidth(sinpks, sint, ones(length(sinpks)); strictbounds=false, relheight=1)
        @test first(widths) !== NaN
    end

    @testset "Min/max width" begin
        sinpks = argmaxima(sint)
        _, proms = peakprom(sinpks, sint)

        @test length(first(peakwidth(sinpks, sint, proms; minwidth=pi*75))) == 2
        @test length(first(peakwidth(sinpks, sint, proms; maxwidth=pi*75))) == 1

        @test_throws ArgumentError peakwidth(1:3, ones(3), ones(3); maxwidth=0.1, minwidth=1)
    end

end
