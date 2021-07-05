a = 3
b = 2
c = 1

fs = 100
T = 1/fs
f1 = 5
f2 = 10
f3 = 30

t = T:T:fs

x1 = a*sin.(2*pi*f1*T*t)+b*sin.(2*pi*f2*T*t)+c*sin.(2*pi*f3*T*t);

@testset "Peak prominence" begin
    @testset "Reciprocity" begin
        maxpks = argmaxima(x1)
        minpks = argminima(-x1)
        _, maxprom = peakprom(maxpks, x1)
        _, minprom = peakprom(minpks, -x1)
        @test maxprom == minprom

        _, maxprom = peakprom(maxpks, x1; strictbounds=false)
        _, minprom = peakprom(minpks, -x1; strictbounds=false)
        @test maxprom == minprom
    end

    @testset "Prominence values" begin
        x2 = sin.(1e-5:1e-5:9*pi)
        i, p = peakprom(argmaxima(x2), x2)

        @test p[[1,5]] ≈ [1., 1.] atol=1e-4
        @test p[[2,3,4]] ≈ [2., 2., 2.] atol=1e-4

        p2 = [1,0,2,0,1]
        @test last(peakprom(argmaxima(p2), p2)) == [2]
        @test last(peakprom(argmaxima(p2; strictbounds=false), p2; strictbounds=false)) == [1,2,1]

        # Prominence should be the same regardless of window size
        p1 = [0,0,3,1,2,0,4,0,0,5]
        @test last(peakprom(argmaxima(p1, 1), p1)) == [3,1,4]
        @test last(peakprom(argmaxima(p1, 2), p1)) == [3,4]

        # A peaks of the same height count as an intersection for reference intervals
        p4 = [0,4,2,4,3,4,0]
        @test last(peakprom(argmaxima(p4), p4)) == [2,1,1]
        @test last(peakprom(argmaxima(p4[2:end-1]), p4[2:end-1])) == [1]
        @test last(peakprom(argmaxima(p4[2:end-1]; strictbounds=false), p4[2:end-1]; strictbounds=false)) == [2,1,1]

        # The presence of a missing/NaN in either bounding interval poisons the prominence
        m4 = [missing; p4; missing]
        n4 = [NaN; p4; NaN]
        @test isequal(last(peakprom(argmaxima(m4), m4)), [missing,1,missing])
        @test isequal(last(peakprom(argmaxima(n4), n4)), [NaN,1.,NaN])
        @test last(peakprom(argmaxima(m4; strictbounds=false), m4; strictbounds=false)) == [2,1,1]
        @test last(peakprom(argmaxima(n4; strictbounds=false), n4; strictbounds=false)) == [2.,1.,1.]

        m5 = [missing, 1, missing]
        n5 = [NaN, 1, NaN]
        @test last(peakprom(argmaxima(m5; strictbounds=false), m5; strictbounds=false)) == [0]
        @test last(peakprom(argmaxima(n5; strictbounds=false), n5; strictbounds=false)) == [0]

        p5 = [-1,6,3,4,2,4,2,5,-2,0]
        @test last(peakprom(argmaxima(p5, 3; strictbounds=false), p5; strictbounds=false)) == [7,3]
        @test last(peakprom(argmaxima(reverse(p5), 3; strictbounds=false), reverse(p5); strictbounds=false)) == [3,7]


    end

    @testset "Minimum prominence" begin
        minprom = 1.5
        x2 = sin.(1e-5:1e-5:9*pi)
        max2 = argmaxima(x2)
        _, p = peakprom(max2, x2)
        _, mp = peakprom(max2, x2; minprom=minprom)
        @test all(x -> x >= minprom, mp)
    end

    # issue #4
    let i, p
        i, p = peakprom(Int[], zeros(10))
        @test isempty(i)
        @test isempty(p)
        i, p = peakprom([1], zeros(10); strictbounds=false)
        @test p == [0.0]
    end

    @test_deprecated peakprom(Minima(), x1)
    @test_deprecated peakprom(Maxima(), x1)
end

