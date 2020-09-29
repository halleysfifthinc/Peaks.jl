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
        pi, pp = peakprom(x1)
        ni, np = peakprom(Minima(), -x1)

        @test pi == ni
        @test pp == np
    end


    @testset "Prominence values" begin
        i, p = peakprom(sin.(1e-5:1e-5:9*pi))

        @test p[[1,5]] ≈ [1., 1.] atol=1e-4
        @test p[[2,3,4]] ≈ [2., 2., 2.] atol=1e-4

        @test last(peakprom([1,0,2,0,1])) == [2]
        @test last(peakprom([1,0,2,0,1]; strictbounds=false)) == [1,2,1]

        # Prominence should be the same regardless of window size
        p1 = [0,0,3,1,2,0,4,0,0,5]
        @test last(peakprom(p1)) == [3,1,4]
        @test last(peakprom(p1, 2)) == [3,4]

        # A peaks of the same height count as an intersection for reference intervals
        p4 = [0,4,2,4,3,4,0]
        @test last(peakprom(p4)) == [2,1,1]
        @test last(peakprom(p4[2:end-1])) == [1]
        @test last(peakprom(p4[2:end-1]; strictbounds=false)) == [2,1,1]

        # The presence of a missing/NaN in either bounding interval poisons the prominence
        m4 = [missing; p4; missing]
        n4 = [NaN; p4; NaN]
        @test isequal(last(peakprom(m4)), [missing,1,missing])
        @test isequal(last(peakprom(n4)), [NaN,1.,NaN])
        @test last(peakprom(m4; strictbounds=false)) == [2,1,1]
        @test last(peakprom(n4; strictbounds=false)) == [2.,1.,1.]

        @test_skip peakprom([missing,1,missing]; strictbounds=false) == [missing]

        p5 = [-1,6,3,4,2,4,2,5,-2,0]
        @test last(peakprom(p5, 3; strictbounds=false)) == [7,3]
        @test last(peakprom(reverse(p5), 3; strictbounds=false)) == [3,7]


    end

    @testset "Minimum prominence" begin
        minprom = 1.5
        i, p = peakprom(sin.(1e-5:1e-5:9*pi))
        mi, mp = peakprom(sin.(1e-5:1e-5:9*pi); minprom=minprom)
        @test all(x -> x >= minprom, mp)
    end

    # issue #4
    let i, p
        i, p = peakprom(zeros(10))
        @test isempty(i)
        @test isempty(p)
    end
end

