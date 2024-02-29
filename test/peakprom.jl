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
        _, maxprom = peakproms(maxpks, x1)
        _, minprom = peakproms(minpks, -x1)
        @test maxprom == minprom

        _, maxprom = peakproms(maxpks, x1; strict=false)
        _, minprom = peakproms(minpks, -x1; strict=false)
        @test maxprom == minprom
    end

    @testset "Prominence values" begin
        x2 = sin.(1e-5:1e-5:9*pi)
        i, p = peakproms(argmaxima(x2), x2)

        @test p[[1,5]] ≈ [1., 1.] atol=1e-4
        @test p[[2,3,4]] ≈ [2., 2., 2.] atol=1e-4

        p2 = [1,0,2,0,1]
        @test last(peakproms(argmaxima(p2), p2)) == [2]
        @test last(peakproms(argmaxima(p2; strict=false), p2; strict=false)) == [1,2,1]

        # Prominence should be the same regardless of window size
        p1 = [0,0,3,1,2,0,4,0,0,5]
        @test last(peakproms(argmaxima(p1, 1), p1)) == [3,1,4]
        @test last(peakproms(argmaxima(p1, 2), p1)) == [3,4]

        # A peaks of the same height count as an intersection for reference intervals
        p4 = [0,4,2,4,3,4,0]
        @test last(peakproms(argmaxima(p4), p4)) == [2,1,1]
        @test last(peakproms(argmaxima(p4[2:end-1]), p4[2:end-1])) == [1]
        @test last(peakproms(argmaxima(p4[2:end-1]; strict=false), p4[2:end-1]; strict=false)) == [2,1,1]

        # The presence of a missing/NaN in either bounding interval poisons the prominence
        m4 = [missing; p4; missing]
        n4 = [NaN; p4; NaN]
        @test isequal(last(peakproms(argmaxima(m4), m4)), [missing,1,missing])
        @test isequal(last(peakproms(argmaxima(n4), n4)), [NaN,1.,NaN])
        @test last(peakproms(argmaxima(m4), m4; strict=false)) == [2,1,1]
        @test last(peakproms(argmaxima(n4), n4; strict=false)) == [2.,1.,1.]

        m5 = [missing, 1, missing]
        n5 = [NaN, 1, NaN]
        @test last(peakproms(argmaxima(m5; strict=false), m5; strict=false)) == [0]
        @test last(peakproms(argmaxima(n5; strict=false), n5; strict=false)) == [0]

        p5 = [-1,6,3,4,2,4,2,5,-2,0]
        @test last(peakproms(argmaxima(p5, 3; strict=false), p5; strict=false)) == [7,3]
        @test last(peakproms(argmaxima(reverse(p5), 3; strict=false), reverse(p5); strict=false)) == [3,7]

        # strict prominence on unstrict peaks give unknown prominence
        p2_unstrict = argmaxima(p2; strict=false)
        @test isequal(last(peakproms(p2_unstrict, p2)), [missing, 2, missing])
    end

    @testset "Min/max prominence" begin
        sint = sin.(T:T:6pi)
        maxs = argmaxima(sint)
        @test length(first(peakproms(maxs, sint; min=1.5))) == 2
        @test length(first(peakproms(maxs, sint; max=1.5))) == 1

        @test_throws ArgumentError peakproms([1,2,3], sint; max=0.1, min=1)

        # TODO: Remove after next breaking release (v0.5)
        @test_logs (:warn, r"renamed") peakproms(maxs, sint; maxprom=1)
        @test_logs (:warn, r"renamed") peakproms(maxs, sint; minprom=1)
        @test_logs (:warn, r"renamed") peakproms!(copy(maxs), copy(sint); maxprom=1)
        @test_logs (:warn, r"renamed") peakproms!(copy(maxs), copy(sint); minprom=1)
    end

    @test_throws ArgumentError peakproms([2], 1:10)

    # issue #4
    let i, p
        i, p = peakproms(Int[], zeros(10))
        @test isempty(i)
        @test isempty(p)
        i, p = peakproms([1], zeros(10); strict=false)
        @test p == [0.0]
    end

    issue91ddaa9 = [1,2,3,2,3,1]
    @test all(!iszero,
        last(peakproms(argminima(issue91ddaa9; strict=false), issue91ddaa9; strict=false)))

end
