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

@testset "Minima/maxima" begin
    @test length(maxima(x1)) == 30
    @test length(minima(x1)) == 30

    @test length(maxima(x1, 1, true)) == 31
    @test length(minima(x1, 1, true)) == 31

    # OffsetArrays
    x1_off = OffsetArray(x1, -200:length(x1)-201)
    @test length(maxima(x1_off)) == 30
    @test length(minima(x1_off)) == 30
    @test maxima(x1_off) == maxima(x1) .- 201
    @test minima(x1_off) == minima(x1) .- 201
    @test length(maxima(x1_off, 1, true)) == 31
    @test length(minima(x1_off, 1, true)) == 31
    @test maxima(x1_off, 1, true) == maxima(x1, 1, true) .- 201
    @test minima(x1_off, 1, true) == minima(x1, 1, true) .- 201

    @test length(maxima(x1, 1000)) == 4
    @test length(minima(x1, 1000)) == 4

    @test length(maxima(x1, 1000, true)) == 5
    @test length(minima(x1, 1000, true)) == 5

    # issue #4
    @test isempty(maxima(zeros(10)))

    @testset "Plateaus" begin
        @test maxima([0,1,1,1,1,0]) == [2]
        @test minima(-[0,1,1,1,1,0]) == [2]

        @test isempty(maxima([0,1,1,1,1,1]))
        @test isempty(minima(-[0,1,1,1,1,1]))

        @test maxima([0,1,1,1,2,0]) == [5]
        @test minima(-[0,1,1,1,2,0]) == [5]
    end

    # A missing or NaN should not occur within the `w` of the peak
    @testset "Missings and NaNs" begin
        m1 = [0,0,1,1,1,missing,missing,0,1,1,1]
        n1 = [0,0,1,1,1,NaN,NaN,0,1,1,1]
        @test isempty(maxima(m1,  1))
        @test isempty(maxima(n1,  1))
        @test isempty(minima(-n1, 1))
        @test isempty(minima(-m1, 1))

        m2 = [0,1,2,1,missing,missing,missing]
        n2 = [0,1,2,1,NaN,NaN,NaN]
        @test maxima(m2, 1) == [3]
        @test isempty(maxima(m2, 2))
        @test maxima(n2, 1) == [3]
        @test isempty(maxima(n2, 2))
        @test maxima(reverse(m2), 1) == [5]
        @test isempty(maxima(reverse(m2), 2))
        @test maxima(reverse(n2), 1) == [5]
        @test isempty(maxima(reverse(n2), 2))

        @test maxima([m2; 1;2;2;0], 1) == [3,9]
        @test maxima([n2; 1;2;2;0], 1) == [3,9]
        @test isempty(maxima([m2; 1;2;2;0], 2))
        @test isempty(maxima([n2; 1;2;2;0], 2))
        @test minima(-[m2; 1;2;2;0], 1) == [3;9]
        @test minima(-[n2; 1;2;2;0], 1) == [3,9]
        @test isempty(minima(-[m2; 1;2;2;0], 2))
        @test isempty(minima(-[n2; 1;2;2;0], 2))
    end
end
