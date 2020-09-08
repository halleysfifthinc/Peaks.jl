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

    @test length(maxima(x1, 1, false)) == 31
    @test length(minima(x1, 1, false)) == 31

    # OffsetArrays
    x1_off = OffsetArray(x1, -200:length(x1)-201)
    @test length(maxima(x1_off)) == 30
    @test length(minima(x1_off)) == 30
    @test maxima(x1_off) == maxima(x1) .- 201
    @test minima(x1_off) == minima(x1) .- 201
    @test length(maxima(x1_off, 1, false)) == 31
    @test length(minima(x1_off, 1, false)) == 31
    @test maxima(x1_off, 1, false) == maxima(x1, 1, false) .- 201
    @test minima(x1_off, 1, false) == minima(x1, 1, false) .- 201

    @test length(maxima(x1, 1000)) == 4
    @test length(minima(x1, 1000)) == 4

    @test length(maxima(x1, 1000, false)) == 5
    @test length(minima(x1, 1000, false)) == 5

    @test maxima([0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],2) == [4,8,12]
    @test maxima([0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],2,false) == [2,6,10,14]

    @testset "Plateaus" begin
        p1 = [0,1,1,1,1,1]
        p2 = [0,0,1,1,1,1]
        p3 = [0,0,0,0,1,1]

        @test maxima( p1, 2, false) == [2]
        @test minima(-p1, 2, false) == [2]
        @test maxima( p2, 2, false) == [3]
        @test minima(-p2, 2, false) == [3]
        @test maxima( p3, 2, false) == [5]
        @test minima(-p3, 2, false) == [5]

        @test isempty(maxima( p1))
        @test isempty(minima(-p1))
        @test isempty(maxima( p2))
        @test isempty(minima(-p2))
        @test isempty(maxima( p3))
        @test isempty(minima(-p3))

        @test maxima(reverse( p1), 2, false) == [1]
        @test minima(reverse(-p1), 2, false) == [1]
        @test maxima(reverse( p2), 2, false) == [1]
        @test minima(reverse(-p2), 2, false) == [1]
        @test maxima(reverse( p3), 2, false) == [1]
        @test minima(reverse(-p3), 2, false) == [1]

        @test isempty(maxima(reverse( p1)))
        @test isempty(minima(reverse(-p1)))
        @test isempty(maxima(reverse( p2)))
        @test isempty(minima(reverse(-p2)))
        @test isempty(maxima(reverse( p3)))
        @test isempty(minima(reverse(-p3)))

        @test maxima( [0,1,1,1,1,0]) == [2]
        @test minima(-[0,1,1,1,1,0]) == [2]
        @test maxima( [0,1,1,1,2,0]) == [5]
        @test minima(-[0,1,1,1,2,0]) == [5]
        @test maxima( [0,1,1,0,2,1], 3, false) == [5]
        @test minima(-[0,1,1,0,2,1], 3, false) == [5]

        # issue #4
        @test isempty(maxima(zeros(10)))
        @test maxima(zeros(10),1,false) == [1]
    end

    # A missing or NaN should not occur within the `w` of the peak
    @testset "Missings and NaNs" begin
        m1 = [0,0,1,1,1,missing,missing,0,0]
        n1 = [0,0,1,1,1,NaN,NaN,0,0]
        @test isempty(maxima( m1, 1, true))
        @test isempty(maxima( n1, 1, true))
        @test isempty(minima(-m1, 1, true))
        @test isempty(minima(-n1, 1, true))
        @test maxima( m1, 1, false) == [3,8]
        @test maxima( n1, 1, false) == [3,8]
        @test minima(-m1, 1, false) == [3,8]
        @test minima(-n1, 1, false) == [3,8]

        @test isempty(maxima(reverse( m1), 1, true))
        @test isempty(maxima(reverse( n1), 1, true))
        @test isempty(minima(reverse(-m1), 1, true))
        @test isempty(minima(reverse(-n1), 1, true))
        @test maxima(reverse( m1), 1, false) == [1,5]
        @test maxima(reverse( n1), 1, false) == [1,5]
        @test minima(reverse(-m1), 1, false) == [1,5]
        @test minima(reverse(-n1), 1, false) == [1,5]

        m2 = [0,1,2,1,missing,missing,missing]
        n2 = [0,1,2,1,NaN,NaN,NaN]
        @test maxima(m2, 1) == [3]
        @test maxima(n2, 1) == [3]
        @test isempty(maxima(m2, 2))
        @test isempty(maxima(n2, 2))
        @test minima(-m2, 1) == [3]
        @test minima(-n2, 1) == [3]
        @test isempty(minima(-m2, 2))
        @test isempty(minima(-n2, 2))
        @test maxima(reverse(m2), 1) == [5]
        @test maxima(reverse(n2), 1) == [5]
        @test isempty(maxima(reverse(m2), 2))
        @test isempty(maxima(reverse(n2), 2))
        @test minima(reverse(-m2), 1) == [5]
        @test minima(reverse(-n2), 1) == [5]
        @test isempty(minima(reverse(-m2), 2))
        @test isempty(minima(reverse(-n2), 2))

        m3 = [0,1,0,1,missing,1,0,1,0]
        n3 = [0,1,0,1,NaN,1,0,1,0]
        @test maxima(m3) == [2,8]
        @test minima(-m3) == [2,8]
        @test maxima(n3) == [2,8]
        @test minima(-n3) == [2,8]

        mn = [1,NaN,missing,1]
        @test isempty(maxima(mn))
        @test isempty(minima(mn))
        @test isempty(maxima(reverse(mn)))
        @test isempty(minima(reverse(mn)))
    end
end
