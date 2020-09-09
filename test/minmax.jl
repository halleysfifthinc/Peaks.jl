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

@testset "argminima/argmaxima" begin
    @test length(argmaxima(x1)) == 30
    @test length(argminima(x1)) == 30

    @test length(argmaxima(x1, 1, false)) == 31
    @test length(argminima(x1, 1, false)) == 31

    # OffsetArrays
    x1_off = OffsetArray(x1, -200:length(x1)-201)
    @test length(argmaxima(x1_off)) == 30
    @test length(argminima(x1_off)) == 30
    @test argmaxima(x1_off) == argmaxima(x1) .- 201
    @test argminima(x1_off) == argminima(x1) .- 201
    @test length(argmaxima(x1_off, 1, false)) == 31
    @test length(argminima(x1_off, 1, false)) == 31
    @test argmaxima(x1_off, 1, false) == argmaxima(x1, 1, false) .- 201
    @test argminima(x1_off, 1, false) == argminima(x1, 1, false) .- 201

    @test length(argmaxima(x1, 1000)) == 4
    @test length(argminima(x1, 1000)) == 4

    @test length(argmaxima(x1, 1000, false)) == 5
    @test length(argminima(x1, 1000, false)) == 5

    @test argmaxima([0,1,0,1,0,1,0,1,0,1,0,1,0,1,0], 2) == [4,8,12]
    @test argmaxima([0,1,0,1,0,1,0,1,0,1,0,1,0,1,0], 2, false) == [2,6,10,14]

    @testset "Plateaus" begin
        p1 = [0,1,1,1,1,1]
        p2 = [0,0,1,1,1,1]
        p3 = [0,0,0,0,1,1]

        @test argmaxima( p1, 2, false) == [2]
        @test argminima(-p1, 2, false) == [2]
        @test argmaxima( p2, 2, false) == [3]
        @test argminima(-p2, 2, false) == [3]
        @test argmaxima( p3, 2, false) == [5]
        @test argminima(-p3, 2, false) == [5]

        @test isempty(argmaxima( p1))
        @test isempty(argminima(-p1))
        @test isempty(argmaxima( p2))
        @test isempty(argminima(-p2))
        @test isempty(argmaxima( p3))
        @test isempty(argminima(-p3))

        @test argmaxima(reverse( p1), 2, false) == [1]
        @test argminima(reverse(-p1), 2, false) == [1]
        @test argmaxima(reverse( p2), 2, false) == [1]
        @test argminima(reverse(-p2), 2, false) == [1]
        @test argmaxima(reverse( p3), 2, false) == [1]
        @test argminima(reverse(-p3), 2, false) == [1]

        @test isempty(argmaxima(reverse( p1)))
        @test isempty(argminima(reverse(-p1)))
        @test isempty(argmaxima(reverse( p2)))
        @test isempty(argminima(reverse(-p2)))
        @test isempty(argmaxima(reverse( p3)))
        @test isempty(argminima(reverse(-p3)))

        @test argmaxima( [0,1,1,1,1,0]) == [2]
        @test argminima(-[0,1,1,1,1,0]) == [2]
        @test argmaxima( [0,1,1,1,2,0]) == [5]
        @test argminima(-[0,1,1,1,2,0]) == [5]
        @test argmaxima( [0,1,1,0,2,1], 3, false) == [5]
        @test argminima(-[0,1,1,0,2,1], 3, false) == [5]

        # issue #4
        @test isempty(argmaxima(zeros(10)))
        @test argmaxima(zeros(10), 1, false) == [1]
    end

    # A missing or NaN should not occur within the `w` of the peak
    @testset "Missings and NaNs" begin
        m1 = [0,0,1,1,1,missing,missing,0,0]
        n1 = [0,0,1,1,1,NaN,NaN,0,0]
        @test isempty(argmaxima( m1, 1, true))
        @test isempty(argmaxima( n1, 1, true))
        @test isempty(argminima(-m1, 1, true))
        @test isempty(argminima(-n1, 1, true))
        @test argmaxima( m1, 1, false) == [3,8]
        @test argmaxima( n1, 1, false) == [3,8]
        @test argminima(-m1, 1, false) == [3,8]
        @test argminima(-n1, 1, false) == [3,8]

        @test isempty(argmaxima(reverse( m1), 1, true))
        @test isempty(argmaxima(reverse( n1), 1, true))
        @test isempty(argminima(reverse(-m1), 1, true))
        @test isempty(argminima(reverse(-n1), 1, true))
        @test argmaxima(reverse( m1), 1, false) == [1,5]
        @test argmaxima(reverse( n1), 1, false) == [1,5]
        @test argminima(reverse(-m1), 1, false) == [1,5]
        @test argminima(reverse(-n1), 1, false) == [1,5]

        m2 = [0,1,2,1,missing,missing,missing]
        n2 = [0,1,2,1,NaN,NaN,NaN]
        @test argmaxima(m2, 1) == [3]
        @test argmaxima(n2, 1) == [3]
        @test isempty(argmaxima(m2, 2))
        @test isempty(argmaxima(n2, 2))
        @test argminima(-m2, 1) == [3]
        @test argminima(-n2, 1) == [3]
        @test isempty(argminima(-m2, 2))
        @test isempty(argminima(-n2, 2))
        @test argmaxima(reverse(m2), 1) == [5]
        @test argmaxima(reverse(n2), 1) == [5]
        @test isempty(argmaxima(reverse(m2), 2))
        @test isempty(argmaxima(reverse(n2), 2))
        @test argminima(reverse(-m2), 1) == [5]
        @test argminima(reverse(-n2), 1) == [5]
        @test isempty(argminima(reverse(-m2), 2))
        @test isempty(argminima(reverse(-n2), 2))

        m3 = [0,1,0,1,missing,1,0,1,0]
        n3 = [0,1,0,1,NaN,1,0,1,0]
        @test argmaxima(m3) == [2,8]
        @test argminima(-m3) == [2,8]
        @test argmaxima(n3) == [2,8]
        @test argminima(-n3) == [2,8]

        mn = [1,NaN,missing,1]
        @test isempty(argmaxima(mn))
        @test isempty(argminima(mn))
        @test isempty(argmaxima(reverse(mn)))
        @test isempty(argminima(reverse(mn)))
    end
end
