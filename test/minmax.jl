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

@testset "maxima/minima" begin
    @test_throws DomainError Peaks.findnextextrema(<, x1, 1, 0, false)
    @test length(argmaxima(x1)) == 30
    @test length(argminima(x1)) == 30

    @test length(argmaxima(x1, 1; strict=false)) == 31
    @test length(argminima(x1, 1; strict=false)) == 31

    @test length(argmaxima(x1, 1000)) == 4
    @test length(argminima(x1, 1000)) == 4

    @test length(argmaxima(x1, 1000; strict=false)) == 5
    @test length(argminima(x1, 1000; strict=false)) == 5

    @test argmaxima([0,1,0,1,0,1,0,1,0,1,0,1,0,1,0], 2) == [4,8,12]
    @test argmaxima([0,1,0,1,0,1,0,1,0,1,0,1,0,1,0], 2; strict=false) == [2,6,10,14]

    @test simplemaxima(x1) == argmaxima(x1)
    @test simpleminima(x1) == argminima(x1)

    # OffsetArrays
    x1_off = OffsetArray(x1, -200:length(x1)-201)
    @test length(argmaxima(x1_off)) == 30
    @test length(argminima(x1_off)) == 30
    @test argmaxima(x1_off) == argmaxima(x1) .- 201
    @test argminima(x1_off) == argminima(x1) .- 201
    @test length(argmaxima(x1_off, 1; strict=false)) == 31
    @test length(argminima(x1_off, 1; strict=false)) == 31
    @test argmaxima(x1_off, 1; strict=false) == argmaxima(x1, 1; strict=false) .- 201
    @test argminima(x1_off, 1; strict=false) == argminima(x1, 1; strict=false) .- 201
    @test simplemaxima(x1_off) == argmaxima(x1_off)
    @test simpleminima(x1_off) == argminima(x1_off)

    @testset "Plateaus" begin
        p1 = [0,1,1,1,1,1]
        p2 = [0,0,1,1,1,1]
        p3 = [0,0,0,0,1,1]

        @test argmaxima( p1, 2; strict=false) == [2]
        @test argminima(-p1, 2; strict=false) == [2]
        @test argmaxima( p2, 2; strict=false) == [3]
        @test argminima(-p2, 2; strict=false) == [3]
        @test argmaxima( p3, 2; strict=false) == [5]
        @test argminima(-p3, 2; strict=false) == [5]

        @test isempty(argmaxima( p1))
        @test isempty(argminima(-p1))
        @test isempty(argmaxima( p2))
        @test isempty(argminima(-p2))
        @test isempty(argmaxima( p3))
        @test isempty(argminima(-p3))

        @test isempty(simplemaxima( p1))
        @test isempty(simpleminima(-p1))
        @test isempty(simplemaxima( p2))
        @test isempty(simpleminima(-p2))
        @test isempty(simplemaxima( p3))
        @test isempty(simpleminima(-p3))

        @test argmaxima(reverse( p1), 2; strict=false) == [1]
        @test argminima(reverse(-p1), 2; strict=false) == [1]
        @test argmaxima(reverse( p2), 2; strict=false) == [1]
        @test argminima(reverse(-p2), 2; strict=false) == [1]
        @test argmaxima(reverse( p3), 2; strict=false) == [1]
        @test argminima(reverse(-p3), 2; strict=false) == [1]

        @test isempty(argmaxima(reverse( p1)))
        @test isempty(argminima(reverse(-p1)))
        @test isempty(argmaxima(reverse( p2)))
        @test isempty(argminima(reverse(-p2)))
        @test isempty(argmaxima(reverse( p3)))
        @test isempty(argminima(reverse(-p3)))

        @test isempty(simplemaxima(reverse( p1)))
        @test isempty(simpleminima(reverse(-p1)))
        @test isempty(simplemaxima(reverse( p2)))
        @test isempty(simpleminima(reverse(-p2)))
        @test isempty(simplemaxima(reverse( p3)))
        @test isempty(simpleminima(reverse(-p3)))

        @test argmaxima( [0,1,1,1,1,0]) == [2]
        @test argminima(-[0,1,1,1,1,0]) == [2]
        @test argmaxima( [0,1,1,1,2,0]) == [5]
        @test argminima(-[0,1,1,1,2,0]) == [5]
        @test argmaxima( [0,1,1,0,2,1], 3; strict=false) == [5]
        @test argminima(-[0,1,1,0,2,1], 3; strict=false) == [5]

        @test simplemaxima( [0,1,1,1,1,0]) == argmaxima( [0,1,1,1,1,0])
        @test simpleminima(-[0,1,1,1,1,0]) == argminima(-[0,1,1,1,1,0])
        @test simplemaxima( [0,1,1,1,2,0]) == argmaxima( [0,1,1,1,2,0])
        @test simpleminima(-[0,1,1,1,2,0]) == argminima(-[0,1,1,1,2,0])

        # issue #4
        @test isempty(argmaxima(zeros(10)))
        @test argmaxima(zeros(10), 1; strict=false) == [1]

        # issue #30
        y30 = [12.452637, 12.389122, 12.452637, 12.512817, 48.756142, 48.410103, 45.00222]
        @test findnextmaxima(y30, 3, 2) === 5

        # issue #44
        y44 = [-0.0761, -0.0799, -0.0845, -0.0897, -0.0948, -0.0989, -0.1009, -0.1004, -0.0969, -0.0899, -0.0814, -0.0726, -0.065, -0.0589, -0.0554, -0.055, -0.0562, -0.0584, -0.0608, -0.063, -0.065, -0.066, -0.065, -0.0618, -0.0569, -0.0509, -0.0443, -0.0384, -0.0333, -0.0299, -0.0282, -0.0299, -0.0339, -0.0393, -0.0454, -0.05, -0.0526, -0.053, -0.0514, -0.0463, -0.0388, -0.0301, -0.0209, -0.0133, -0.0086, -0.0065, -0.005, -0.004, -0.0043, -0.006, -0.0084, -0.009, -0.0078, -0.0041, -0.0002, 0.0034, 0.0078, 0.0132, 0.0194, 0.0261, 0.032, 0.0358, 0.0379, 0.038, 0.038, 0.0365, 0.0341, 0.0317, 0.0301, 0.0312, 0.0355, 0.0429, 0.0526, 0.0638, 0.075, 0.0865, 0.0975, 0.1066, 0.1127, 0.1173, 0.1216, 0.1256, 0.1294, 0.1329, 0.1367, 0.1415, 0.1471, 0.153, 0.1591, 0.1659, 0.1728, 0.1794, 0.1858, 0.1928, 0.2008, 0.2098, 0.2184, 0.2251, 0.2295, 0.2318, 0.2319, 0.23, 0.2272, 0.2238, 0.2203, 0.2176, 0.2164, 0.217, 0.2158, 0.2117, 0.2039, 0.1924, 0.1772, 0.159, 0.139, 0.1192, 0.1016, 0.086, 0.0713, 0.0586, 0.0475, 0.0369, 0.0257, 0.0147, 0.0045, -0.0047, -0.0133, -0.0218, -0.0299, -0.0374, -0.0453, -0.0532, -0.0601, -0.0662, -0.0716, -0.0763, -0.0798, -0.0825, -0.084, -0.0836, -0.0823, -0.0836, -0.088, -0.0946, -0.1011, -0.1069, -0.112, -0.1164, -0.1196, -0.121, -0.12, -0.1175, -0.1131, -0.1078, -0.1025, -0.0979, -0.0938, -0.0902, -0.0861, -0.0816, -0.0778, -0.0756, -0.0747, -0.0749, -0.0767, -0.0795, -0.0831, -0.087, -0.0909, -0.0936, -0.0946, -0.0942, -0.094, -0.0938, -0.0928, -0.0914, -0.089, -0.0864, -0.0837, -0.082, -0.0838, -0.0868, -0.0891, -0.0907, -0.091, -0.0894, -0.0863, -0.0814, -0.0753, -0.0685, -0.0628, -0.0594, -0.0592, -0.062, -0.067, -0.0724, -0.0773, -0.0808, -0.082, -0.081, -0.0785]
        @test findnextmaxima(y44, 64, 50) == 101
        @test findnextmaxima(y44, 64, 50; strict=false) == 101
        @test findnextmaxima([0,0,1,1,0], 3, 2) == findnextmaxima([0,0,1,0,0], 3, 2)

    end

    # A missing or NaN should not occur within the `w` of the peak
    @testset "Missings and NaNs" begin
        m1 = [0,0,1,1,1,missing,missing,0,0]
        n1 = [0,0,1,1,1,NaN,NaN,0,0]

        @test isempty(argmaxima( m1, 1; strict=true))
        @test isempty(argmaxima( n1, 1; strict=true))
        @test isempty(argminima(-m1, 1; strict=true))
        @test isempty(argminima(-n1, 1; strict=true))
        @test argmaxima( m1, 1; strict=false) == [3,8]
        @test argmaxima( n1, 1; strict=false) == [3,8]
        @test argminima(-m1, 1; strict=false) == [3,8]
        @test argminima(-n1, 1; strict=false) == [3,8]

        @test_throws MethodError simplemaxima(m1)
        @test_throws MethodError simpleminima(m1)
        @test isempty(simplemaxima( n1))
        @test isempty(simpleminima(-n1))

        @test isempty(argmaxima(reverse( m1), 1; strict=true))
        @test isempty(argmaxima(reverse( n1), 1; strict=true))
        @test isempty(argminima(reverse(-m1), 1; strict=true))
        @test isempty(argminima(reverse(-n1), 1; strict=true))
        @test argmaxima(reverse( m1), 1; strict=false) == [1,5]
        @test argmaxima(reverse( n1), 1; strict=false) == [1,5]
        @test argminima(reverse(-m1), 1; strict=false) == [1,5]
        @test argminima(reverse(-n1), 1; strict=false) == [1,5]

        @test isempty(simplemaxima(reverse( n1)))
        @test isempty(simpleminima(reverse(-n1)))

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

        @test simplemaxima(n2) == argmaxima(n2, 1)
        @test simpleminima(-n2) == argminima(-n2, 1)
        @test simplemaxima(reverse(n2)) == argmaxima(reverse(n2), 1)
        @test simpleminima(reverse(-n2)) == argminima(reverse(-n2), 1)

        m3 = [0,1,0,1,missing,1,0,1,0]
        n3 = [0,1,0,1,NaN,1,0,1,0]
        @test argmaxima(m3) == [2,8]
        @test argminima(-m3) == [2,8]
        @test argmaxima(n3) == [2,8]
        @test argminima(-n3) == [2,8]

        @test simpleminima(n3) == argminima(n3)
        @test simpleminima(-n3) == argminima(-n3)

        mn = [1,NaN,missing,1]
        @test isempty(argmaxima(mn))
        @test isempty(argminima(mn))
        @test isempty(argmaxima(reverse(mn)))
        @test isempty(argminima(reverse(mn)))
    end

    @testset "is(minima|maxima|plateau)" begin
        isx =  [0,0,3,1,2,0,4,4,0,5]
        ispk = [0,0,1,0,1,0,1,0,0,1] # Not strict

        pks = argmaxima(isx)
        @test ismaxima(10, isx; strict=true) == false
        @test ismaxima(10, isx; strict=false) == true
        @test ismaxima.(eachindex(isx), Ref(isx); strict=false) == ispk

        pks = argminima(-isx)
        @test isminima(10, -isx; strict=true) == false
        @test isminima(10, -isx; strict=false) == true
        @test isminima.(eachindex(isx), Ref(-isx); strict=false) == ispk

        @test isplateau(7, isx) == true
        @test isplateau(8, isx) == false
        @test isplateau(3, isx) == false

        @test isplateau(10, isx; strict=true) == false
        @test isplateau(10, isx; strict=false) == false
    end

    pks, vals = @test_nowarn findmaxima(x1)
    @test x1[pks] == vals
    @test x1[pks] == maxima(x1)
    pks, vals = @test_nowarn findminima(x1)
    @test x1[pks] == vals
    @test x1[pks] == minima(x1)
end
