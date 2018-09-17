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

    @test length(maxima(x1, 1000)) == 4
    @test length(minima(x1, 1000)) == 4

    @test length(maxima(x1, 1000, true)) == 5
    @test length(minima(x1, 1000, true)) == 5
end
