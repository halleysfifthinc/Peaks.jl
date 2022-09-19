let
    @testset "Plotting" begin
        t = 0:1/100:1
        y = 2*sin.(5*t)+3*sin.(10*t)+2*sin.(30*t)
        pks, vals = findmaxima(y)
        pks, proms = peakproms!(pks, y; minprom=1)

        plt = peaksplot(t, y, peaks=pks)

        @test plt isa Plots.Plot
        
        savepath_png = abspath(joinpath(@__DIR__, "..", "docs","src","assets","images","peaks_prom_width.png"))
        savefig(plt, savepath_png)
        @info "Plots saved to <$savepath_png>"
    end
end