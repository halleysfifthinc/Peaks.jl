let
    @testset "Plotting" begin
        # generate data
        t = 0:1/100:1
        y = 2 * sin.(5 * t) + 3 * sin.(10 * t) + 2 * sin.(30 * t)

        # find and plot maxima
        pks, vals = findmaxima(y)
        pks, proms = peakproms!(pks, y; minprom=1)

        plt = peaksplot(t, y, peaks=pks)
        @test plt isa Plots.Plot

        savepath_png = abspath(joinpath(@__DIR__, "..", "docs", "src", "assets", "images", "peaks_prom_width.png"))
        savefig(plt, savepath_png)
        @info "Plots saved to <$savepath_png>"

        # add minima to plot
        pks, vals = findminima(y)
        pks, proms = peakproms!(pks, y; minprom=1)

        plt = peaksplot!(t, y, peaks=pks, maxima=false)
        @test plt isa Plots.Plot

        savepath_png = abspath(joinpath(@__DIR__, "..", "docs", "src", "assets", "images", "extrema_prom_width.png"))
        savefig(plt, savepath_png)
        @info "Plots saved to <$savepath_png>"
    end
end
