let
    @testset "Plotting" begin
        # generate data
        t = 0:1/100:1
        y = 2 * sin.(5 * t) + 3 * sin.(10 * t) + 2 * sin.(30 * t)

        # find and plot maxima
        pks, vals = findmaxima(y)
        pks, proms = peakproms!(pks, y; minprom=1)

        plt = plotpeaks(t, y, peaks=pks, prominences=true, widths=true)
        @test plt isa Plots.Plot

        savepath_png = abspath(joinpath(@__DIR__, "..", "docs", "src", "assets", "images", "maxima_prom_width.png"))
        savefig(plt, savepath_png)

        # add minima to plot
        pks, vals = findminima(y)
        pks, proms = peakproms!(pks, y; minprom=1)
        plt = peaksplot!(t, y, peaks=pks, prominences=true, widths=true)
        @test plt isa Plots.Plot
        savepath_png = abspath(joinpath(@__DIR__, "..", "docs", "src", "assets", "images", "extrema_prom_width.png"))
        savefig(plt, savepath_png)

        plt = plotpeaks(t, y, peaks=pks, prominences=true, widths=true)
        savepath_png = abspath(joinpath(@__DIR__, "..", "docs", "src", "assets", "images", "minima_prom_width.png"))
        savefig(plt, savepath_png)
    end
end
