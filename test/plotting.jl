@testset "Plotting" begin
    # generate data
    t = 0:1/100:1;
    y = 2 * sin.(5 * t) + 3 * sin.(10 * t) + 2 * sin.(30 * t);

    # find and plot maxima
    pks, vals = findmaxima(y);
    pks, proms = peakproms!(pks, y; min=1)
    _, _, edges... = peakwidths!(pks, y, proms)
    pks_nt  = findmaxima(y) |> peakproms!(; min=1) |> peakwidths!

    @testset "Plots.jl recipe" begin
        # first peak isn't an extrema
        @test_throws ArgumentError plotpeaks(t, y, [2])
        @test_throws ArgumentError plotpeaks(t, y, [2]; show_prominences=2)
        @test_throws ArgumentError plotpeaks(t, y, [2]; show_widths=2)
        @test_warn r"deprecated" plotpeaks(t, y; peaks=pks)
        @test_warn r"renamed" plotpeaks(t, y, pks; prominences=false)
        @test_warn r"renamed" plotpeaks(t, y, pks; widths=false)

        # show_prominences=false, show_widths=false
        plt = plotpeaks(t, y, pks; show_prominences=false, show_widths=false)
        @test_reference "references/plots/onlypeaks.png" plt
        plt = plotpeaks(t, y; peaks=pks, show_prominences=false, show_widths=false)
        @test_reference "references/plots/onlypeaks.png" plt
        plt = plotpeaks(t, pks_nt; show_prominences=false, show_widths=false)
        @test_reference "references/plots/onlypeaks.png" plt
        plt = plotpeaks(y, pks; show_prominences=false, show_widths=false)
        @test_reference "references/plots/onlypeaks_nox.png" plt
        plt = plotpeaks(pks_nt; show_prominences=false, show_widths=false)
        @test_reference "references/plots/onlypeaks_nox.png" plt

        # show_prominences=true, show_widths=false
        plt = plotpeaks(t, y, pks; show_widths=false)
        @test_reference "references/plots/peaks_and_proms.png" plt
        plt = plotpeaks(t, y, pks; proms, show_widths=false)
        @test_reference "references/plots/peaks_and_proms.png" plt
        plt = plotpeaks(t, pks_nt; show_widths=false)
        @test_reference "references/plots/peaks_and_proms.png" plt

        # show_prominences=false, show_widths=true
        plt = plotpeaks(t, y, pks; show_prominences=false, show_widths=true)
        @test_reference "references/plots/peaks_and_widths.png" plt
        plt = plotpeaks(t, pks_nt; show_prominences=false, show_widths=true)
        @test_reference "references/plots/peaks_and_widths.png" plt

        # show_prominences=true, show_widths=true
        plt = plotpeaks(t, y, pks)
        @test_reference "references/plots/everything.png" plt
        plt = plotpeaks(t, y, pks; edges)
        @test_reference "references/plots/everything.png" plt
        plt = plotpeaks(t, y, pks; edges=collect(zip(edges...)))
        @test_reference "references/plots/everything.png" plt
        plt = plotpeaks(t, pks_nt)
        @test_reference "references/plots/everything.png" plt

        # add minima to plot
        minpks = findminima(y) |> peakproms!(; min=1) |> peakwidths!
        @test_reference "references/plots/everything_minima.png" plotpeaks(t, minpks)
    end

    @testset "Makie.jl recipe" begin
        # first peak isn't an extrema
        @test_throws ArgumentError peaksplot(t, y, [2])
        @test_throws ArgumentError peaksplot(t, pks; show_prominences=2)
        @test_throws ArgumentError peaksplot(t, pks; show_widths=2)

        # show_prominences=false, show_widths=false
        plt = peaksplot(t, y, pks; show_prominences=false, show_widths=false)
        @test_reference "references/makie/onlypeaks.png" plt
        plt = peaksplot(t, pks_nt; show_prominences=false, show_widths=false)
        @test_reference "references/makie/onlypeaks.png" plt
        plt = peaksplot(y, pks; show_prominences=false, show_widths=false)
        @test_reference "references/makie/onlypeaks_nox.png" plt
        plt = peaksplot(pks_nt; show_prominences=false, show_widths=false)
        @test_reference "references/makie/onlypeaks_nox.png" plt

        # show_prominences=true, show_widths=false
        plt = peaksplot(t, y, pks; show_widths=false)
        @test_reference "references/makie/peaks_and_proms.png" plt
        plt = peaksplot(t, y, pks; proms, show_widths=false)
        @test_reference "references/makie/peaks_and_proms.png" plt
        plt = peaksplot(t, pks_nt; show_widths=false)
        @test_reference "references/makie/peaks_and_proms.png" plt

        # show_prominences=false, show_widths=true
        plt = peaksplot(t, y, pks; show_prominences=false, show_widths=true)
        @test_reference "references/makie/peaks_and_widths.png" plt
        plt = peaksplot(t, pks_nt; show_prominences=false, show_widths=true)
        @test_reference "references/makie/peaks_and_widths.png" plt

        # show_prominences=true, show_widths=true
        plt = peaksplot(t, y, pks)
        @test_reference "references/makie/everything.png" plt
        plt = peaksplot(t, y, pks; edges)
        @test_reference "references/makie/everything.png" plt
        plt = peaksplot(t, y, pks; edges=collect(zip(edges...)))
        @test_reference "references/makie/everything.png" plt
        plt = peaksplot(t, pks_nt)
        @test_reference "references/makie/everything.png" plt

        # add minima to plot
        pks = findminima(y) |> peakproms!(; min=1) |> peakwidths!
        @test_reference "references/makie/everything_minima.png" peaksplot(t, pks)
    end
end
