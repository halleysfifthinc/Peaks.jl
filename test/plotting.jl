@testset "Plotting" begin
    # generate data
    t = 0:1/100:1;
    y = 2 * sin.(5 * t) + 3 * sin.(10 * t) + 2 * sin.(30 * t);

    # find and plot maxima
    pks, vals = findmaxima(y);
    pks, proms = peakproms!(pks, y; min=1)
    _, _, edges... = peakwidths!(pks, y, proms)

    plt = plotpeaks(t, y, pks; prominences=false, widths=false)
    @test_reference "references/onlypeaks.png" plt

    plt = plotpeaks(t, y, pks; prominences=true, widths=false)
    @test_reference "references/peaks_and_proms.png" plt
    plt = plotpeaks(t, y, pks; prominences=proms, widths=false)
    @test_reference "references/peaks_and_proms.png" plt

    plt = plotpeaks(t, y, pks; prominences=false, widths=true)
    @test_reference "references/peaks_and_widths.png" plt

    plt = plotpeaks(t, y, pks; prominences=true, widths=true)
    @test_reference "references/everything.png" plt
    plt = plotpeaks(t, y, pks; prominences=true, edges)
    @test_reference "references/everything.png" plt
    plt = plotpeaks(t, y, pks; prominences=true, edges=collect(zip(edges...)))
    @test_reference "references/everything.png" plt


    # add minima to plot
    pks = findminima(y) |> peakproms!(; min=1) |> peakwidths!
    @test_reference "references/everything_minima.png" plotpeaks(t, y, pks.indices)

    # first peak isn't an extrema
    @test_throws ArgumentError plotpeaks(t, y, [2])
end
