module MakieExt

using Peaks: ismaxima, isminima, peakproms, peakwidths, interp, drop_irrelevant_side
import Peaks: peaksplot, peaksplot!
using Makie: Makie, ScatterLines, @recipe, automatic, documented_attributes,
             lines!, linesegments!, scatter!

"""
    peaksplot([x,] y, indices)
    peaksplot([x,] pks::NamedTuple)
    peaksplot!(ax, [x,] y, indices)
    peaksplot!(ax, [x,] pks::NamedTuple)

Plot the peaks of a line using Makie

See also: [`findmaxima`](@ref Peaks.findmaxima), [`peakproms`](@ref), [`peakwidths`](@ref),
"""
@recipe PeaksPlot begin
    "Show the reference level for each peak's prominence"
    show_prominences = true

    "Show the left/right bounds of each peak's width"
    show_widths = true

    "The prominence of each peak"
    proms = nothing

    """The edges of each peak. Can either be a zipped vector of left and right edges (e.g.\
    `[(4.5, 6.25)]`) or a tuple of the left and right edges vectors (e.g. `([4.5],\
    [6.25])`)"""
    edges = nothing

    documented_attributes(ScatterLines)...

    "Sets the peak marker"
    marker = @inherit marker
    "Sets the color of peak markers. These default to `color`"
    markercolor = :red
    "Sets the colormap for peak markers. This defaults to `colormap`"
    markercolormap = automatic
    "Sets the colorrange for peak markers. This defaults to `colorrange`"
    markercolorrange = automatic
    "Sets the size of the markers (peak and width if relevant)."
    markersize = @inherit markersize
end

function Makie.convert_arguments(::Type{<:PeaksPlot},
    x::AbstractVector, y::AbstractVector, indices::AbstractVector{Int}
)
    return (x, (;data=y, indices))
end
function Makie.convert_arguments(::Type{<:PeaksPlot},
    y::AbstractVector, indices::AbstractVector{Int}
)
    return (eachindex(y), (;data=y, indices))
end
Makie.convert_arguments(::Type{<:PeaksPlot}, pks::NamedTuple) = (eachindex(pks.data), pks)

Makie.argument_names(::Type{<:PeaksPlot}, N::Integer) = (:x, :pks)

function Makie.plot!(p::PeaksPlot{<:Tuple{<:AbstractVector,<:NamedTuple}})
    show_prominences, show_widths = p.show_prominences[], p.show_widths[]
    show_prominences isa Bool || throw(ArgumentError("Keyword argument `show_prominences` must be a Bool"))
    show_widths isa Bool || throw(ArgumentError("Keyword argument `show_widths` must be a Bool"))

    pks = p.pks[]
    x, y, peaks = p.x[], pks.data, pks.indices

    if ismaxima(first(peaks), y; strict=false)
        maxima = true
    elseif isminima(first(peaks), y; strict=false)
        maxima = false
    else
        throw(ArgumentError("The first peak in `peaks` is not a local extrema"))
    end
    sgn = maxima ? -1 : +1

    xvals = x[peaks]
    yvals = hasproperty(pks, :heights) ? pks.heights : y[peaks]

    # Plot raw data
    lines!(p, x, y;
        color = p.color,
        depth_shift = p.depth_shift,
        linestyle = p.linestyle,
        linewidth = p.linewidth,
        linecap = p.linecap,
        joinstyle = p.joinstyle,
        miter_limit = p.miter_limit,
        colormap = p.colormap,
        colorscale = p.colorscale,
        colorrange = p.colorrange,
        inspectable = p.inspectable,)

    # markercolor is the same as linecolor if left automatic
    map!(p, [:colormap, :markercolormap], :real_markercolormap) do colormap, markercolormap
        return markercolormap === automatic ? colormap : markercolormap
    end
    map!(p, [:colorrange, :markercolorrange], :real_markercolorrange) do colorrange, markercolorrange
        return markercolorrange === automatic ? colorrange : markercolorrange
    end

    if show_prominences || show_widths
        proms = if hasproperty(pks, :proms)
            pks.proms
        elseif !isnothing(p.proms[])
            p.proms[]
        else
            last(peakproms(peaks, y; strict=false))
        end
    end

    if show_prominences
        promlinesx = repeat(xvals; inner=2)
        promlinesy = vec([yvals+sgn*proms yvals]')
        linesegments!(p, promlinesx, promlinesy;
            color=:blue,
            depth_shift = p.depth_shift,
            linewidth = p.linewidth,)

        _, _, plower, pupper = peakwidths(peaks, y, proms;
            strict=false, relheight=prevfloat(1.0))

        plower_x = interp.(Ref(x), drop_irrelevant_side.(plower, peaks, Ref(y), maxima))
        pupper_x = interp.(Ref(x), drop_irrelevant_side.(pupper, peaks, Ref(y), maxima))

        promwidthlinesx = vec([plower_x pupper_x]')
        promwidthlinesy = repeat(yvals + sgn * proms; inner=2)
        linesegments!(p, promwidthlinesx, promwidthlinesy;
            color=:blue,
            linestyle = :dash,
            linewidth = p.linewidth,
            depth_shift = p.depth_shift,)
    end

    if show_widths
        wlower, wupper = if hasproperty(pks, :edges)
            first.(pks.edges), last.(pks.edges)
        elseif p.edges[] isa Vector{Tuple{Float64,Float64}}
            first.(p.edges[]), last.(p.edges[])
        elseif p.edges[] isa Tuple{Vector{Float64},Vector{Float64}}
            p.edges[]
        else # show_widths must be true if we're here
            _, _, ledge, redge = peakwidths(peaks, y, proms; strict=false)
            ledge, redge
        end

        wlower_x = interp(x, wlower)
        wupper_x = interp(x, wupper)

        halfwidthlinesx = vec([wlower_x wupper_x]')
        halfwidthlinesy = repeat(yvals + sgn * proms .* 0.5; inner=2)
        linesegments!(p, halfwidthlinesx, halfwidthlinesy;
            color=:gray,
            linestyle = :dash,
            linewidth = p.linewidth,
            depth_shift = p.depth_shift,)
        scatter!(p, halfwidthlinesx, halfwidthlinesy;
            color=:gray,
            depth_shift = p.depth_shift,)
    end

    scatter!(p, xvals, yvals;
        color = p.markercolor,
        depth_shift = p.depth_shift,
        strokecolor = p.strokecolor,
        strokewidth = p.strokewidth,
        marker = p.marker,
        markersize = p.markersize,
        colormap = p.real_markercolormap,
        colorscale = p.colorscale,
        colorrange = p.real_markercolorrange,
        inspectable = p.inspectable,)

    return p
end

end # module MakieExt
