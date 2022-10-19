using RecipesBase

"Get `x[idx::Float]` using linear interpolation."
function interp(x::AbstractVector{<:Real}, idx::Real)
    isinteger(idx) && return x[Int(idx)]
    prev = Int(floor(idx))
    next = Int(ceil(idx))
    return x[prev] + ((idx - prev) * (x[next] - x[prev])) / (next - prev)
end

function interp(x::AbstractVector{<:Real}, idx::AbstractVector{<:Real})
    return interp.(Ref(x), idx)
end

@shorthands peaksplot
@recipe function f(::Type{Val{:peaksplot}}, x, y, z; peaks::AbstractVector{Int}, prominences=false, widths=false)
    if all(ismaxima.(peaks, Ref(y); strict=false))
        maxima = true
    elseif all(isminima.(peaks, Ref(y); strict=false))
        maxima = false
    else
        throw(error("The peaks in `peaks` must *either* local maxima *or* minima."))
    end
    sgn = maxima ? -1 : +1
    ext_color = maxima ? :Red : :Green
    ext_label = maxima ? "maxima" : "minima"

    xvals = x[peaks]
    yvals = y[peaks]

    @series begin  # plot raw data
        linewidth --> 2
        seriestype := :path
        label := "signal"
        linecolor --> :Black
        x := x
        y := y
    end

    if prominences || widths
        _, proms = peakproms(peaks, y; strict=false)
    end

    if prominences
        nans = fill(NaN, (length(peaks), 1))
        promlinesx = vec(vcat(reshape(repeat(xvals; inner=2), 2, length(peaks)), nans'))
        promlinesy = vec([yvals + sgn * proms yvals nans]')

        _, _, lower, upper = peakwidths(peaks, y, proms;
            strict=false, relheight=prevfloat(1.0))

        lower_x = interp(x, lower)
        upper_x = interp(x, upper)

        @series begin  # plot vertical prominence lines
            seriestype := :path
            label := "prominence"
            linewidth := 1
            linecolor := :Blue
            x := promlinesx
            y := promlinesy
        end

        promwidthlinesx = vec([lower_x upper_x nans]')
        promwidthlinesy = vec(vcat(
            reshape(repeat(yvals + sgn * proms; inner=2), 2, length(peaks)),
            nans'))

        @series begin  # plot horizontal prominence lines
            seriestype := :path
            linewidth := 1
            label := nothing
            linecolor := :Blue
            linealpha := 0.5
            linestyle := :dash
            x := promwidthlinesx
            y := promwidthlinesy
        end
    end

    if widths
        nans = fill(NaN, (length(peaks), 1))
        _, _, lower, upper = peakwidths(peaks, y, proms; strict=false)
        lower_x = interp(x, lower)
        upper_x = interp(x, upper)

        halfwidthlinesx = vec([lower_x upper_x nans]')
        halfwidthlinesy = vec(vcat(reshape(repeat(yvals + sgn * proms .* 0.5; inner=2), 2, length(peaks)), nans'))

        @series begin  # plot width lines
            seriestype := :path
            linewidth := 1
            label := "width"
            linewidth := 1
            linecolor := :Gray
            linestyle := :dash
            x := halfwidthlinesx
            y := halfwidthlinesy
        end
        @series begin  # plot width points
            seriestype := :scatter
            markercolor := :Gray
            label := nothing
            x := halfwidthlinesx
            y := halfwidthlinesy
        end
    end

    @series begin  # plot extrema points
        seriestype := :scatter
        label := ext_label
        markercolor := ext_color
        x := xvals
        y := yvals
    end
end
