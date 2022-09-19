using RecipesBase

"Get `x[idx::Float]` using linear interpolation."
function interp(x::AbstractVector{<:Real}, idx::Real)
    isinteger(idx) && return x[Int(idx)]
    idx1 = Int(floor(idx))
    idx2 = Int(ceil(idx))
    return x[idx1] + ((idx - idx1) * (x[idx2] - x[idx1])) / (idx2 - idx1)
end

function interp(x::AbstractVector{<:Real}, idx::AbstractVector{<:Real})
    return interp.(Ref(x), idx)
end

@shorthands peaksplot
@recipe function f(::Type{Val{:peaksplot}}, x, y, z; peaks::AbstractArray{<:Integer}, maxima=true)
    sgn = maxima ? -1 : +1
    ext_color = maxima ? :Red : :Green
    ext_label = maxima ? "maxima" : "minima"

    pks, proms = peakproms(peaks, y)
    vals = y[pks]
    promlinesx = vec(vcat(reshape(repeat(x[pks]; inner=2), 2, length(pks)), fill(NaN, (1, length(pks)))))
    promlinesy = vec([vals + sgn * proms vals fill(NaN, (length(pks), 1))]')

    _, widths, lower, upper = peakwidths(pks, y, proms; strict=true)
    lower_x = interp(x, lower)
    upper_x = interp(x, upper)
    halfwidthlinesx = vec([lower_x upper_x fill(NaN, (length(pks), 1))]')
    halfwidthlinesy = vec(vcat(reshape(repeat(vals + sgn * proms .* 0.5; inner=2), 2, length(pks)), fill(NaN, (1, length(pks)))))

    _, widths, lower, upper = peakwidths(pks, y, proms; strict=true, relheight=prevfloat(1.0))
    lower_x = interp(x, lower)
    upper_x = interp(x, upper)
    fullwidthlinesx = vec([lower_x upper_x fill(NaN, (length(pks), 1))]')
    fullwidthlinesy = vec(vcat(reshape(repeat(vals + sgn * proms .* 1; inner=2), 2, length(pks)), fill(NaN, (1, length(pks)))))

    @series begin  # plot raw data
        linewidth --> 2
        seriestype := :path
        label := "signal"
        linecolor --> :Black
        x, y
    end
    @series begin  # plot vertical prominence lines
        seriestype := :path
        label := "prominence"
        linewidth := 1
        linecolor := :Blue
        x := promlinesx
        y := promlinesy
    end
    @series begin  # plot horizontal prominence lines
        seriestype := :path
        linewidth := 1
        label := nothing
        linecolor := :Blue
        x := fullwidthlinesx
        y := fullwidthlinesy
    end
    @series begin  # plot extrema points
        seriestype := :scatter
        label := ext_label
        markercolor := ext_color
        x := x[pks]
        y := vals
    end
    @series begin  # plot width lines
        seriestype := :path
        linewidth := 1
        label := "width"
        linewidth := 1
        linecolor := :Gray
        linestyle := :dot
        x := halfwidthlinesx
        y := halfwidthlinesy
    end
    @series begin  # plot width points
        seriestype := :scatter
        markercolor := :Gray
        markerstrokecolor := nothing
        label := nothing
        x := halfwidthlinesx
        y := halfwidthlinesy
    end
end
