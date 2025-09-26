module PlotsExt

using Peaks: Peaks, ismaxima, isminima, peakproms, peakwidths, interp, drop_irrelevant_side
import Peaks: plotpeaks, plotpeaks!
using RecipesBase: RecipesBase, @recipe, @series, @userplot

@userplot struct PlotPeaks{T,X<:AbstractVector,Y<:AbstractVector{T}}
    x::X
    data::Y
    indices::Vector{Int}
    proms::Union{Nothing,Vector{T}}
    edges::Union{Nothing,Vector{Tuple{Float64,Float64}}}

    function PlotPeaks(args)
        if 1 ≤ length(args) ≤ 2
            pks = last(args)
            if pks isa NamedTuple
                if issubset((:indices, :data), propertynames(pks))
                    data = pks.data
                    indices = pks.indices
                    proms = hasproperty(pks, :proms) ? pks.proms : nothing
                    edges = hasproperty(pks, :edges) ? pks.edges : nothing
                    x = if length(args) === 2
                        length(first(args)) === length(data) || throw(ArgumentError("length(x) doesn't match length(data)"))
                        first(args)
                    else
                        eachindex(data)
                    end
                    T = eltype(data)
                    X = typeof(x)
                    Y = typeof(data)
                    return new{T,X,Y}(x, data, indices, proms, edges)
                else
                    throw(ArgumentError("an incomplete `pks` NamedTuple was provided. `pks` must have at least :indices and :data fields"))
                end
            elseif length(args) === 2
                x = first(args)
                data = pks
                T = eltype(data)
                X = typeof(x)
                Y = typeof(data)
                return new{T,X,Y}(x, data, Int[], nothing, nothing)
            else
                throw(ArgumentError("not enough arguments provided"))
            end
        elseif length(args) === 3
            T = eltype(args[2])
            X = typeof(args[1])
            Y = typeof(args[2])
            return new{T,X,Y}(args..., nothing, nothing)
        end
    end
end
@recipe function f(_peaks::PlotPeaks{T};
    show_prominences=true, show_widths=true, proms=nothing, edges=nothing,
    prominences=nothing, widths=nothing, peaks=nothing
) where T
    show_prominences isa Bool ||
        throw(ArgumentError("Keyword argument `show_prominences` must be a Bool"))
    show_widths isa Bool ||
        throw(ArgumentError("Keyword argument `show_widths` must be a Bool"))
    if !isnothing(edges)
        edges isa Union{Vector{Tuple{Float64,Float64}},Tuple{Vector{Float64},Vector{Float64}}} ||
            throw(ArgumentError("the type of `edges` is not acceptable. Accepted types include `Vector{Tuple{Float64,Float64}}` and `Tuple{Vector{Float64},Vector{Float64}}`"))
    end
    if !isnothing(peaks)
        Base.depwarn("Providing `peaks` as a kwarg is deprecated. Please use `plotpeaks(x, y, peaks; kwargs...)`.", :plotpeaks)
        !isempty(_peaks.indices) && throw(ArgumentError("combining new `plotpeaks` API with deprecated kwargs is not permitted. Got `plotpeaks(x, y, peaks; peaks, kwargs...)`."))
    end
    if !isnothing(prominences)
        Base.depwarn("The `prominences` kwarg has been renamed to `show_prominences`.", :plotpeaks)
        show_prominences = prominences
    end
    if !isnothing(widths)
        Base.depwarn("The `widths` kwarg has been renamed to `show_widths`.", :plotpeaks)
        show_widths = widths
    end

    if isnothing(peaks) && isempty(_peaks.indices)
        y = _peaks.x
        peaks = _peaks.data
        x = eachindex(y)
    else
        x, y = _peaks.x, _peaks.data
        if isnothing(peaks)
            peaks = _peaks.indices
        end
    end

    if ismaxima(first(peaks), y; strict=false)
        maxima = true
    elseif isminima(first(peaks), y; strict=false)
        maxima = false
    else
        throw(ArgumentError("The first value in `peaks` is not a local extrema"))
    end
    sgn = maxima ? -1 : +1
    ext_color = maxima ? :Red : :Green
    ext_label = maxima ? "maxima" : "minima"

    xvals = x[peaks]
    yvals = y[peaks]

    @series begin  # plot raw data
        linewidth --> 2
        seriestype := :path
        label --> "signal"
        linecolor --> :Black
        x, y
    end

    @series begin  # plot extrema points
        seriestype := :scatter
        label := ext_label
        markercolor --> ext_color
        xvals, yvals
    end

    if show_prominences || show_widths
        proms = if !isnothing(_peaks.proms)
            _peaks.proms
        elseif !isnothing(proms)
            proms
        else
            last(peakproms(peaks, y; strict=false))
        end
    end

    if show_prominences
        nans = fill(NaN, (length(peaks), 1))
        promlinesx = vec(vcat(reshape(repeat(xvals; inner=2), 2, length(peaks)), nans'))
        promlinesy = vec([yvals + sgn * proms yvals nans]')

        _, _, lower, upper = peakwidths(peaks, y, proms;
            strict=false, relheight=prevfloat(1.0))

        lower_x = interp.(Ref(x), drop_irrelevant_side.(lower, peaks, Ref(y), maxima))
        upper_x = interp.(Ref(x), drop_irrelevant_side.(upper, peaks, Ref(y), maxima))

        promwidthlinesx = vec([lower_x upper_x nans]')
        promwidthlinesy = vec(vcat(
            reshape(repeat(yvals + sgn * proms; inner=2), 2, length(peaks)),
            nans'))

        @series begin  # plot vertical prominence lines
            seriestype := :path
            label := "prominence"
            linewidth := 1
            linecolor := :Blue
            promlinesx, promlinesy
        end
        @series begin  # plot horizontal prominence lines
            seriestype := :path
            linewidth := 1
            label := nothing
            linecolor := :Blue
            linealpha := 0.5
            linestyle := :dash
            promwidthlinesx, promwidthlinesy
        end
    end

    if show_widths
        lower, upper = if !isnothing(_peaks.edges)
            first.(_peaks.edges), last.(_peaks.edges)
        elseif edges isa Vector{Tuple{Float64,Float64}}
            first.(edges), last.(edges)
        elseif edges isa Tuple{Vector{Float64},Vector{Float64}}
            edges
        else # show_widths must be true if we're here
            _, _, ledge, redge = peakwidths(peaks, y, proms; strict=false)
            ledge, redge
        end

        lower_x = interp(x, lower)
        upper_x = interp(x, upper)

        nans = fill(NaN, (length(peaks), 1))
        halfwidthlinesx = vec([lower_x upper_x nans]')
        halfwidthlinesy = vec(vcat(reshape(repeat(yvals + sgn * proms .* 0.5; inner=2), 2, length(peaks)), nans'))

        @series begin  # plot width lines
            seriestype := :path
            linewidth := 1
            label := "width"
            linewidth := 1
            linecolor := :Gray
            linestyle := :dash
            halfwidthlinesx, halfwidthlinesy
        end
        @series begin  # plot width points
            seriestype := :scatter
            markercolor := :Gray
            label := nothing
            halfwidthlinesx, halfwidthlinesy
        end
    end

    nothing
end

end # module PlotsExt
