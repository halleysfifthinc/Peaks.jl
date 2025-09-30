module Peaks

using SIMD: SIMD, Vec, VecTypes, vload

export argmaxima, simplemaxima, argminima, simpleminima, maxima, minima, findmaxima,
    findminima, findnextmaxima, findnextminima, peakproms, peakproms!, peakwidths,
    peakwidths!, peakheights, peakheights!, findpeaks, ismaxima, isminima, isplateau,
    filterpeaks!, plotpeaks, plotpeaks!, peaksplot, peaksplot!

include("simple.jl")
include("minmax.jl")
include("utils.jl")
include("peakprom.jl")
include("peakwidth.jl")
include("peakheight.jl")

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if (exc.f === simplemaxima || exc.f === simpleminima) && eltype(argtypes[1]) >: Missing
            if exc.f === simplemaxima
                cmp = "max"
            else
                cmp = "min"
            end
            printstyled(io, "\nsimple$(cmp)ima"; color=:cyan)
            print(io, " does not support vectors with ")
            printstyled(io, "missing"; color=:cyan)
            print(io, "s. Use the ")
            printstyled(io, "arg$(cmp)ima"; color=:cyan)
            print(io, " function to find $(cmp)ima when ")
            printstyled(io, "missing"; color=:cyan)
            print(io, "s are present.")
        end
    end
end

"""
    plotpeaks(x, y, indices; show_prominences=true, show_widths=true, kwargs...)
    plotpeaks(x, pks::NamedTuple; show_prominences=true, show_widths=true, kwargs...)
    plotpeaks!(plt, x, y, indices; show_prominences=true, show_widths=true, kwargs...)
    plotpeaks!(plt, x, pks::NamedTuple; show_prominences=true, show_widths=true, kwargs...)

Plot the peaks of a line; optionally, also show the reference level for each peak's
prominence, and the left/right bounds of each peak's width. This is a Plots.jl recipe.

See also: [`findmaxima`](@ref), [`peakproms`](@ref), [`peakwidths`](@ref),

#### Keyword arguments:
- `proms`: The prominence of each peak. The `proms` keyword will be ignored if a `pks::NamedTuple` is given and has a defined `proms` property.
- `edges`: The edges of each peak. Can either be a zipped vector of left and right
  edges (e.g. `[(4.5, 6.25)]`) or a tuple of the left and right edges vectors (e.g. `([4.5],
  [6.25])`). The `edges` keyword will be ignored if a `pks::NamedTuple` is given and has a defined `edges` property.
"""
function plotpeaks end
function plotpeaks! end
@doc (@doc plotpeaks) plotpeaks!

function peaksplot end
function peaksplot! end

end # module Peaks
