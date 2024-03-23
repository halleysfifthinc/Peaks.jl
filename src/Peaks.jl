module Peaks

using Compat

export argmaxima, simplemaxima, argminima, simpleminima, maxima, minima, findmaxima,
    findminima, findnextmaxima, findnextminima, peakproms, peakproms!, peakwidths,
    peakwidths!, peakheights, peakheights!, ismaxima, isminima, isplateau, filterpeaks!

include("minmax.jl")
include("utils.jl")
include("peakprom.jl")
include("peakwidth.jl")
include("peakheight.jl")
include("plot.jl")

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if exc.f == simplemaxima || exc.f == simpleminima
            if exc.f == simplemaxima
                cmp = "max"
            else
                cmp = "min"
            end
            printstyled(io, "\nsimple$(cmp)ima"; color=:cyan)
            print(io, " does not support vectors with ")
            printstyled(io, "missing"; color=:cyan)
            print(io, "s. See the ")
            printstyled(io, "arg$(cmp)ima"; color=:cyan)
            print(io, " function to find $(cmp)ima when ")
            printstyled(io, "missing"; color=:cyan)
            print(io, "s are present.")
        end
    end
end

end # module Peaks
