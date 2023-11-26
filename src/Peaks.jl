module Peaks

using Compat

export argmaxima, argminima, maxima, minima, findmaxima, findminima, findnextmaxima, findnextminima, peakproms, peakproms!, peakwidths, peakwidths!, peakheights, peakheights!, ismaxima, isminima, findpeaks, filterpeaks!

include("minmax.jl")
include("utils.jl")
include("peakprom.jl")
include("peakwidth.jl")
include("peakheight.jl")
include("api_rework.jl")
include("plot.jl")

end # module Peaks
