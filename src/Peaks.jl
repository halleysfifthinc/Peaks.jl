module Peaks

using Compat

# Old exports:
# export argmaxima, argminima, maxima, minima, findmaxima, findminima, findnextmaxima, findnextminima, peakproms, peakproms!, peakwidths, peakwidths!, peakheights, peakheights!, ismaxima, isminima, findpeaks
# Proposed exports:
export findpeaks, peakproms!, peakwidths!, peakheights!

include("minmax.jl")
include("peakprom.jl")
include("peakwidth.jl")
include("peakheight.jl")
include("api_rework.jl")
include("plot.jl")

end # module Peaks
