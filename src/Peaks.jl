module Peaks

using Compat

export argmaxima, argminima, maxima, minima, findmaxima, findminima, findnextmaxima,
    findnextminima, peakproms, peakproms!, peakwidths, peakwidths!, peakheights,
    peakheights!

include("minmax.jl")
include("peakprom.jl")
include("peakwidth.jl")
include("peakheight.jl")
include("plot.jl")

end # module Peaks
