module Peaks

export argmaxima, argminima, maxima, minima, findmaxima, findminima, findnextmaxima,
    findnextminima, peakproms, peakproms!, peakwidths, peakwidths!, peakheights,
    peakheights!, ismaxima, isminima, filterpeaks!

include("minmax.jl")
include("utils.jl")
include("peakprom.jl")
include("peakwidth.jl")
include("peakheight.jl")
include("plot.jl")

end # module Peaks
