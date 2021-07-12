module Peaks

using Compat

export argmaxima, argminima, findmaxima, findminima, findnextmaxima, findnextminima,
    peakproms, peakproms!, peakwidths, peakwidths!

export Maxima, Minima

include("minmax.jl")
include("peakprom.jl")
include("peakwidth.jl")

end # module Peaks
