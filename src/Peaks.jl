module Peaks

using Compat

export argmaxima, argminima, findmaxima, findminima, findnextmaxima, findnextminima,
    peakprom

export Maxima, Minima

include("minmax.jl")
include("peakprom.jl")

end # module Peaks
