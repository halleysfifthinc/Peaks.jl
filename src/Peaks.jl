module Peaks

using Compat

export argmaxima, argminima, findmaxima, findminima, peakprom

export Maxima, Minima

include("minmax.jl")
include("peakprom.jl")

end # module Peaks
