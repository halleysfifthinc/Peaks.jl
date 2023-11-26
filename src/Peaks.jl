module Peaks

using Compat

# Function related to locating peaks
export argmaxima, argminima
export maxima, minima
export findmaxima, findminima
export findnextmaxima, findnextminima
export ismaxima, isminima

# Functions related to working with a set of peaks
export peakproms,   peakproms!
export peakwidths,  peakwidths!
export peakheights, peakheights!
export filterpeaks!

include("minmax.jl")
include("utils.jl")
include("peakprom.jl")
include("peakwidth.jl")
include("peakheight.jl")
include("plot.jl")

end # module Peaks
