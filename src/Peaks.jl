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

function bla(x)
    pks = findpeaks(x)
    pks = peakproms!(pks, min=0.5)
    pks = peakwidths!(pks)
    pks = peakheights!(pks)
    return pks
end

function bla2(x)
    pks = findpeaks(x) |> peakproms!(min=0.5) |> peakwidths!() |> peakheights!()
    return pks
end

export bla, bla2
end # module Peaks
