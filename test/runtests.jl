using Peaks
using Test, OffsetArrays, Plots, Documenter

DocMeta.setdocmeta!(Peaks, :DocTestSetup, :(using Peaks); recursive=true)

include("minmax.jl")
include("peakprom.jl")
include("peakwidth.jl")
include("peakheight.jl")
include("plotting.jl")

@testset "Doctests" begin
    doctest(Peaks; manual=false)
end
