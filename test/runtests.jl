using Peaks
using Test, OffsetArrays, Plots, Documenter

DocMeta.setdocmeta!(Peaks, :DocTestSetup, :(using Peaks); recursive=true)

@testset verbose=true "Peaks" begin
    include("simd.jl")
    include("minmax.jl")
    include("peakprom.jl")
    include("peakwidth.jl")
    include("peakheight.jl")
    include("plotting.jl")
    include("utils.jl")

    @testset "Doctests" begin
        doctest(Peaks; manual=false)
    end
end
