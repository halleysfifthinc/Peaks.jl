using Peaks
using Test, OffsetArrays, Plots

@testset verbose=true "Peaks" begin
    include("simd.jl")
    include("minmax.jl")
    include("peakprom.jl")
    include("peakwidth.jl")
    include("peakheight.jl")
    include("plotting.jl")
    include("utils.jl")
end
