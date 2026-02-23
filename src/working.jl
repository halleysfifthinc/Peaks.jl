using Peaks, Test
using Peaks: findall_offset, _simd_extrema!

# ENV["JULIA_DEBUG"] = "Peaks"

function _simdmaxima(x)
    pks = BitVector(undef, length(x))
    fill!(pks, false)

    _simd_extrema!(pks, <, x)
    return pks
end

for pklen in 3:96
    arr = repeat([0,1]; outer=cld(pklen,2))
    @test @views argmaxima(arr[1:pklen]) == simplemaxima(arr[1:pklen])
    # @test argmaxima(@view(arr[1:pklen])) == findall_offset(_simdmaxima(@view(arr[1:pklen])), 0)
    @test @views argmaxima(arr[1:pklen]) == findall_offset(_simdmaxima(arr[1:pklen]), 0)

    arr = repeat([1,0]; outer=cld(pklen,2))
    @test @views argmaxima(arr[1:pklen]) == simplemaxima(arr[1:pklen])
    # @test argmaxima(@view(arr[1:pklen])) == findall_offset(_simdmaxima(@view(arr[1:pklen])), 0)
    @test @views argmaxima(arr[1:pklen]) == findall_offset(_simdmaxima(arr[1:pklen]), 0)

end
println("All tests passed!")
