using Peaks: lowest_set_bits, highest_set_bits, matching_bit_runs_mask_highest_bit,
    matching_bit_runs_mask_lowest_bit, _simd_extrema!, _simpleextrema_base

_simplemaxima(x) = _simpleextrema_base(<, x, Val(:packed))
_simpleminima(x) = _simpleextrema_base(>, x, Val(:packed))

function _simdmaxima(x)
    pks = BitVector(undef, length(x))
    fill!(pks, false)

    _simd_extrema!(pks, <, x)
    return pks
end

function _simdminima(x)
    pks = BitVector(undef, length(x))
    fill!(pks, false)

    _simd_extrema!(pks, >, x)
    return pks
end

collect_bits(x::Integer) = ntuple(i -> !iszero((x >> (i-1)) & one(typeof(x))), sizeof(x)*8 )

function next_zero(x::Integer, bit)
    bit == sizeof(x)*8 && return bit+1
    while !iszero(x & 2^bit) && bit < sizeof(x)*8
        bit += 1
    end
    return bit
end

function next_one(x::Integer, bit)
    bit == sizeof(x)*8 && return bit+1
    while iszero(x & 2^bit) && bit < sizeof(x)*8
        bit += 1
    end
    return bit
end

function prev_zero(x::Integer, bit)
    bit < 1 && return -1
    while bit >= 0 && !iszero(x & 2^bit)
        bit -= 1
    end
    return bit
end

function prev_one(x::Integer, bit)
    bit < 1 && return -1
    while bit >= 0 && iszero(x & 2^bit)
        bit -= 1
    end
    return bit
end

function manual_test_matching_bit_runs_mask_lowest_bit(bits, mask)
    res = zero(typeof(bits))
    shift = trailing_zeros(mask)+1
    to = 0
    two = oftype(bits, 2)
    nbits = sizeof(bits)*8
    while shift < nbits && to < nbits
        # by definition (starting from trailing_zeros(mask)+1, and continuing from
        # `next_one(mask, to)`) the current bit `mask & 2^shift` is set
        to = next_zero(bits, shift)
        res |= bits & (widen(two)^to - two^(shift-1)) % typeof(bits)
        shift = next_one(mask, to)
    end
    return res
end

function manual_test_matching_bit_runs_mask_highest_bit(bits, mask)
    res = zero(typeof(bits))
    to = sizeof(bits)*8
    shift = to - leading_zeros(mask)-1
    two = oftype(bits, 2)
    while shift >= 0 && to > 0
        # by definition (starting from trailing_zeros(mask)+1, and continuing from
        # `next_one(mask, to)`) the current bit `mask & 2^shift` is set
        to = prev_zero(bits, shift)
        res |= bits & (widen(two)^(shift+1) - two^(to+1)) % typeof(bits)
        shift = prev_one(mask, to)
    end
    return res
end

function test_bit_run_mask_permutations(::Type{T}, mask_type=:lowest) where {T}
    mask_type in (:lowest, :highest) ||
        throw(ArgumentError("`mask_type` must be :lowest or :highest; got $mask_type"))

    if mask_type == :lowest
        mask_func = lowest_set_bits
        reference_matching_func = manual_test_matching_bit_runs_mask_lowest_bit
        test_func = matching_bit_runs_mask_lowest_bit
    else
        mask_func = highest_set_bits
        reference_matching_func = manual_test_matching_bit_runs_mask_highest_bit
        test_func = matching_bit_runs_mask_highest_bit
    end

    for i in 0:typemax(T)
        base_mask = mask_func(i)
        @test test_func(i, base_mask) == reference_matching_func(i, base_mask)
        zeroing = next_one(base_mask, 0)
        for _ in 1:count_ones(base_mask)
            mask = base_mask & ~(2^zeroing) # zero one bit-run

            pass = @test test_func(i, mask) == reference_matching_func(i, mask)
            if !(pass isa Test.Pass)
                @error i, mask
            end

            zeroing = next_one(base_mask, zeroing+1)
        end
    end
end

@testset "maxima/minima: simd edition" begin
    @test lowest_set_bits(0b01110101) == 0b00010101
    @test highest_set_bits(0b01110101) == 0b01000101

    # Obvious manually created bits/mask
    @test matching_bit_runs_mask_lowest_bit(0b01110101, 0b00010001) == 0b01110001
    @test matching_bit_runs_mask_highest_bit(0b01110101, 0b01000001) == 0b01110001

    test_bit_run_mask_permutations(UInt16, :lowest)
    test_bit_run_mask_permutations(UInt16, :highest)

    # Single element peak that is moved through packed peaks chunks and (at some point)
    # crosses a chunk boundary
    arr = zeros(64*2)
    for pki in 2:(2*64-3)
        arr[pki] = 1
        @test simplemaxima(arr) == [pki]
        pass = @test _simdmaxima(arr) == _simplemaxima(arr)
        if !(pass isa Test.Pass)
            @error pki
        end
        arr[pki] = 0
    end

    arr = zeros(64*3)
    for pki in 2:(3*64-3)
        for plat_len in 1:(64+16) # max plateau length longer than 64-bit chunk
            arr[pki:min(pki+plat_len,lastindex(arr)-1)] .= 1
            @test argmaxima(arr) == [pki]
            @test simplemaxima(arr) == [pki]
            pass = @test _simdmaxima(arr) == _simplemaxima(arr)
            if !(pass isa Test.Pass)
                @error pki, plat_len
            end
            arr[min(pki+plat_len,lastindex(arr)-1)] = 2
            @test simplemaxima(arr) == [min(pki+plat_len,lastindex(arr)-1)]
            @test argmaxima(arr) == [min(pki+plat_len,lastindex(arr)-1)]
            pass = @test _simdmaxima(arr) == _simplemaxima(arr)
            if !(pass isa Test.Pass)
                @error pki, min(pki+plat_len,lastindex(arr)-1), plat_len
            end
            arr[pki:min(pki+plat_len,lastindex(arr)-1)] .= 0
        end
    end

    for pklen in 3:96
        arr = repeat([0,1]; outer=cld(pklen,2))
        @test @views argmaxima(arr[1:pklen]) == simplemaxima(arr[1:pklen])
        arr = repeat([1,0]; outer=cld(pklen,2))
        @test @views argmaxima(arr[1:pklen]) == simplemaxima(arr[1:pklen])
    end
end
