# copied and lightly modified from Base.findall(::BitVector)
function findall_offset(B::BitVector, offset::Int)
    nnzB = count(B)
    I = Vector{eltype(keys(B))}(undef, nnzB)
    nnzB == 0 && return I
    Bc = B.chunks
    Bs = size(B)
    Bi = i1 = i = 1
    irest = ntuple(one, ndims(B) - 1)
    c = Bc[1]
    @inbounds while true
        while c == 0
            Bi == length(Bc) && return I
            i1 += 64
            Bi += 1
            c = Bc[Bi]
        end

        tz = trailing_zeros(c)
        c = Base._blsr(c)

        i1, irest = Base._overflowind(i1 + tz, irest, Bs)
        I[i] = Base._toind(i1, irest) - offset
        i += 1
        i1 -= tz
    end
end

_simpleextrema(@nospecialize(f), cmp, x::AbstractVector{>:Missing}) = throw(MethodError(f, Tuple{typeof(x)}))
function _simpleextrema(@nospecialize(f), cmp::F, x::AbstractVector{T}) where {F,T}
    if typeof(axes(x,1)) <: AbstractUnitRange && length(x) > 32
        if T <: VecTypes && hasmethod(vload, Tuple{Type{Vec{8,T}}, typeof(x), Int})
            pks = BitVector(undef, length(x))
            fill!(pks, false)

            _simd_extrema!(pks, cmp, x)
        else
            pks = _simpleextrema_base(cmp, x, Val(:packed))
        end
        return findall_offset(pks, -firstindex(x)+1)
    else
        return _simpleextrema_base(cmp, x, Val(:push))
    end
end

function _simpleextrema_base(cmp::F, x::AbstractVector{T}, peaktype::Val) where {F,T}
    if peaktype === Val(:packed)
        pks = BitVector(undef, length(x))
        fill!(pks, false)
    else
        pks = Int[]
    end
    lasti = lastindex(x)
    ioffs = firstindex(x)-1

    i = firstindex(x) + 1
    @inbounds while i < lasti
        xi = x[i]
        xim1 = x[i-1]

        if cmp(xim1, xi) # rise
            xip1 = x[i+1]
            if cmp(xip1, xi) # fall
                if peaktype === Val(:packed)
                    pks[i-ioffs] = true
                else
                    push!(pks, i)
                end
                i += 2
                continue # skip i += 1 at the end of the loop
            elseif xip1 == xi # plateau
                j = something(findnext(Base.Fix2(!=, xi), x, i+2), lasti+1)
                if j ≤ lasti && cmp(x[j], xi) # fall
                    if peaktype === Val(:packed)
                        pks[i-ioffs] = true
                    else
                        push!(pks, i)
                    end
                    i = j+1
                continue # skip i += 1 at the end of the loop
                end
            end
        end
        i += 1
    end

    if peaktype === Val(:packed)
        return pks
    else
        return pks
    end
end

"""
    lowest_set_bits(x)

Return LSBs of each run of consecutively set bits (including single-bit runs).

# Examples
```jldoctest; setup = :(using Peaks: lowest_set_bits)
julia> lowest_set_bits(0b01110101) |> bitstring
"00010101"
```
"""
function lowest_set_bits(x)
    return x & ~(x << 1)
end

"""
    highest_set_bits(x)

Return MSBs of each run of consecutively set bits (including single-bit runs).

# Examples
```jldoctest; setup = :(using Peaks: highest_set_bits)
julia> highest_set_bits(0b01110101) |> bitstring
"01000101"
```
"""
function highest_set_bits(x)
    return x & ~(x >> 1)
end

"""
    matching_runs_mask_lsb(bits, mask)

Return an integer containing only bit runs whose LSB is set in mask.

The bit mask is coerced to only contain LSB of bit runs (i.e. `mask & lowest_set_bits(bits)`).

# Examples
```jldoctest; setup = :(using Peaks: matching_runs_mask_lsb)
julia> matching_runs_mask_lsb(0b01110101, 0b00010001) |> bitstring
"01110001"
```
"""
function matching_runs_mask_lsb(bits, mask)
    matches = lowest_set_bits(bits) & mask
    return _unchecked_matching_runs_mask_lsb(bits, matches)
end

_unchecked_matching_runs_mask_lsb(bits, mask) = bits & (~(bits + mask))

"""
    matching_runs_mask_msb(bits, mask)

Return an integer containing only bit runs whose MSB is set in mask.

The bit mask is coerced to only contain MSB of bit runs (i.e. `mask & highest_set_bits(bits)`).

# Examples
```jldoctest; setup = :(using Peaks: matching_runs_mask_msb)
julia> matching_runs_mask_msb(0b01110101, 0b01000001) |> bitstring
"01110001"
```
"""
function matching_runs_mask_msb(bits, mask)
    matches = highest_set_bits(bits) & mask
    let ⊕=reversecarry_add
        return bits & (~(bits ⊕ matches))
    end
end

"""
    reversecarry_add(x, y)

Add x and y with reversed (left-to-right) carry semantics.

Based on a [Kogge-Stone adder](https://en.wikipedia.org/wiki/Kogge%E2%80%93Stone_adder).

# Examples
```jldoctest; setup = :(using Peaks: reversecarry_add)
julia> reversecarry_add(0b01110101, 0b01000001) |> bitstring
"00001100"
```
"""
function reversecarry_add(x::T, y::T) where T<:Union{UInt8,UInt16,UInt32,UInt64}
    P = x ⊻ y # bits that will cause carry's to propagate
    G = x & y # additions that generate carry bits

    # First stage: merge adjacent 1-bit runs
    G |= P & (G >> 1) # (max run length in G is 2 at this point)
    P &= P >> 1

    # Second stage: merge adjacent 2-bit runs
    G |= P & (G >> 2) # (max run length in G is 4 at this point)
    P &= P >> 2

    # Third stage: merge adjacent 4-bit runs, etc
    G |= P & (G >> 4)
    if sizeof(T) > 1 # updated P and higher stages aren't needed for UInt8
        P &= P >> 4

        G |= P & (G >> 8)
        if sizeof(T) > 2
            P &= P >> 8

            G |= P & (G >> 16)
            if sizeof(T) > 4
                P &= P >> 16

                G |= P & (G >> 32)
            end
        end
    end

    return x ⊻ y ⊻ (G >> 1)
end

"""
    lsb_of_runs_mask_msb(bits, mask)

Return the LSB of the bit runs selected by matching MSB in the mask.

The bit mask is coerced to only contain MSB of runs (i.e. `mask & highest_set_bits(bits)`).

# Examples
```jldoctest; setup = :(using Peaks: lsb_of_runs_mask_msb)
julia> lsb_of_runs_mask_msb(0b01110101, 0b01000001) |> bitstring
"00010001"
```
"""
function lsb_of_runs_mask_msb(bits::T, mask::T) where T <: Union{UInt8,UInt16,UInt32,UInt64}
    G = mask & highest_set_bits(bits)
    P = bits ⊻ G

    G |= P & (G >> 1)
    P &= P >> 1

    G |= P & (G >> 2)
    P &= P >> 2

    G |= P & (G >> 4)
    if sizeof(T) > 1
        P &= P >> 4

        G |= P & (G >> 8)
        if sizeof(T) > 2
            P &= P >> 8

            G |= P & (G >> 16)
            if sizeof(T) > 4
                P &= P >> 16

                G |= P & (G >> 32)
            end
        end
    end

    return lowest_set_bits(G)
end

top_bit_set(x) = top_bit_set(typeof(x), x)
function top_bit_set(::Type{T}, x) where T
    return !iszero(x & typemin(signed(T)))
end

@noinline _unreachable_assertion() = throw(AssertionError("this is logically unreachable"))

function _simd_extrema!(pks::BitVector, cmp::F, x::AbstractVector{T}) where {F,T <: VecTypes}
    # Fear this hideous monstrosity...

    lasti = lastindex(x)

    if length(x) > 10
        ovrflw_plateau = false

        j = 0
        i = firstindex(x)-1

        carry = UInt64(0)
        ovrflw_plateau_chunk = chunk = 1
        ovrflw_plateau_mask = UInt64(0)

        _rise, _flat = UInt64(0), UInt64(0)
        @inbounds while i < lasti-10
            local _pk, _fall

            shift = (j & 0x3f) + 1 # manual equivalent of (j % 64) + 1

            # full-width comparisons with the top bits of `_flat` and `_rise`
            # that were not resolved in this chunk
            (pk, rise, flat, fall) = carry, _rise >> 7, _flat >> 7, UInt64(0)

            # in 64 elements blocks, create bitmasks for `x[i-1] < x[i]` (aka `rise`),
            # `x[i+1] < x[i]` (aka `fall`) and `x[i+1] == x[i]` (aka `flat`)
            # "rise" and "fall" as in towards/away from a peak (whether maxima or minima)
            for _ in 1:8
                # @debug "" i, j, shift
                x_pre = vload(Vec{8,T}, x, i+1)
                x_curr = vload(Vec{8,T}, x, i+2)
                x_post = vload(Vec{8,T}, x, i+3)

                _rise = UInt64(bitmask(cmp(x_pre, x_curr)))
                _fall = UInt64(bitmask(cmp(x_post, x_curr)))
                _pk = _rise & _fall

                _flat = UInt64(bitmask(x_post == x_curr))

                # shift will never be >= 64; the assert allows the compiler to remove bounds
                # checking for the shift amount
                let shift=shift
                    shift < 64 || _unreachable_assertion()

                    pk |= _pk << shift
                    flat |= _flat << shift
                    rise |= _rise << shift

                    # `fall` is only kept around for verifying plateau ends, where the *next*
                    # comparison after the last plateau element must be true. So although the -1
                    # shift is now incorrect/mismatched to compare against `rise`, it is correct
                    # to compare against `flat`, which is the only remaining use for it. (Since
                    # we already compared `_rise` & `_fall`)
                    fall |= _fall << (shift-1)
                end

                j += 8
                i += 8
                shift += 8
                i+11 < lasti || break
            end
            # equivalent to chunk, r = divrem(j, 64)
            r = j & 0x3f
            chunk = j >> 6

            r_is_zero = r == 0
            chunk += Int(!r_is_zero)
            # @debug "" flat, fall, rise, _fall, _flat

            # (come back to the comments on this section after reading the next commented
            # sections)
            if ovrflw_plateau && r_is_zero # only handle the ovrflw_plateau for a complete chunk
                # The plateau began before this chunk so at least the first (lowest) bit of
                # `flat` must be set
                t1s = trailing_ones(flat)

                # if `t1s == 64` then the plateau does not end in this chunk
                if t1s < 64
                    t1s -= 1
                    # give the compiler enough information to safely elide the shift bounds check
                    (t1s & 63) == t1s  || _unreachable_assertion()
                    plat_end_mask = UInt64(0x1) << t1s

                    # The top bit of the (first) plateau in `flat` must match a (set) bit in
                    # `fall` for the plateau to be confirmed
                    if !iszero(fall & plat_end_mask)
                        pks.chunks[ovrflw_plateau_chunk] |= ovrflw_plateau_mask
                        # @debug "(confirmed) Plateau ends at $(plat_begin+t1s)"
                    end

                    # else this isn't a plateau; either way, reset `continuing_plat`
                    ovrflw_plateau = false
                end
            end

            # skip unnecessary plateau checks
            if !iszero(flat)
                # plateaus must begin with a rise bit (i.e. a plateau must begin with an element
                # less than the plateau value)
                # The only set bits in `rise` are either unset or the LSB of a run in `flat`
                flat = _unchecked_matching_runs_mask_lsb(flat, flat & rise)

                # plateaus must end with a fall bit, but the plateau beginning (i.e. LSB
                # in the run) is considered the peak location
                pk |= lsb_of_runs_mask_msb(flat, fall)
            end

            pks.chunks[chunk] |= pk
            # When `shift == 57`, the MSB of _pk is shifted (overflows) into bit 65.
            # Carry it for bit 0 of the next chunk (where shift=1 leaves bit 0 unused).
            carry = UInt64(r_is_zero & top_bit_set(UInt8, _pk))

            # @debug "" !continuing_plat, top_bit_set(UInt8, _flat), !top_bit_set(fall)

            if !ovrflw_plateau && # can't already be in the middle of an overflowing plateau
                    top_bit_set(UInt8, _flat) && # MUST overflow[^1]
                    !top_bit_set(fall) # MUST NOT have finished in this chunk[^2]
                # [^1]: The top bit of _flat (and _rise) is shifted out of (overflows) the top of chunk
                # [^2]: The MSB of `_flat` is *after* the MSB of `fall`. It's possible for
                # `top_bit_set(UInt8, _flat) & top_bit_set(fall) == true` if a plateau ended
                # in the MSB of `fall` with two consecutively equal values (that are less
                # than the plateau height)

                # Logically, the highest set bit in `rise` corresponds to the LSB of the
                # highest (overflowing) run in `flat`
                l0s = leading_zeros(rise)
                if l0s < 64
                    ovrflw_plateau = true
                    ovrflw_plateau_chunk = chunk

                    # give the compiler enough information to safely elide the shift bounds check
                    (l0s & 63) == l0s  || _unreachable_assertion()
                    ovrflw_plateau_mask = (typemin(Int64) % UInt64) >> l0s
                elseif top_bit_set(UInt8, _rise) && i >= lasti-10
                    # Plateau starts entirely in the overflow bit (bit 65) and we are exiting the SIMD loop
                    # Set ovrflw_plateau* stuff to set up for the scalar loop
                    ovrflw_plateau = true
                    ovrflw_plateau_chunk = chunk + 1
                    ovrflw_plateau_mask = UInt64(0x1)
                end
            end
        end

        # Flush any carry from the last full chunk
        if carry != 0
            pks.chunks[chunk + 1] |= carry
        end

        # @debug "Ending plateau status" continuing_plateau, plat_chunk, plat_mask
        if ovrflw_plateau
            plat_begin = ovrflw_plateau_chunk*64 + trailing_zeros(ovrflw_plateau_mask) - 63
            j = plat_begin
            i = plat_begin + firstindex(x) - 1
        else
            j = j + 2
            i = i + 2
        end
        # @debug j, i
    else
        j = 2
        i = firstindex(x)+1
    end

    # We now return to our regularly scheduled (intelligible) programming
    @inbounds while i < lasti
        xi = x[i]
        xim1 = x[i-1]

        if cmp(xim1, xi) # rise
            xip1 = x[i+1]
            if cmp(xip1, xi) # fall
                pks[j] = true
                i += 2
                j += 2
                continue # skip i += 1 at the end of the loop
            elseif xip1 == xi # plateau
                k = something(findnext(Base.Fix2(!=, xi), x, i+2), lasti+1)
                # @debug "" i,j,k k ≤ lasti, cmp(x[k], xi)
                if k ≤ lasti && cmp(x[k], xi) # fall
                    pks[j] = true
                    j = j+(k-i)+1
                    i = k+1
                continue # skip i += 1 at the end of the loop
                end
            end
        end
        i += 1
        j += 1
    end

    return pks
end

