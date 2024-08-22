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

function _simpleextrema(@nospecialize(f), cmp::F, x::AbstractVector{T}) where {F,T}
    T >: Missing && throw(MethodError(f, Tuple{typeof(x)}))

    if typeof(axes(x,1)) <: AbstractUnitRange && length(x) > 25
        if T <: SIMD.VecTypes && hasmethod(vload, Tuple{Type{Vec{8,T}}, typeof(x), Int})
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

        if cmp(xim1, xi) # pre
            xip1 = x[i+1]
            if cmp(xip1, xi) # post
                if peaktype === Val(:packed)
                    pks[i-ioffs] = true
                else
                    push!(pks, i)
                end
                i += 2
                continue # skip i += 1 at the end of the loop
            elseif xip1 == xi # plateau
                j = @something findnext(Base.Fix2(!=, xi), x, i+2) lasti+1
                if j ≤ lasti && cmp(x[j], xi) # post
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

Return a bitmask where the lowest bit in each series of consecutive, non-zero bits is set.

# Examples
```julia-repl
julia> lowest_set_bits(0b01110101) |> bitstring
"0b00010101"
```
"""
function lowest_set_bits(x)
    return x & ~(x << 1)
end

"""
    highest_set_bits(x)

Return a bitmask where the highest bit in each series of consecutive, non-zero bits is set.

# Examples
```julia-repl
julia> lowest_set_bits(0b01110101) |> bitstring
"0b01000101"
```
"""
function highest_set_bits(x)
    return x & ~(x >> 1)
end

"""
    matching_bit_runs_mask_lowest_bit(bits, mask)

Return an integer where bits are zeroed unless the lowest set bit of a bit-run (consecutive
non-zero bits) is also set in the bit mask.

The bit mask must be a partial subset of bits set in `lowest_set_bits(bits)`.
(i.e. `iszero(mask ⊻ lowest_set_bits(bits)) && mask ≤ lowest_set_bits(bits)`)
It is the callers responsibility to ensure that `mask` is correct, and this function is not
correct or tested for incorrect masks.

# Examples
```julia-repl
julia> matching_bit_runs_mask_lowest_bit(0b01110101, 0b00010001) |> bitstring
"0b01110001"
```
"""
function matching_bit_runs_mask_lowest_bit(bits, mask)
    matches = lowest_set_bits(bits) & mask
    return bits & (~(bits + matches))
end

"""
    matching_bit_runs_mask_highest_bit(bits, mask)

Return an integer where bits are zeroed unless the highest set bit of a bit-run (consecutive
non-zero bits) is also set in the bit mask.

The bit mask must be a subset of bits set in `highest_set_bits(bits)`.
(i.e. `iszero(mask ⊻ highest_set_bits(bits)) && mask ≤ highest_set_bits(bits)`)
It is the callers responsibility to ensure that `mask` is correct, and this function is not
correct or tested for incorrect masks.

# Examples
```julia-repl
julia> matching_bit_runs_mask_highest_bit(0b01110101, 0b01000001) |> bitstring
"0b01110001"
```
"""
function matching_bit_runs_mask_highest_bit(bits, mask)
    # TODO: Figure out how to do this directly. ~7% of simd_extrema runtime spent here (the
    # bitreverse's specifically)
    matches = highest_set_bits(bits) & mask
    return bitreverse(matching_bit_runs_mask_lowest_bit(bitreverse(bits), bitreverse(matches)))
end

top_bit_set(x) = top_bit_set(typeof(x), x)
function top_bit_set(::Type{T}, x) where T
    return !iszero(x & typemin(signed(T)))
end

function _simd_extrema!(pks::BitVector, cmp::F, x::AbstractVector{T}) where {F,T <: SIMD.VecTypes}
    # Fear this hideous monstrosity...

    lasti = lastindex(x)

    continuing_plat = false
    plat_begin = 0

    if length(x) > 10
        j = 0
        i = firstindex(x)-1

        _c = UInt64(0)
        _plat = UInt64(0)
        _pre = UInt64(0)
        _post = UInt64(0)
        plat = UInt64(0)
        pre = UInt64(0)
        @inbounds while i < lasti-10
            pk = UInt64(0)
            post = UInt64(0)
            shift = 1
            # in 64 elements blocks, create bitmasks for `x[i-1] < x[i]` (aka `pre`),
            # `x[i+1] < x[i]` (aka `post`) and `x[i+1] == x[i]` (aka `plat`)
            for _ in 1:8
                # @debug "" i, j, shift
                xpre = vload(Vec{8,T}, x, i+1)
                xcurr = vload(Vec{8,T}, x, i+2)
                xpost = vload(Vec{8,T}, x, i+3)

                _pre = sum(convert(Vec{8,UInt64}, cmp(xpre, xcurr)) << Vec((0,1,2,3,4,5,6,7)))
                _post = sum(convert(Vec{8,UInt64}, cmp(xpost, xcurr)) << Vec((0,1,2,3,4,5,6,7)))
                _c = _pre & _post

                _plat = sum(convert(Vec{8,UInt64}, xpost == xcurr) << Vec((0,1,2,3,4,5,6,7)))

                pk |= _c << shift
                plat |= _plat << shift
                pre |= _pre << shift

                # `post` is only kept around for verifying plateau ends, where the *next*
                # comparison after the last plateau element must be true. So although the -1
                # shift is now incorrect/mismatched to compare against `pre`, it is correct
                # to compare against `plat`, which is the only remaining use for it. (Since
                # we already compared `_pre` & `_post`)
                post |= _post << (shift-1)

                j += 8
                i += 8
                shift += 8
                i+11 < lasti || break
            end
            pks_j, r = divrem(j, 64)
            pks_j += r > 0
            # @debug "" plat, post, pre, _post, _plat

            # (come back to the comments on this section after reading the next commented
            # sections)
            if continuing_plat && r == 0
                # The plateau began before this chunk so at least the first (lowest) bit of
                # `plat` must be set
                t1s = trailing_ones(plat)
                @assert t1s > 0

                # if `t1s == 64` then the plateau does not end in this chunk
                if t1s < 64
                    plat_end_mask = (2^(t1s-1))
                    # The top bit of the (first) plateau in `plat` must match a (set) bit in
                    # `post` for the plateau to be confirmed
                    if !iszero(post & plat_end_mask)
                        pks[plat_begin] = true

                        # Remove the plateau from consideration in the forthcoming ops to
                        # create this chunk's `pk` (by unsetting the relevent `post` bit)
                        post ⊻= plat_end_mask
                        # @debug "(confirmed) Plateau ends at $(plat_begin+t1s)"
                    end

                    # else this isn't a plateau; either way, reset `continuing_plat`
                    continuing_plat = false
                end
            end

            # plateaus must begin with a pre bit (i.e. a plateau must begin with an element
            # less than the plateau value)
            plat = matching_bit_runs_mask_lowest_bit(plat, pre)

            # plateaus must end with a post bit, but the plateau beginning (i.e. lowest bit
            # in the run) is considered the peak location
            pk |= lowest_set_bits(matching_bit_runs_mask_highest_bit(plat, post))

            # Set the chunk of the bitvector (using OR because if the top bit of _c is set,
            # we need to set the first bit of the next chunk. This is needed because the
            # shifting is always +1 of a byte range, so the top bit of `_c` is always
            # shifted beyond `pk` (i.e. non-existent bit 65).
            pks.chunks[pks_j] |= pk
            # @debug "" r, _c
            if r == 0 && top_bit_set(UInt8, _c)
                pks.chunks[pks_j+1] |= 1
            end

            # @debug "" !continuing_plat, top_bit_set(UInt8, _plat), !top_bit_set(post)

            # A plateau might begin in one chunk and not end until later. We only begin a
            # `continuing_plat` if we're not already in the middle of one. The last element
            # of this chunk must be a plateau and the plateau must not end in this chunk
            # (i.e. top bit of `post` is not set).
            # `_plat` and not `_post` because `_post`/`post` aren't shifted + 1 (so the top
            # bit of `_post` isn't shifted beyond `post`, unlike the other bitmasks)
            if !continuing_plat && top_bit_set(UInt8, _plat) && !top_bit_set(post)
                l1s = leading_ones(plat)
                # @debug "" _plat, _pre l1s > 0, top_bit_set(UInt8, _plat & _pre)

                # `plat` has leading ones, or the top-bit of `_plat` (aka bit 65 of `plat)
                # is set. We store the index of the beginning of the plateau to set (or not)
                # that index in `pks` if this plateau is confirmed
                if l1s > 0 || top_bit_set(UInt8, _plat & _pre)
                    continuing_plat = true
                    plat_begin = j-leading_ones(plat)+1
                    # @debug "Plateau begins @" plat_begin

                # This was an incomplete chunk; we calculate the index of the beginning of
                # the plateau by counting the leading zeros after zero'ing the upper bits of
                # the plateau (if multiple elements long)
                # We still set `continuing_plat` and `plat_begin` for the last, non-SIMD
                # peak finding
                elseif r > 0
                    l0z = leading_zeros(lowest_set_bits(plat))
                    continuing_plat = true
                    plat_begin = j-l0z+8 # plus 8 because idky but it's correct
                end
            end

            # Shift in the top bits of `_plat` and `_pre` which were not resolved in this
            # chunk
            # (`_plat` and `_pre` are UInt64s but only the lower 8 bits might be filled (and
            # only the 8th bit is intended to be kept)
            plat = _plat >> 7
            pre = _pre >> 7
        end

        # @debug "Ending plateau status" continuing_plat, plat_begin
        j = ifelse(continuing_plat,plat_begin,j-8)
        i = ifelse(continuing_plat,plat_begin-(firstindex(x)-1),i-8)
        # @debug j, i
    else
        j = 2
        i = firstindex(x)+1
    end

    # We now return to our regularly scheduled (intelligible) programming
    @inbounds while i < lasti
        xi = x[i]
        xim1 = x[i-1]

        if cmp(xim1, xi) # pre
            xip1 = x[i+1]
            if cmp(xip1, xi) # post
                pks[j] = true
                i += 2
                j += 2
                continue # skip i += 1 at the end of the loop
            elseif xip1 == xi # plateau
                k = @something findnext(Base.Fix2(!=, xi), x, i+2) lasti+1
                # @debug "" i,j,k k ≤ lasti, cmp(x[k], xi)
                if k ≤ lasti && cmp(x[k], xi) # post
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

