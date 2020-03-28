export maxima,
       minima

for (funcname, comp) in ((:maxima, :<),
                         (:minima, :>))

    @eval begin
        function ($funcname)(x::AbstractVector,
                        w::Integer=1,
                        includebounds::Bool=false)
            w > 0 || throw(ArgumentError("window cannot be negative"))
            xlen = length(x)
            idxs = Int[]

            # There can't be more than one peak every `w` elements
            maxN = includebounds ? cld(xlen,2*w) : fld(xlen,2*w)
            N = 0

            ww = [ collect(-w:-1); collect(1:w) ]

            sizehint!(idxs,maxN)

            firsti = firstindex(x)
            lasti = lastindex(x)
            i = firsti
            maxI = lasti-w

            if includebounds
                @inbounds while i <= firsti + w
                    peak = true
                    if ismissing(x[i])
                        i = something(findnext(!ismissing, x, i+1), lasti+1)
                    elseif isnan(x[i])
                        i = something(findnext(!isnan, x, i+1), lasti+1)
                    else
                        for j in max(i-w, firsti):(i+w)
                            i === j && continue
                            if coalesce(($comp)(x[i], x[j]), true) || isnan(x[j])
                                peak &= false
                                break # No reason to continue checking if the rest of the elements qualify
                            elseif x[i] === x[j] && j > i
                                _i = findnext(y -> x[i] !== y, x, j+1)
                                if isnothing(_i) # x is constant till the end, not a peak
                                    peak &= false
                                    i = lasti+1
                                    break
                                elseif coalesce(($comp)(x[i], x[_i]), true) || isnan(x[j]) # x moves towards a peak
                                    peak &= false
                                    break
                                else # Push new peak here to shift the right number of elements
                                    peak &= false
                                    push!(idxs,i)
                                    N += 1
                                    i = max(_i,i+w)
                                    break
                                end
                            end
                        end

                        if peak
                            push!(idxs,i)
                            N += 1
                            i += w # There can't be another peak for at least `w` more elements
                        end
                        i += 1
                    end
                end
            end

            i = firsti + w

            @inbounds while i <= maxI
                peak = true
                if ismissing(x[i])
                    i = something(findnext(!ismissing, x, i+1), lasti+1)
                elseif isnan(x[i])
                    i = something(findnext(!isnan, x, i+1), lasti+1)
                else
                    for j in ww # For all elements within the window
                        if coalesce(($comp)(x[i], x[i+j]), true) || isnan(x[i+j])
                            peak &= false
                            break # No reason to continue checking if the rest of the elements qualify
                        elseif x[i] === x[i+j] && i+j > i
                            _i = findnext(y -> x[i] !== y, x, i+j+1)
                            if isnothing(_i) # x is constant till the end, not a peak
                                peak &= false
                                i = lasti+1
                                break
                            elseif coalesce(($comp)(x[i], x[_i]), true) || isnan(x[_i])
                                # x moves towards a peak or ends with a missing or NaN
                                peak &= false
                                break
                            else # Push first element of plateau as peak here to shift the correct number of elements
                                peak &= false
                                push!(idxs,i)
                                N += 1
                                i = max(_i,i+w)
                                break
                            end
                        end
                    end

                    if peak
                        push!(idxs,i)
                        N += 1
                        i += w # There can't be another peak for at least `w` more elements
                    end
                    i += 1
                end
            end

            if includebounds
                @inbounds while i <= lasti
                    peak = true
                    if ismissing(x[i])
                        i = something(findnext(!ismissing, x, i+1), lasti+1)
                    elseif isnan(x[i])
                        i = something(findnext(!isnan, x, i+1), lasti+1)
                    else
                        for j in (i-w):min((i+w),lasti)
                            i === j && continue
                            if coalesce(($comp)(x[i], x[j]), true)
                                peak &= false
                                break # No reason to continue checking if the rest of the elements qualify
                            elseif x[i] === x[j] && j > i
                                _i = findnext(y -> x[i] !== y, x, j+1)
                                if isnothing(_i) # x is constant till the end, not a peak
                                    peak &= false
                                    i = lasti+1
                                    break
                                elseif coalesce(($comp)(x[i], x[_i]), true) || isnan(x[j]) # x moves towards a peak
                                    peak &= false
                                    break
                                else # Push new peak here to shift the right number of elements
                                    peak &= false
                                    push!(idxs,i)
                                    N += 1
                                    i = max(_i,i+w)
                                    break
                                end
                            end
                        end

                        if peak
                            push!(idxs,i)
                            N += 1
                            i += w # There can't be another peak for at least `w` more elements
                        end
                        i += 1
                    end
                end
            end

            resize!(idxs,N)
            sizehint!(idxs,N)

            return idxs
        end
    end
end

@doc """
    maxima(x[, w=1, includebounds=false])

Find the local maxima of `x`.

`w` sets the minimum allowed distance between maxima. If `includebounds` is `true`, maxima
are allowed to be less than `w` away from the bounds of `x`, otherwise, maxima may not be
located any closer to the ends of `x` than `w`.
"""
maxima

@doc """
    minima(x[, w=1, includebounds=false])

Find the local minima of `x`.

`w` sets the minimum allowed distance between minima. If `includebounds` is `true`, minima
are allowed to be less than `w` away from the bounds of `x`, otherwise, minima may not be
located any closer to the ends of `x` than `w`.
"""
minima

