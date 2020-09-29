for (funcname, comp, notcomp, val, TT) in ((:argmaxima, :<, :>, :i, Int),
    (:argminima, :>, :<, :i, Int))
    # FIXME: Enable maxima/minima => values after deprecation cycle
    # (:maxima, :<, :>, :(x[i]), :T),
    # (:minima, :>, :<, :(x[i]), :T))

    @eval begin
        function ($funcname)(x::AbstractVector{T},
                        w::Integer=1;
                        strictbounds::Bool=true) where T
            w > 0 || throw(ArgumentError("window cannot be negative"))
            xlen = length(x)
            out = Base.nonmissingtype($TT)[]

            # There can't be more than one peak every `w` elements, but a peak is an element as well
            maxN = strictbounds ? max(0,fld(xlen-w,w+1)) : cld(xlen,w+1)
            N = 0

            ww = [ collect(-w:-1); collect(1:w) ]

            sizehint!(out,maxN)

            firsti = firstindex(x)
            lasti = lastindex(x)
            i = firsti
            maxI = lasti - w

            if !strictbounds
                while i <= firsti + w
                    xi = x[i]
                    peak = true
                    if ismissing(xi) || isnan(xi)
                        i = something(findnext(x -> !isnan(x) & !ismissing(x), x, i+1), lasti+1)
                    else
                        for j in firsti:min(i+w,lasti)
                            i === j && continue
                            xj = x[j]
                            if coalesce(($comp)(xi, xj), strictbounds) || (strictbounds && isnan(xj))
                                peak &= false
                                break # No reason to continue checking if the rest of the elements qualify
                            elseif xi === xj
                                k = findnext(y -> xi !== y, x, j+1)
                                if isnothing(k) # x is constant till the end, not a peak
                                    push!(out, $val)
                                    N += 1
                                    i = lasti+1
                                    peak &= false
                                    break
                                else
                                    xk = x[k]
                                    if coalesce(($comp)(xi, xk), strictbounds) # x moves towards a peak
                                        peak &= false
                                        break
                                    else
                                        if k < i + w # if the plateau is within w there could be a larger peak afterwards
                                            j = k
                                        else # Push new peak here to shift the right number of elements
                                            push!(out, $val)
                                            N += 1
                                            i = k
                                            peak &= false
                                            break
                                        end
                                    end
                                end
                            end
                        end

                        if peak
                            push!(out, $val)
                            N += 1
                            i += w # There can't be another peak for at least `w` more elements
                        end
                        i += 1
                    end
                end
            end

            i = max(i, firsti + w, something(findfirst(y -> y !== x[firsti], x), firsti))

            while i <= maxI
                xi = x[i]
                peak = true
                if ismissing(xi) || isnan(xi)
                    i = something(findnext(x -> !isnan(x) & !ismissing(x), x, i+1), lasti+1)
                else
                    for j in ww # For all elements within the window
                        xj = x[i+j]
                        if coalesce(($comp)(xi, xj), strictbounds) || (strictbounds && isnan(xj))
                            peak &= false
                            break # No reason to continue checking if the rest of the elements qualify
                        elseif xi === xj
                            xi1 = x[i-1]
                            if coalesce(($notcomp)(xi, xi1), !strictbounds) || (!strictbounds && isnan(xi1))
                                k = findnext(y -> xi !== y, x, i+j+1)
                                if isnothing(k)
                                    if !strictbounds # x is constant till the end, only a peak for strictbounds false
                                        push!(out, $val)
                                        N += 1
                                        i = lasti+1
                                        peak &= false
                                        break
                                    else
                                        i = lasti+1
                                        peak &= false
                                        break
                                    end
                                else
                                    xk = x[k]
                                    if coalesce(($comp)(xi, xk), strictbounds) || (strictbounds && isnan(xk))
                                        # x moves towards a peak or ends with a missing or NaN
                                        peak &= false
                                        break
                                    elseif j === 1 # Push first element of plateau as peak here to shift the correct number of elements
                                        if k < i + w
                                            j = k - i
                                        else
                                            push!(out, $val)
                                            N += 1
                                            i = k
                                            peak &= false
                                            break
                                        end
                                    end
                                end
                            else
                                peak &= false
                                break
                            end
                        end
                    end

                    if peak
                        push!(out, $val)
                        N += 1
                        i += w # There can't be another peak for at least `w` more elements
                    end
                    i += 1
                end
            end

            if !strictbounds
                while i <= lasti
                    xi = x[i]
                    peak = true
                    if ismissing(xi) || isnan(xi)
                        i = something(findnext(x -> !isnan(x) & !ismissing(x), x, i+1), lasti+1)
                    else
                        for j in (i-w):lasti
                            i === j && continue
                            xj = x[j]
                            if coalesce(($comp)(xi, xj), strictbounds) || (strictbounds && isnan(xj))
                                peak &= false
                                break # No reason to continue checking if the rest of the elements qualify
                            elseif xi === xj
                                xi1 = x[i-1]
                                if coalesce(($notcomp)(xi, xi1), !strictbounds) || (!strictbounds && isnan(xi1))
                                    k = findnext(y -> xi !== y, x, j+1)
                                    if isnothing(k)
                                        push!(out, $val)
                                        N += 1
                                        i = lasti+1
                                        peak &= false
                                        break
                                    else
                                        xk = x[k]
                                        if coalesce(($comp)(xi, xk), strictbounds) || (strictbounds && isnan(xk)) # x moves towards a peak
                                            peak &= false
                                            break
                                        else # Push new peak here to shift the right number of elements
                                            if k < lasti
                                                j = k
                                            else
                                                push!(out, $val)
                                                N += 1
                                                i = max(k,i+w)
                                                peak &= false
                                                break
                                            end
                                        end
                                    end
                                else
                                    peak &= false
                                    break
                                end
                            end
                        end

                        if peak
                            push!(out, $val)
                            N += 1
                            i += w # There can't be another peak for at least `w` more elements
                        end
                        i += 1
                    end
                end
            end

            resize!(out,N)
            sizehint!(out,N)

            return out
        end
    end
end

@doc """
    argmaxima(x[, w=1; strictbounds=true])

Find the indices of the local maxima of `x` where each maxima is either the maximum of
`x[-w:w]` or the first index of a plateau.

If `strictbounds` is `true`, all elements of `x[-w:w]` must exist and may not be `missing`
or `NaN`. If `strictbounds` is `false`, elements of `x[-w:w]` may not exist (eg peaks may
be less than `w` indices from either end of `x`), or may be `missing` or `NaN`.
"""
argmaxima

function findmaxima(x, w::Integer=1; strictbounds::Bool=true)
    idxs = argmaxima(x, w; strictbounds=strictbounds)
    return (idxs, x[idxs])
end

@doc """
    argminima(x[, w=1; strictbounds=false])

Find the indices of the local minima of `x` where each minima is either the minimum of
`x[-w:w]` or the first index of a plateau.

If `strictbounds` is `true`, all elements of `x[-w:w]` must exist and may not be `missing`
or `NaN`. If `strictbounds` is `false`, elements of `x[-w:w]` may not exist (eg peaks may
be less than `w` indices from either end of `x`), or may be `missing` or `NaN`.
"""
argminima

function findminima(x, w::Integer=1; strictbounds::Bool=true)
    idxs = argminima(x, w; strictbounds=strictbounds)
    return (idxs, x[idxs])
end

# Deprecations
@deprecate maxima(x, w=1, strictbounds=true) argmaxima(x, w; strictbounds=strictbounds)
@deprecate minima(x, w=1, strictbounds=true) argminima(x, w; strictbounds=strictbounds)
