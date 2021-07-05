function findnextextrema(cmp, x::AbstractVector, i::Int, w::Int, strictbounds::Bool)
    strictbounds && i - w < firstindex(x) && (i = w + 1) # First peak can't be closer than w
    maxlast = strictbounds ? lastindex(x) - w : lastindex(x)

    while i ≤ maxlast
        peak = true
        xi = x[i]
        if ismissing(xi) || isnan(xi)
            i = something(findnext(p -> !isnan(p) & !ismissing(p), x, i+1), lastindex(x)+1)
        else
            for j in max(firstindex(x),i-w):min(i+w,lastindex(x))
                j === i && continue
                xj = x[j]
                if coalesce(cmp(xi, xj), strictbounds) || (strictbounds && isnan(xj))
                    peak &= false
                    i = j > i ? j : i + 1
                    break # No reason to continue checking if the rest of the elements qualify
                elseif xi === xj # plateau
                    k = findnext(p -> xi !== p, x, j+1)
                    # @show k i j
                    if strictbounds
                        if !isnothing(k) && i-1 ≥ firstindex(x)
                            xi1 = x[i-1]
                            xk = x[k]
                            # @show xi1 xi xk i j k
                            if !ismissing(xi1) && !isnan(xi1) && cmp(xi1, xi) &&
                                !ismissing(xk) && !isnan(xk) && cmp(xk, xi)
                                # plateau begins and ends with elements less than the current
                                return i
                            else
                                peak &= false
                                i = k > i ? k : i+1
                                break
                            end
                        else
                            peak &= false
                            i = isnothing(k) ? lastindex(x)+1 :
                                                    k > i ? k : i + 1
                            break
                        end
                    else
                        xk = isnothing(k) ? NaN : x[k]
                        xi1 = i-1 ≥ firstindex(x) ? x[i-1] : NaN
                        if something(k, lastindex(x)+1) ≥ i+w &&
                            (ismissing(xi1) || isnan(xi1) || cmp(xi1, xi)) &&
                            (ismissing(xk) || isnan(xk) || cmp(xk, xi))
                            return i
                        elseif !isnothing(k)
                            if k > i
                                peak &= false
                                i = k
                                break
                            end
                        else
                            peak &= false
                            i = lastindex(x)+1
                            break
                        end
                    end
                end
            end
            peak && return i
        end
    end

    return lastindex(x) + 1
end

"""
    findnextmaxima(x, i, w, strictbounds)

Find the next maxima in `x` at index `i` or greater. Returns `lastindex(x) + 1` if no maxima
occur between `i` and the end of `x`.
"""
findnextmaxima(x, i, w=1; strictbounds=true) = findnextextrema(<, x, i, w, strictbounds)

"""
    argmaxima(x[, w=1; strictbounds=true])

Find the indices of the local maxima of `x` where each maxima is either the maximum of
`x[-w:w]` or the first index of a plateau.

If `strictbounds` is `true`, all elements of `x[-w:w]` must exist and may not be `missing`
or `NaN`. If `strictbounds` is `false`, elements of `x[-w:w]` may not exist (eg peaks may
be less than `w` indices from either end of `x`), or may be `missing` or `NaN`.
"""
function argmaxima(x::AbstractVector{T}, w::Int=1; strictbounds::Bool=true) where T
    w > 0 || throw(ArgumentError("window cannot be negative"))
    pks = Int[]

    i = strictbounds ? w + firstindex(x) : firstindex(x)
    i = findnextmaxima(x, i, w; strictbounds=strictbounds)
    while i ≤ lastindex(x)
        push!(pks, i)
        i = findnextmaxima(x, i+w+1, w; strictbounds=strictbounds)
    end

    return pks
end

function findmaxima(x, w::Int=1; strictbounds::Bool=true)
    idxs = argmaxima(x, w; strictbounds=strictbounds)
    return (idxs, x[idxs])
end

"""
    findnextminima(x, i, w, strictbounds)

Find the next minima in `x` at index `i` or greater. Returns `lastindex(x) + 1` if no minima
occur between `i` and the end of `x`.
"""
findnextminima(x, i, w=1; strictbounds=true) = findnextextrema(>, x, i, w, strictbounds)

"""
    argminima(x[, w=1; strictbounds=false])

Find the indices of the local minima of `x` where each minima is either the minimum of
`x[-w:w]` or the first index of a plateau.

If `strictbounds` is `true`, all elements of `x[-w:w]` must exist and may not be `missing`
or `NaN`. If `strictbounds` is `false`, elements of `x[-w:w]` may not exist (eg peaks may
be less than `w` indices from either end of `x`), or may be `missing` or `NaN`.
"""
function argminima(x::AbstractVector{T}, w::Int=1; strictbounds::Bool=true) where T
    w > 0 || throw(ArgumentError("window cannot be negative"))
    pks = Int[]

    i = strictbounds ? w + firstindex(x) : firstindex(x)
    i = findnextminima(x, i, w; strictbounds=strictbounds)
    while i ≤ lastindex(x)
        push!(pks, i)
        i = findnextminima(x, i+w+1, w; strictbounds=strictbounds)
    end

    return pks
end

function findminima(x, w::Int=1; strictbounds::Bool=true)
    idxs = argminima(x, w; strictbounds=strictbounds)
    return (idxs, x[idxs])
end

# Deprecations
@deprecate maxima(x::AbstractVector, w::Int=1, strictbounds::Bool=true) argmaxima(x, w; strictbounds=strictbounds)
@deprecate minima(x::AbstractVector, w::Int=1, strictbounds::Bool=true) argminima(x, w; strictbounds=strictbounds)
