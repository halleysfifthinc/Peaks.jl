_smallwerr(w) = throw(DomainError(w,"window size `w` must be greater than zero"))

function findnextextrema(cmp, x::AbstractVector, i::Int, w::Int, strict::Bool)
    w > 0 || _smallwerr(w)
    strict && i - w < firstindex(x) && (i = w + 1) # First peak can't be closer than w
    maxlast = strict ? lastindex(x) - w : lastindex(x)

    while i ≤ maxlast
        peak = true
        xi = x[i]
        if ismissing(xi) || isnan(xi)
            i = something(findnext(p -> !isnan(p) & !ismissing(p), x, i+1), lastindex(x)+1)
        else
            for j in max(firstindex(x),i-w):min(i+w,lastindex(x))
                j === i && continue
                xj = x[j]
                if coalesce(cmp(xi, xj), strict) || (strict && isnan(xj))
                    peak &= false
                    i = j > i ? j : i + 1
                    break # No reason to continue checking if the rest of the elements qualify
                elseif xi === xj # plateau
                    k = findnext(p -> xi !== p, x, j+1)
                    # @show k i j
                    if strict
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
    findnextmaxima(x, i[, w=1; strict=true])

Find the index of the next maxima in `x` after or including `i`, where the maxima `i` is
either the maximum of `x[i-w:i+w]` or the first index of a plateau. Returns `lastindex(x) +
1` if no maxima occur after `i`.

If `strict` is `true`, all elements in `x[i-w:i+w]` or the bounds of a plateau must
exist and must not be `missing` or `NaN`. For `strict == false`, a maxima is the
maximum or first element of all existing and non-NaN or missing elements in `x[i-w:i+w]` or
the bounds of a plateau.

See also: [`argmaxima`](@ref)

# Examples
```jldoctest
julia> findnextmaxima([0,2,0,1,1,0], 2)
2

julia> findnextmaxima([0,2,0,1,1,0], 3)
4

```
"""
findnextmaxima(x, i, w=1; strict=true) = findnextextrema(<, x, i, w, strict)

"""
    argmaxima(x[, w=1; strict=true])

Find the indices of the local maxima of `x` where each maxima `i` is either the maximum of
`x[i-w:i+w]` or the first index of a plateau.

If `strict` is `true`, all elements of `x[i-w:i+w]` or the bounds of a plateau must
exist and must not be `missing` or `NaN`. For `strict == false`, a maxima is the
maximum or first element of all existing and non-NaN or missing elements in `x[i-w:i+w]` or
the bounds of a plateau.

See also: [`findmaxima`](@ref), [`findnextmaxima`](@ref)

# Examples
```jldoctest
julia> argmaxima([0,2,0,1,1,0])
2-element Vector{Int64}:
 2
 4

julia> argmaxima([2,0,1,1])
Int64[]

julia> argmaxima([2,0,1,1]; strict=false)
2-element Vector{Int64}:
 1
 3
```
"""
function argmaxima(
    x::AbstractVector{T}, w::Int=1; strict::Bool=true, strictbounds=nothing
) where T
    if !isnothing(strictbounds)
        Base.depwarn("Keyword `strictbounds` has been renamed to `strict`", :argmaxima)
        strict=strictbounds
    end

    w > 0 || _smallwerr(w)
    pks = Int[]

    i = strict ? w + firstindex(x) : firstindex(x)
    i = findnextmaxima(x, i, w; strict=strict)
    while i ≤ lastindex(x)
        push!(pks, i)
        i = findnextmaxima(x, i+w+1, w; strict=strict)
    end

    return pks
end

"""
    findmaxima(x[, w=1; strict=true]) -> (idxs, vals)

Find the indices and values of the local maxima of `x` where each maxima `i` is either the maximum of
`x[i-w:i+w]` or the first index of a plateau.

See also: [`argmaxima`](@ref), [`findnextmaxima`](@ref)
"""
function findmaxima(x, w::Int=1; strict::Bool=true, strictbounds=nothing)
    if !isnothing(strictbounds)
        Base.depwarn("Keyword `strictbounds` has been renamed to `strict`", :findmaxima)
        strict=strictbounds
    end

    idxs = argmaxima(x, w; strict=strict)
    return (idxs, x[idxs])
end

"""
    findnextminima(x, i[, w=1, strict=true])

Find the index of the next minima in `x` after or including `i`, where the minima `i` is
either the minimum of `x[i-w:i+w]` or the first index of a plateau. Returns `lastindex(x) +
1` if no minima occur after `i`.

If `strict` is `true`, all elements in `x[i-w:i+w]` or the bounds of a plateau must
exist and must not be `missing` or `NaN`. For `strict == false`, a minima is the
minimum or first element of all existing and non-NaN or missing elements in `x[i-w:i+w]` or
the bounds of a plateau.

See also: [`argminima`](@ref)

# Examples
```jldoctest
julia> findnextminima([3,2,3,1,1,3], 2)
2

julia> findnextminima([3,2,3,1,1,3], 3)
4

```
"""
findnextminima(x, i, w=1; strict=true) = findnextextrema(>, x, i, w, strict)

"""
    argminima(x[, w=1; strict=false])


Find the indices of the local minima of `x` where each minima `i` is either the minimum of
`x[i-w:i+w]` or the first index of a plateau.

If `strict` is `true`, all elements of `x[i-w:i+w]` or the bounds of a plateau must
exist and must not be `missing` or `NaN`. For `strict == false`, a minima is the
minimum or first element of all existing and non-NaN or missing elements in `x[i-w:i+w]` or
the bounds of a plateau.

See also: [`findminima`](@ref), [`findnextminima`](@ref)

# Examples
```jldoctest
julia> argminima([3,2,3,1,1,3])
2-element Vector{Int64}:
 2
 4

julia> argminima([2,3,1,1])
Int64[]

julia> argminima([2,3,1,1]; strict=false)
2-element Vector{Int64}:
 1
 3
```
"""
function argminima(
    x::AbstractVector{T}, w::Int=1; strict::Bool=true, strictbounds=nothing
) where T
    if !isnothing(strictbounds)
        Base.depwarn("Keyword `strictbounds` has been renamed to `strict`", :argminima)
        strict=strictbounds
    end

    w > 0 || _smallwerr(w)
    pks = Int[]

    i = strict ? w + firstindex(x) : firstindex(x)
    i = findnextminima(x, i, w; strict=strict)
    while i ≤ lastindex(x)
        push!(pks, i)
        i = findnextminima(x, i+w+1, w; strict=strict)
    end

    return pks
end

"""
    findminima(x[, w=1; strict=true]) -> (idxs, vals)

Find the indices and values of the local minima of `x` where each minima `i` is either the minimum of
`x[i-w:i+w]` or the first index of a plateau.

See also: [`argminima`](@ref), [`findnextminima`](@ref)
"""
function findminima(x, w::Int=1; strict::Bool=true, strictbounds=nothing)
    if !isnothing(strictbounds)
        Base.depwarn("Keyword `strictbounds` has been renamed to `strict`", :findminima)
        strict=strictbounds
    end

    idxs = argminima(x, w; strict=strict)
    return (idxs, x[idxs])
end

# Deprecations
@deprecate maxima(x::AbstractVector, w::Int=1, strictbounds::Bool=true) argmaxima(x, w; strict=strictbounds)
@deprecate minima(x::AbstractVector, w::Int=1, strictbounds::Bool=true) argminima(x, w; strict=strictbounds)
