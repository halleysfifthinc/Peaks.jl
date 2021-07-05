"""
    peakprom(peaks, x; strictbounds, minprom, maxprom) -> (peaks, proms)

Calculate the prominences of `peaks` in `x`, returning peaks and prominences between `minprom` and `maxprom`.

Peak prominence is calculated as the distance between the current peak and the larger of
the extremums of the left and right bounding intervals. Bounding intervals extend from the
next/previous index from the current peak to the first element larger than or equal to
the current peak, or the end of the signal, whichever comes first.

When bounding intervals contain `NaN` or `missing`, the reported prominence for that peak
will be `NaN` or `missing`, respectively, if `strictbounds = true`; if `strictbounds = false`,
the reference level for the contaminated bounding interval will be the maximum non-`NaN` or
`missing` value.

See also: [`argminima`](@ref), [`argmaxima`](@ref), [`findminima`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = rand(1000);

julia> ma, pa = peakprom(x);

julia> mi, pi = peakprom(Minima(), -x);

julia> @assert (mi == ma) && (pa == pi)
```
"""
function peakprom(peaks::AbstractVector{Int}, x::AbstractVector{T};
    strictbounds=true, minprom=nothing, maxprom=nothing
) where T
    if !isnothing(minprom) || !isnothing(maxprom)
        _peaks = copy(peaks)
    else
        # peaks will not be modified
        _peaks = peaks
    end
    return peakprom!(_peaks, x; strictbounds, minprom, maxprom)
end

"""
    peakprom!(peaks, x; strictbounds, minprom, maxprom) -> (peaks, proms)

Calculate the prominences of `peaks` in `x` and delete peaks with prominences outside of
`minprom` and `maxprom`. Returns the modified `peaks` and the prominences.

See also: [`peakprom`](@ref),[`argminima`](@ref), [`argmaxima`](@ref), [`findminima`](@ref), [`findmaxima`](@ref)
"""
function peakprom!(peaks::AbstractVector{Int}, x::AbstractVector{T};
    strictbounds=true, minprom=nothing, maxprom=nothing
) where T
    proms = similar(peaks,Union{Missing,T})
    isempty(peaks) && return peaks, proms

    fp = first(peaks)
    if fp > 1 && ((x[fp] < x[fp-1]) === true)
        pktype = :minima
    else
        pktype = :maxima
    end
    cmp = pktype === :maxima ? (≥) : (≤)
    exm = pktype === :maxima ? minimum : maximum
    exa = pktype === :maxima ? max : min

    if !strictbounds
        # The extremum search space in the bounding intervals can be reduced by
        # restricting the search space to known peaks/reverse peaks. The cost of
        # finding all peaks/reverse peaks should be mitigated by the fact that
        # the same peaks/reverse peaks will be the pivotal elements for
        # numerous peaks.
        if pktype === :maxima
            peaks′ = argmaxima(x, 1; strictbounds=false)
            notm = argminima(x, 1; strictbounds=false)
        else
            peaks′ = argminima(x, 1; strictbounds=false)
            notm = argmaxima(x, 1; strictbounds=false)
        end
    end

    if strictbounds
        lbegin, lend = firstindex(x), lastindex(x)

        @inbounds for i in eachindex(peaks, proms)
            # Find left and right bound (self-intersections)
            lb = something(findprev(y -> cmp(y, x[peaks[i]]) === true, x, peaks[i] - 2),
                lbegin)
            rb = something(findnext(y -> cmp(y, x[peaks[i]]) === true, x, peaks[i] + 2),
                lend)

            # Find extremum of left and right bounds
            if isempty(lb:(peaks[i] - 1))
                lref = missing
            else
                lref = exm(view(x, lb:(peaks[i] - 1)))
            end

            if isempty((peaks[i] + 1):rb)
                rref = missing
            else
                rref = exm(view(x, (peaks[i] + 1):rb))
            end

            proms[i] = abs(x[peaks[i]] - exa(lref, rref))
        end
    else
        peaks′val = x[peaks′]
        notmval = x[notm]

        j = something(findfirst(y -> y === peaks[1], peaks′), firstindex(peaks′))
        k = something(findfirst(>(peaks[1]), notm), firstindex(notm))
        @inbounds for i in eachindex(peaks, proms)
            j = something(findnext(y -> y === peaks[i], peaks′, j), j)
            k = something(findnext(>(peaks[i]), notm, k), lastindex(notm))

            # Find left and right bounding peaks
            _lb = something(findprev(y -> cmp(y, peaks′val[j]) === true, peaks′val, j - 1), firstindex(peaks′))
            _rb = something(findnext(y -> cmp(y, peaks′val[j]) === true, peaks′val, j + 1), lastindex(peaks′))

            # Find left and right reverse peaks just inside the bounding peaks
            lb = something(findprev(<(peaks′[_lb]+2), notm, k-1), firstindex(notm))
            rb = something(findnext(>(peaks′[_rb]-2), notm, k), lastindex(notm))

            if isempty(lb:(k-1))
                lref = missing
            else
                lref = exm(view(notmval, lb:(k - 1)))
            end

            if isempty(k:rb)
                rref = missing
            else
                rref = exm(view(notmval, k:rb)) # Slice corollary upper side
            end

            proms[i] = abs(x[peaks[i]] - exa(coalesce(lref, rref), coalesce(rref, lref)))
        end
    end

    if !isnothing(minprom) || !isnothing(maxprom)
        lo = something(minprom, zero(eltype(x)))
        up = something(maxprom, typemax(nonmissingtype(eltype(x))))
        matched = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), proms)
        deleteat!(peaks, matched)
        deleteat!(proms, matched)
    end

    return peaks, proms
end

struct Maxima; end
struct Minima; end

@deprecate peakprom(x::AbstractVector, w::Int=1; strictbounds=true, minprom=nothing) peakprom(argmaxima(x, w; strictbounds), x; strictbounds, minprom)
@deprecate peakprom(m::Minima, x::AbstractVector, w::Int=1; strictbounds=true, minprom=nothing) peakprom(argminima(x, w; strictbounds), x; strictbounds, minprom)
@deprecate peakprom(m::Maxima, x::AbstractVector, w::Int=1; strictbounds=true, minprom=nothing) peakprom(argmaxima(x, w; strictbounds), x; strictbounds, minprom)

