"""
    peakprom(peaks, x;
        strictbounds=true,
        minprom=nothing,
        maxprom=nothing
    ) -> (peaks, proms)

Calculate the prominences of `peaks` in `x`, filtering peaks and prominences less than
`minprom` and greater than `maxprom`, if either are given.

Peak prominence is the absolute height difference between the current peak and the larger of
the two adjacent smallest magnitude points between the current peak and adjacent larger
peaks or signal ends.

The prominence for a peak with a `NaN` or `missing` between the current peak and either
adjacent larger peaks will be `NaN` or `missing` if `strictbounds == true`, or it will be
the larger of the smallest non-`NaN` or `missing` values between the current peak and
adjacent larger peaks for `strictbounds == false`.

See also: [`argminima`](@ref), [`argmaxima`](@ref), [`findminima`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> xpks = argmaxima(x)
3-element Vector{Int64}:
 2
 4
 7

julia> peakprom(xpks, x)
([2, 4, 7], Union{Missing, Int64}[5, 1, 3])

julia> x = [missing,5,2,3,3,1,4,0];

julia> peakprom(xpks, x)
([2, 4, 7], Union{Missing, Int64}[missing, 1, 3])

julia> peakprom(xpks, x; strictbounds=false)
([2, 4, 7], Union{Missing, Int64}[4, 1, 3])
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
    return peakprom!(_peaks, x; strictbounds=strictbounds, minprom=minprom, maxprom=maxprom)
end

"""
    peakprom!(peaks, x; strictbounds, minprom, maxprom) -> (peaks, proms)

Calculate the prominences of `peaks` in `x`, deleting `peaks` with prominences less than
`minprom` and greater than `maxprom`, if either are given. Returns the modified `peaks` and
their prominences.

See also: [`peakprom`](@ref), [`argminima`](@ref), [`argmaxima`](@ref), [`findminima`](@ref), [`findmaxima`](@ref)
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
        up = something(maxprom, typemax(Base.nonmissingtype(eltype(x))))
        matched = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), proms)
        deleteat!(peaks, matched)
        deleteat!(proms, matched)
    end

    return peaks, proms
end

struct Maxima; end
struct Minima; end

@deprecate peakprom(x::AbstractVector, w::Int=1; strictbounds=true, minprom=nothing) peakprom(argmaxima(x, w; strictbounds=strictbounds), x; strictbounds=strictbounds, minprom=minprom)
@deprecate peakprom(m::Minima, x::AbstractVector, w::Int=1; strictbounds=true, minprom=nothing) peakprom(argminima(x, w; strictbounds=strictbounds), x; strictbounds=strictbounds, minprom=minprom)
@deprecate peakprom(m::Maxima, x::AbstractVector, w::Int=1; strictbounds=true, minprom=nothing) peakprom(argmaxima(x, w; strictbounds=strictbounds), x; strictbounds=strictbounds, minprom=minprom)

