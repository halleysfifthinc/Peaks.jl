"""
    peakproms!(peaks, x;
        strict=true,
        minprom=nothing,
        maxprom=nothing
    ) -> (peaks, proms)

Calculate the prominences of `peaks` in `x`, and removing `peaks` with prominences less than
`minprom` and/or greater than `maxprom`. Returns the modified arrays peaks and their
prominences.

See also: [`peakproms`](@ref), [`findminima`](@ref), [`findmaxima`](@ref)
"""
function peakproms!(peaks::AbstractVector{Int}, x::AbstractVector{T};
    strict=true, minprom=nothing, maxprom=nothing, 
    min=minprom, max=maxprom
) where {T}
    if !isnothing(minprom)
        Base.depwarn("Keyword `minprom` has been renamed to `min`", :peakproms!)
    end
    if !isnothing(maxprom)
        Base.depwarn("Keyword `maxprom` has been renamed to `max`", :peakproms!)
    end
    if !isnothing(min) && !isnothing(max)
        min < max || throw(ArgumentError("minimal prominence must be less than maximal prominence"))
    end
    all(∈(eachindex(x)), peaks) ||
        throw(ArgumentError("peaks contains invalid indices to x"))
    isempty(peaks) && return peaks, T[]

    # if peaks was calculated with strict=false, first(peaks) could be minima at firstindex
    fp = length(peaks) > 1 ? peaks[2] : first(peaks)
    if fp > 1 && ((x[fp] < x[fp-1]) === true)
        pktype = :minima
    else
        pktype = :maxima
    end
    cmp = pktype === :maxima ? (≥) : (≤)
    exm = pktype === :maxima ? minimum : maximum
    exa = pktype === :maxima ? Base.max : Base.min

    _ref = Missing <: T ? missing :
           Float64 <: T ? NaN :
           Float32 <: T ? NaN32 :
           Float16 <: T ? NaN16 :
                          missing

    proms = similar(peaks, promote_type(T, typeof(_ref)))

    if strict
        lbegin, lend = firstindex(x), lastindex(x)

        @inbounds for i in eachindex(peaks, proms)
            # Find left and right bound (self-intersections)
            lb = something(findprev(y -> cmp(y, x[peaks[i]]) === true, x, peaks[i] - 2),
                lbegin)
            rb = something(findnext(y -> cmp(y, x[peaks[i]]) === true, x, peaks[i] + 2),
                lend)

            # Find extremum of left and right bounds
            if isempty(lb:(peaks[i]-1))
                lref = _ref
            else
                lref = exm(view(x, lb:(peaks[i]-1)))
            end

            if isempty((peaks[i]+1):rb)
                rref = _ref
            else
                rref = exm(view(x, (peaks[i]+1):rb))
            end

            proms[i] = abs(x[peaks[i]] - exa(lref, rref))
        end
    else
        # The extremum search space in the bounding intervals can be reduced by
        # restricting the search space to known peaks/reverse peaks. The cost of
        # finding all peaks/reverse peaks should be mitigated by the fact that
        # the same peaks/reverse peaks will be the pivotal elements for
        # numerous peaks.
        if pktype === :maxima
            peaks′ = argmaxima(x, 1; strict=false)
            notm = argminima(x, 1; strict=false)
        else
            peaks′ = argminima(x, 1; strict=false)
            notm = argmaxima(x, 1; strict=false)
        end

        notmval = x[notm]

        for i in eachindex(peaks, proms)
            j = searchsorted(peaks′, peaks[i])

            # Find left and right bounding peaks
            _lb = findprev(y -> cmp(x[y], x[peaks[i]]) === true, peaks′, first(j) - 1)
            peaks′[j] === peaks[i] && (j += 1)
            _rb = findnext(y -> cmp(x[y], x[peaks[i]]) === true, peaks′, last(j) + 1)

            # Find left and right reverse peaks just inside the bounding peaks
            lb = isnothing(_lb) ? firstindex(notm) :
                 searchsortedfirst(notm, peaks′[_lb])
            rb = isnothing(_rb) ? lastindex(notm) :
                 searchsortedlast(notm, peaks′[_rb])

            k = searchsortedfirst(notm, peaks[i])

            if isempty(lb:(k-1))
                lref = missing
            else
                lref = exm(view(notmval, lb:(k-1)))
            end

            rb > lastindex(notm) && (rb -= 1)
            if isempty(k:rb)
                rref = missing
            else
                rref = exm(view(notmval, k:rb)) # Slice corollary upper side
            end

            proms[i] = abs(x[peaks[i]] - exa(coalesce(lref, rref), coalesce(rref, lref)))
        end
    end

    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(x)))
        up = something(max, typemax(Base.nonmissingtype(eltype(x))))
        matched = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), proms)
        deleteat!(peaks, matched)
        deleteat!(proms, matched)
    end

    return peaks, proms
end
export peakproms!

"""
    peakproms(peaks, x;
        strict=true,
        minprom=nothing,
        maxprom=nothing
    ) -> (peaks, proms)

Calculate the prominences of `peaks` in `x`, and removing peaks with prominences less than
`minprom` and/or greater than `maxprom`.

Peak prominence is the absolute height difference between the current peak and the larger of
the two adjacent smallest magnitude points between the current peak and adjacent larger
peaks or signal ends.

The prominence for a peak with a `NaN` or `missing` between the current peak and either
adjacent larger peaks will be `NaN` or `missing` if `strict == true`, or it will be
the larger of the smallest non-`NaN` or `missing` values between the current peak and
adjacent larger peaks for `strict == false`.

See also: [`findminima`](@ref), [`findmaxima`](@ref), [`peakproms!`](@ref)

# Examples
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> xpks = argmaxima(x)
3-element Vector{Int64}:
 2
 4
 7

julia> peakproms(xpks, x)
([2, 4, 7], Union{Missing, Int64}[5, 1, 3])

julia> x = [missing,5,2,3,3,1,4,0];

julia> peakproms(xpks, x)
([2, 4, 7], Union{Missing, Int64}[missing, 1, 3])

julia> peakproms(xpks, x; strict=false)
([2, 4, 7], Union{Missing, Int64}[5, 1, 3])
```
"""
function peakproms(peaks::AbstractVector{Int}, x::AbstractVector{T};
    strict=true, minprom=nothing, maxprom=nothing,
    min=minprom, max=maxprom
) where {T}
    if !isnothing(minprom)
        Base.depwarn("Keyword `minprom` has been renamed to `min`", :peakproms)
    end
    if !isnothing(maxprom)
        Base.depwarn("Keyword `maxprom` has been renamed to `max`", :peakproms)
    end
    if !isnothing(min) || !isnothing(max)
        _peaks = copy(peaks)
    else
        # peaks will not be modified
        _peaks = peaks
    end
    return peakproms!(_peaks, x; strict=strict, min=min, max=max)
end
export peakproms

