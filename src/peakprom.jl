"""
    peakproms(indices, x; [strict=true, min, max]) -> (indices, proms)
    peakproms(pks::NamedTuple; [strict=true, min, max]) -> NamedTuple

Calculate the prominences of peak `indices` in `x`, and remove peaks with prominences less
than `min` and/or greater than `max`.

Peak prominence is the absolute height (value) difference between the current peak and the
larger of the two adjacent smallest magnitude points between the current peak and adjacent
larger peaks or signal ends.

If a NamedTuple `pks` is given, a new NamedTuple is returned with filtered copies of fields
from `pks`. `pks` must have `:indices` and `:heights` fields. If `pks` has a `:proms` field,
prominences will only be filtered, and not be recalculated. The fields `:widths` and
`:edges` will also be filtered if present, and any remaining fields will be copied
unmodified.

If `strict == true`, the prominence for a peak with a `NaN` or `missing` between the current
peak and either adjacent larger peaks will be `NaN` or `missing`, otherwise, it will be the
larger of the smallest non-`NaN` or `missing` values between the current peak and adjacent
larger peaks for `strict == false`.

See also: [`peakproms!`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> pks = findmaxima([0,5,2,3,3,1,4,0]);

julia> pks = peakproms(pks; min=2)
(indices = [2, 7], heights = [5, 4], data = [0, 5, 2, 3, 3, 1, 4, 0], proms = Union{Missing, Int64}[5, 3])

julia> inds, proms = peakproms(pks.indices, pks.data; max=4)
([7], Union{Missing, Int64}[3])
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

peakproms(pks::NamedTuple; kwargs...) = peakproms!(deepcopy(pks); kwargs...)

"""
    peakproms(; [strict, min, max]) -> Function

Create a function, `f(pks::NamedTuple)`, that calculates and filters the peak
prominences of a copy of its argument, `pks`, using any given keyword arguments.

# Examples
```jldoctest
julia> findmaxima([0,5,2,3,3,1,4,0]) |> peakproms(; min=2)
(indices = [2, 7], heights = [5, 4], data = [0, 5, 2, 3, 3, 1, 4, 0], proms = Union{Missing, Int64}[5, 3])
```
"""
peakproms(; kwargs...) = function _curried_peakproms(pks)
    return peakproms(pks; kwargs...)
end

function _strict_inner_promscalcloop!(cmp::C, extremum::M, extrema::A, _ref::MT, x::AbstractVector{T}, peaks::AbstractVector{Int}, proms::AbstractVector{MT}) where {C,M,A,T,MT}
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
                lref = extremum(view(x, lb:(peaks[i]-1)))
            end

            if isempty((peaks[i]+1):rb)
                rref = _ref
            else
                rref = extremum(view(x, (peaks[i]+1):rb))
            end

            proms[i] = abs(x[peaks[i]] - extrema(lref, rref))
        end

        return nothing
end

"""
    peakproms!(indices, x; [strict=true, min, max]) -> (indices, proms)
    peakproms!(pks::NamedTuple; [strict=true, min, max]) -> NamedTuple

Calculate the prominences of peak `indices` in `x`, and remove peaks with prominences less
than `min` and/or greater than `max`.

If a NamedTuple `pks` is given, a new NamedTuple is returned with the same fields
(references) from `pks`. `pks` must have `:indices` and `:heights` fields. If `pks` has a
`:proms` field, prominences will only be filtered, and not be recalculated. The fields
`:widths` and `:edges` will also be filtered (mutated) if present, and any remaining fields
will be copied unmodified.

See also: [`peakproms`](@ref), [`findmaxima`](@ref)
#
# Examples
```jldoctest
julia> pks = findmaxima([0,5,2,3,3,1,4,0]);

julia> pks = peakproms!(pks; min=2)
(indices = [2, 7], heights = [5, 4], data = [0, 5, 2, 3, 3, 1, 4, 0], proms = Union{Missing, Int64}[5, 3])

julia> inds, proms = peakproms!(pks.indices, pks.data; max=4)
([7], Union{Missing, Int64}[3])
```
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
        min < max || throw(ArgumentError("Keyword `min` must be less than `max`"))
    end
    all(∈(eachindex(x)), peaks) ||
        throw(ArgumentError("peaks contains invalid indices to x"))
    isempty(peaks) && return peaks, T[]

    # if peaks was calculated with strict=false, first(peaks) could be minima at firstindex
    if ismaxima(first(peaks), x; strict=false)
        maxima = true
    elseif isminima(first(peaks), x; strict=false)
        maxima = false
    else
        throw(ArgumentError("The first peak in `indices` is not a local extrema"))
    end
    cmp = maxima ? (≥) : (≤)
    exm = maxima ? minimum : maximum
    exa = maxima ? Base.max : Base.min

    _ref = Missing <: T ? missing :
           Float64 <: T ? NaN :
           Float32 <: T ? NaN32 :
           Float16 <: T ? NaN16 :
                          missing

    proms = similar(peaks, promote_type(T, typeof(_ref)))

    if strict
        # Add a function barrier and manually union-split the cmp/exm/exa functions
        if maxima
            _strict_inner_promscalcloop!(≥, minimum, Base.max, _ref, x, peaks, proms)
        else
            _strict_inner_promscalcloop!(≤, maximum, Base.min, _ref, x, peaks, proms)
        end
    else
        # The extremum search space in the bounding intervals can be reduced by
        # restricting the search space to known peaks/reverse peaks. The cost of
        # finding all peaks/reverse peaks should be mitigated by the fact that
        # the same peaks/reverse peaks will be the pivotal elements for
        # numerous peaks.
        if maxima
            peaks′ = argmaxima(x, 1; strict=false)
            notm = argminima(x, 1; strict=false)
        else
            peaks′ = argminima(x, 1; strict=false)
            notm = argmaxima(x, 1; strict=false)
        end

        notmval = x[notm]

        for i in eachindex(peaks, proms)
            j = only(searchsorted(peaks′, peaks[i]))

            # Find left and right bounding peaks
            if maxima # cmp = ≥, manual union-splitting
                _lb = findprev(y -> ≥(x[y], x[peaks[i]]) === true, peaks′, j - 1)
                _rb = findnext(y -> ≥(x[y], x[peaks[i]]) === true, peaks′, j + 1)
            else # cmp = ≤
                _lb = findprev(y -> ≤(x[y], x[peaks[i]]) === true, peaks′, j - 1)
                _rb = findnext(y -> ≤(x[y], x[peaks[i]]) === true, peaks′, j + 1)
            end

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

            # exa(coalesce(lref, rref), coalesce(rref, lref)))
            # we manually union-split this for better type-inference
            exa_ref = if ismissing(lref)
                rref
            elseif ismissing(rref)
                lref
            else
                exa(lref, rref)
            end

            proms[i] =  abs(x[peaks[i]] - exa_ref)
        end
    end

    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(x)))
        hi = something(max, typemax(Base.nonmissingtype(eltype(x))))
        matched = findall(x -> !ismissing(x) && !(lo ≤ x ≤ hi), proms)
        deleteat!(peaks, matched)
        deleteat!(proms, matched)
    end

    return peaks, proms
end

function peakproms!(pks::NamedTuple{names,Ts}; strict=true, min=nothing, max=nothing) where {names,Ts}
    if !hasproperty(pks, :proms)
        # Wait to filter until after merging `pks`
        _, proms = peakproms(pks.indices, pks.data; strict)
        pks = NamedTuple{(names..., :proms)}((values(pks)..., proms))
    end
    filterpeaks!(pks, :proms; min, max)
    return pks
end

"""
    peakproms!(; [strict, min, max]) -> Function

Create a function, `f(pks::NamedTuple)`, that calculates and filters (mutates) the peak
prominences and other fields of its argument, `pks`, using any given keyword arguments.

# Examples
```jldoctest
julia> findmaxima([0,5,2,3,3,1,4,0]) |> peakproms!(; min=2)
(indices = [2, 7], heights = [5, 4], data = [0, 5, 2, 3, 3, 1, 4, 0], proms = Union{Missing, Int64}[5, 3])
```
"""
peakproms!(; kwargs...) = function _curried_peakproms!(pks)
    return peakproms!(pks; kwargs...)
end

