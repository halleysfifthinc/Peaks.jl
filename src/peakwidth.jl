"""
    peakwidths(peaks, x, proms;
        strict=true,
        relheight=0.5,
        minwidth=nothing,
        maxwidth=nothing
    ) -> (peaks, widths, leftedge, rightedge)

Calculate the widths of `peaks` in `x` at a reference level based on `proms` and
`relheight`, and removing peaks with widths less than `minwidth` and/or greater than
`maxwidth`. Returns the peaks, widths, and the left and right edges at the reference level.

Peak width is the distance between the signal crossing a reference level before and after
the peak. Signal crossings are linearly interpolated between indices. The reference level is
the difference between the peak height and `relheight` times the peak prominence. Width
cannot be calculated for a `NaN` or `missing` prominence.

The width for a peak with a gap in the signal (e.g. `NaN`, `missing`) at the reference level
will match the value/type of the signal gap if `strict == true`. For `strict ==
false`, the signal crossing will be linearly interpolated between the edges of the gap.

See also: [`peakprom`](@ref), [`findminima`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = [0,1,0,-1.];

julia> xpks = argmaxima(x)
1-element Vector{Int64}:
 2

julia> peakwidths(xpks, x, [1])
([2], [1.0], [1.5], [2.5])

julia> x[3] = NaN;

julia> peakwidths(xpks, x, [1])
([2], [NaN], [1.5], [NaN])

julia> peakwidths(xpks, x, [1]; strict=false)
([2], [1.0], [1.5], [2.5])
```
"""
function _peakwidths(
    peaks::AbstractVector{Int}, x::AbstractVector, proms::AbstractVector;
    strict=true, relheight=0.5, minwidth=nothing, maxwidth=nothing
)
    if !isnothing(minwidth) || !isnothing(maxwidth)
        _peaks = copy(peaks)
    else
        # peaks will not be modified
        _peaks = peaks
    end
    _peakwidths!(_peaks, x, proms; strict=strict, relheight=relheight,
        minwidth=minwidth, maxwidth=maxwidth)
end

"""
    peakwidths!(peaks, x, proms;
        strict=true,
        relheight=0.5,
        minwidth=nothing,
        maxwidth=nothing
    ) -> (peaks, widths, leftedge, rightedge)

Calculate the widths of `peaks` in `x` at a reference level based on `proms` and
`relheight`, removing peaks with widths less than `minwidth` and/or greater than `maxwidth`.
Returns the modified peaks, widths, and the left and right edges at the reference level.

See also: [`peakwidths`](@ref), [`peakproms`](@ref), [`findminima`](@ref), [`findmaxima`](@ref)
"""
function _peakwidths!(
    peaks::AbstractVector{Int}, x::AbstractVector{T}, proms::AbstractVector{U};
    strict=true, relheight=0.5, minwidth=nothing, maxwidth=nothing
) where {T,U}
    if !isnothing(minwidth) && !isnothing(maxwidth)
        minwidth < maxwidth || throw(ArgumentError("maxwidth must be greater than minwidth"))
    end
    all(∈(eachindex(x)), peaks) ||
        throw(ArgumentError("peaks contains invalid indices to x"))

    # if peaks was calculated with strict=false, first(peaks) could be minima at firstindex
    fp = length(peaks) > 1 ? peaks[2] : first(peaks)
    if fp > 1 && ((x[fp] < x[fp-1]) === true)
        pktype = :minima
    else
        pktype = :maxima
    end
    cmp = pktype === :maxima ? (≤) : (≥)
    op = pktype === :maxima ? (-) : (+)

    V1 = promote_type(T, U)
    _bad = Missing <: V1 ? missing : float(Int)(NaN)

    V = promote_type(V1, typeof(_bad))
    ledge = similar(proms, typeof(one(V)/one(V)))  # typeof(one(V)/one(V)) because the 
    redge = similar(proms, typeof(one(V)/one(V)))  # vector eltype need to survive division

    if strict
        lst, fst = _bad, _bad
    else
        lst = lastindex(x)
        fst = firstindex(x)
    end

    for i in eachindex(peaks, ledge, redge)
        prom = proms[i]
        if ismissing(prom) || isnan(prom)
            redge[i] = _bad
            ledge[i] = _bad
        else
            ht = op(x[peaks[i]], relheight * proms[i])
            lo = findprev(v -> !ismissing(v) && cmp(v, ht), x, peaks[i])
            up = findnext(v -> !ismissing(v) && cmp(v, ht), x, peaks[i])

            if !strict
                if !isnothing(lo)
                    lo1 = findnext(v -> !ismissing(v) && cmp(ht, v), x, lo + 1)
                    lo += (ht - x[lo]) / (x[lo1] - x[lo]) * (lo1 - lo)
                end
                if !isnothing(up)
                    up1 = findprev(v -> !ismissing(v) && cmp(ht, v), x, up - 1)
                    up -= (ht - x[up]) / (x[up1] - x[up]) * (up - up1)
                end
            else
                !isnothing(lo) && (lo += (ht - x[lo]) / (x[lo+1] - x[lo]))
                !isnothing(up) && (up -= (ht - x[up]) / (x[up-1] - x[up]))
            end
            redge[i] = something(up, lst)
            ledge[i] = something(lo, fst)
        end
    end

    widths::Vector{V} = redge - ledge

    if !isnothing(minwidth) || !isnothing(maxwidth)
        lo = something(minwidth, zero(eltype(widths)))
        up = something(maxwidth, typemax(Base.nonmissingtype(eltype(widths))))
        matched = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), widths)
        deleteat!(peaks, matched)
        deleteat!(ledge, matched)
        deleteat!(redge, matched)
        deleteat!(widths, matched)
    end

    return peaks, widths, ledge, redge
end

