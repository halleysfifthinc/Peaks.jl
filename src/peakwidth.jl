function peakwidth(
    peaks::AbstractVector{Int}, x::AbstractVector, proms::AbstractVector;
    strictbounds=true, relheight=0.5, minwidth=nothing, maxwidth=nothing,
)
    if !isnothing(minwidth) && !isnothing(maxwidth)
        minwidth < maxwidth || throw(ArgumentError("maxwidth must be greater than minwidth"))
    end
    all(∈(eachindex(x)), peaks) ||
        throw(ArgumentError("peaks contains invalid indices to x"))

    if !isnothing(minwidth) || !isnothing(maxwidth)
        _peaks = copy(peaks)
    else
        # peaks will not be modified
        _peaks = peaks
    end
    peakwidth!(_peaks, x, proms; strictbounds=strictbounds, relheight=relheight,
        minwidth=minwidth, maxwidth=maxwidth)
end

function peakwidth!(
    peaks::AbstractVector{Int}, x::AbstractVector{T}, proms::AbstractVector{U};
    strictbounds=true, relheight=0.5, minwidth=nothing, maxwidth=nothing,
) where {T, U}
    if !isnothing(minwidth) && !isnothing(maxwidth)
        minwidth < maxwidth || throw(ArgumentError("maxwidth must be greater than minwidth"))
    end
    all(∈(eachindex(x)), peaks) ||
        throw(ArgumentError("peaks contains invalid indices to x"))

    fp = first(peaks)
    if fp > 1 && ((x[fp] < x[fp-1]) === true)
        pktype = :minima
    else
        pktype = :maxima
    end
    cmp = pktype === :maxima ? (≤) : (≥)
    op = pktype === :maxima ? (-) : (+)

    V1 = promote_type(T,U)
    _bad = Missing <: V1 ? missing :
           Float64 <: V1 ? NaN :
           Float32 <: V1 ? NaN32 :
           Float16 <: V1 ? NaN16 :
                           missing

    V = promote_type(V, typeof(_bad))
    lower = similar(proms, V)
    upper = similar(proms, V)

    if strictbounds
        lst, fst = _bad, _bad
    else
        lst = lastindex(x)
        fst = firstindex(x)
    end

    for i in eachindex(peaks, lower, upper)
        prom = proms[i]
        if ismissing(prom) || isnan(prom)
            upper[i] = _bad
            lower[i] = _bad
        else
            ht = op(x[peaks[i]], relheight*proms[i])
            lo = findprev(v -> !ismissing(v) && cmp(v,ht), x, peaks[i])
            up = findnext(v -> !ismissing(v) && cmp(v,ht), x, peaks[i])

            if !strictbounds
                if !isnothing(lo)
                    lo1 = findnext(v -> !ismissing(v) && cmp(ht,v), x, lo+1)
                    lo += (ht - x[lo])/(x[lo1] - x[lo])*(lo1 - lo)
                end
                if !isnothing(up)
                    up1 = findprev(v -> !ismissing(v) && cmp(ht,v), x, up-1)
                    up -= (ht - x[up])/(x[up1] - x[up])*(up - up1)
                end
            else
                !isnothing(lo) && (lo += (ht - x[lo])/(x[lo+1] - x[lo]))
                !isnothing(up) && (up -= (ht - x[up])/(x[up-1] - x[up]))
            end

            upper[i] = something(up, lst)
            lower[i] = something(lo, fst)
        end
    end

    widths::Vector{V} = upper - lower

    if !isnothing(minwidth) || !isnothing(maxwidth)
        lo = something(minwidth, zero(eltype(widths)))
        up = something(maxwidth, typemax(Base.nonmissingtype(eltype(widths))))
        matched = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), widths)
        deleteat!(peaks, matched)
        deleteat!(lower, matched)
        deleteat!(upper, matched)
        deleteat!(widths, matched)
    end

    return peaks, widths, lower, upper
end

