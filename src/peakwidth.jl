function peakwidth(
    peaks::AbstractVector{Int}, x::AbstractVector, proms::AbstractVector;
    strictbounds=true, relheight=0.5, minwidth=nothing, maxwidth=nothing,
)
    peakwidth!(copy(peaks), x, proms; strictbounds, relheight, minwidth, maxwidth)
end

function peakwidth!(
    peaks::AbstractVector{Int}, x::AbstractVector{T}, proms::AbstractVector;
    strictbounds=true, relheight=0.5, minwidth=nothing, maxwidth=nothing,
) where T
    widths = similar(proms, Union{Missing, T})

    fp = first(peaks)
    if fp > 1 && ((x[fp] < x[fp-1]) === true)
        pktype = :minima
    else
        pktype = :maxima
    end
    cmp = pktype === :maxima ? (≤) : (≥)
    op = pktype === :maxima ? (-) : (+)

    if strictbounds
        lst, fst = missing, missing
    else
        lst = lastindex(x)+1
        fst = firstindex(x)-1
    end

    for i in eachindex(peaks, widths)
        ht = op(x[peaks[i]], relheight*proms[i])
        up = findnext(cmp(ht), x, peaks[i])
        lo = findprev(cmp(ht), x, peaks[i])

        !isnothing(lo) && (lo += (ht - x[lo])/(x[lo+1] - x[lo]))
        !isnothing(up) && (up -= (ht - x[up])/(x[up-1] - x[up]))

        widths[i] = something(up, lst) - something(lo, fst)
    end

    if !isnothing(minwidth) || !isnothing(maxwidth)
        lo = something(minwidth, zero(eltype(widths)))
        up = something(maxwidth, typemax(nonmissingtype(eltype(widths))))
        matched = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), widths)
        deleteat!(peaks, matched)
        deleteat!(widths, matched)
    end

    return peaks, widths
end

