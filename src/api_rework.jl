# I am putting everything in here for now. The contents of this file should be moved around in the future.
"""
    findpeaks(x) -> NamedTuple
    findpeaks(x, w=1; strict=true) -> NamedTuple

Find the peaks in a vector `x`, where each maxima i is either 
the maximum of x[i-w:i+w] or the first index of a plateau.
A `NamedTuple` is returned with the original vector 
in the field `data`, and the indices of the peaks 
in the field `indices`.

This function serves as the entry-point for other 
functions such as `peakproms!` and `peakwidths!`
"""
function findpeaks(x::AbstractVector, w::Int=1; strict=true)
    indices, heights = findmaxima(x, w; strict)
    return (data=x, indices=indices, heights=heights)
end

"""
    filterpeaks!(pks, mask) -> Nothing

Given a NamedTuple `pks` (as returned by `findpeaks`), filter 
the peaks according to the given `mask` (A `BitVector` or `Vector{Bool}`. The functions `peakheights`, `peakproms` and `peakwidths` allow 
filtering for a max/min limit of the relevant peak feature.
This function can be used to perform more complicated filtering, such 
as keeping a peak if it has a certain height _or_ a certain width.

# Examples:
ToDo: Make example
"""
function filterpeaks!(pks, mask::BitVector)
    features_to_filter = (:indices, :proms, :heights, :widths, :edges)

    # Check lengths first to avoid a dimension mismatch 
    # after having filtered some features.
    for field in features_to_filter
        if length(mask) != length(getfield(pks, field))
            throw(DimensionMismatch(
            "Length of `mask` is ($(length(mask))), but the length 
            of `pks.$field` is $(length(getfield(pks, field))). 
            This means that the given mask can not be used to filter 
            the field `$field`."))
        end
    end

    for field in features_to_filter  # Only risk mutating fields added by this package
        hasfield(pks, field) || continue  # Do nothing if field is not present
        v_to_be_mutated = getfield(pks, field)
        deleteat!(v_to_be_mutated, mask)
    end
    return nothing
end
export filterpeaks!

"""
    peakproms!(pks::NamedTuple; min=0, max=Inf, strict=true)
    peakproms!(; min=0, max=Inf, strict=true)

Find the prominences of the peaks in `pks`, and filter out any peak 
with a prominence smaller than `min` or greater than `max`.
The prominences are returned in the field `:proms` of the returned named tuple.

If the positional argument `pks` is omitted, an anonymus function is returned 
that performs the action (adding field `:proms` and filtering) to its input.
"""
function peakproms!(pks::NamedTuple; minprom=nothing, maxprom=nothing, min=minprom, max=maxprom, strict=true)
    if !hasproperty(pks, :proms)
        # Avoid filtering by min/max/strict here, so that it always happens outside if-statement.
        # Pro: one less edge case. Con: More internal allocations
        _, proms = _peakproms(pks.indices, pks.data)
        pks = merge(pks, (; proms))
    end
    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(x)))
        up = something(max, typemax(Base.nonmissingtype(eltype(x))))
        mask = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), pks.proms)
        filterpeaks!(pks, mask)
    end
    return pks
end
peakproms!(; kwargs...) = pks -> peakproms!(pks; kwargs...)


"""
    peakwidths!(pks::NamedTuple; min=0, max=Inf, relheight=0.5, strict=true)
    peakwidths!(; min=0, max=Inf, relheight=0.5, strict=true)

Find the widths of the peaks in `pks`, and filter out any peak 
with a width smaller than `min` or greater than `max`.
The widths are returned in the field `:widths` of the returned named tuple.
The edges of the peaks are also added in the field `:edges`.

If the positional argument `pks` is omitted, an anonymus function is returned 
that performs the action (adding fields `:widths` and `:edges` and filtering) to its input.

Note: If `pks` does not have a field `proms`, it is added. This is 
because it is needed to calculate the peak width.
"""
function peakwidths!(pks::NamedTuple; minwidth=nothing, maxwidth=nothing, min=minwidth, max=maxwidth, relheight=0.5, strict=true)
    if !hasproperty(pks, :proms)  # Add proms if needed
        pks = peakproms!(pks)
    end
    if xor(hasproperty(pks, :widths), hasproperty(pks, :edges))
        throw(ArgumentError("The named tuple `pks` (first argument to `peakwidths!` is expected have both the fields `:widths` and `:edges`, or to have neither of them. The provided `pks` only has one of them."))
    end
    if !hasproperty(pks, :widths)
        # Avoid filtering by min/max/strict here, so that it always happens outside if-statement.
        # Pro: one less edge case. Con: More internal allocations
        _, widths, leftedges, rightedges = _peakwidths(pks.indices, pks.data, pks.proms; relheight, strict)
        pks = merge(pks, (; widths, edges=collect(zip(leftedges, rightedges))))
    end
    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(x)))
        up = something(max, typemax(Base.nonmissingtype(eltype(x))))
        mask = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), pks.widths)
        filterpeaks!(pks, mask)
    end
    return pks
end
peakwidths!(; kwargs...) = pks -> peakwidths!(pks; kwargs...)

"""
    peakheights!(pks::NamedTuple; min=0, max=Inf)
    peakheights!(; min=0, max=Inf)

Find the heights of the peaks in `pks`, and filter out any peak 
with a heights smaller than `min` or greater than `max`.
Note that because the peaks returned by `findpeaks` already have 
the feature `heights` calculated, this function is mainly useful to 
filter peaks by a minimum and/or maximum height.

If the positional argument `pks` is omitted, an anonymus function is returned that performs the action (adding field `heights` and filtering) to its input.
"""
function peakheights!(pks::NamedTuple; minheight=nothing, maxheight=nothing, min=minheight, max=maxheight)
    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(x)))
        up = something(max, typemax(Base.nonmissingtype(eltype(x))))
        mask = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), pks.heights)
        filterpeaks!(pks, mask)
    end
    return pks
end
peakheights!(; kwargs...) = pks -> peakheights!(pks; kwargs...)