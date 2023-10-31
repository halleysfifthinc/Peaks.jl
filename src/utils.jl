"""
    filterpeaks!(pks, mask) -> pks

Given a NamedTuple `pks` (as returned by `findpeaks`), filter 
the peaks according to the given `mask` (A `BitVector` or `Vector{Bool}`.
This means if element `i` of `mask` is false, then the peak 
at index `i` will be removed from `pks`.

The functions `peakheights`, `peakproms` and `peakwidths` already 
allow filtering by maximal and minimal values for different peak features.
This function can be used to perform more complicated filtering, such 
as keeping a peak if it has a certain height _or_ a certain width.

# Examples
ToDo: Make examples
"""
function filterpeaks!(pks::NamedTuple, mask::Union{BitVector, Vector{Bool}})
    features_to_filter = (:indices, :proms, :heights, :widths, :edges)

    # Check lengths first to avoid a dimension mismatch 
    # after having filtered some features.
    for field in features_to_filter
        hasproperty(pks, field) || continue  # Do nothing if field is not present
        if length(mask) != length(getfield(pks, field))
            throw(DimensionMismatch(
            "Length of `mask` is ($(length(mask))), but the length 
            of `pks.$field` is $(length(getfield(pks, field))). 
            This means that the given mask can not be used to filter 
            the field `$field`."))
        end
    end

    for field in features_to_filter  # Only risk mutating fields added by this package
        hasproperty(pks, field) || continue  # Do nothing if field is not present
        v_to_be_mutated = getfield(pks, field)
        deleteat!(v_to_be_mutated, mask)
    end
    return nothing
end
export filterpeaks!



#====================================================================
We store a version of findpeak here, as it might be implemented soon, 
and I did not want to throw away this implementation
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
=#