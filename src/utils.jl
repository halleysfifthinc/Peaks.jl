const known_fields = (:indices, :proms, :heights, :widths, :edges)

function check_known_fields_equal_length(pks::NamedTuple)
    features_to_filter = known_fields

    feature_lengths = [length(pks[feature]) for feature in features_to_filter if hasproperty(pks, feature)]

    # We refrain from using `allequal` to support Julia < 1.8
    if !all(first(feature_lengths) == feature_lengths[i] for i in eachindex(feature_lengths))
        length_pairs = [feature=>length(pks[feature]) for feature in features_to_filter if hasproperty(pks, feature)]
        throw(DimensionMismatch("Expected all known fields of `pks` to be of equal length. Instead found the following pairs of known field and length:\n$length_pairs
        This should not happen, and indicates that the argument `pks` been created or modified by something outside Peaks.jl"))
    end
    return nothing
end

function check_has_known_field(pks::NamedTuple)
    if !any(hasproperty(pks, prop) for prop in known_fields)
        throw(ArgumentError(
            "Attempting to filter a named tuple `pks` that contains none of the known fields $known_fields. Because 
            this is thought to be an error, this error is thrown to help catch the original error."
        ))
    end
    return nothing
end

"""
    filterpeaks!(pks, mask) -> pks

Filter the known fields of `pks` using `mask`, by removing elements of vectors 
in `pks` fields for which `mask[i]` is `false`. Known fields 
of `pks` are `:indices`, `:proms`, `:heights`, `:widths`, `:edges`.

The functions `peakheights`, `peakproms` and `peakwidths` already 
allow filtering by maximal and minimal values for different peak features.
This function can be used to perform more complicated filtering, such 
as keeping a peak if it has a certain height _or_ a certain width.

If you find it inconvenient to define the the mask, see also the 
version of `filterpeaks!` that takes a function as its first argument.

# Examples
julia> data = [1, 2, 3, 2, 3, 4, 0];

julia> pks = findmaxima(data)
(indices = [3, 6], heights = [3, 4], data = [1, 2, 3, 2, 3, 4, 0])

julia> pks = peakwidths!(pks)
(indices = [3, 6], heights = [3, 4], data = [1, 2, 3, 2, 3, 4, 0], proms = Union{Missing, Int64}[1, 3], widths = [1.0, 1.875], edges = [(2.5, 3.5), (4.5, 6.375)])

julia> # We can demand that the peak height is greater than 3.5 AND that the width is smaller than 1.5 
julia> # with the following mask. Note that with this data, that leaves no peaks.

julia> my_mask = [pks.heights[i] > 3.5  &&  pks.widths[i] < 1.5 for i in eachindex(pks.indices)]
2-element Vector{Bool}:
 0
 0

julia> filterpeaks!(pks, my_mask)
(indices = Int64[], heights = Int64[], data = [1, 2, 3, 2, 3, 4, 0], proms = Union{Missing, Int64}[], widths = Float64[], edges = Tuple{Float64, Float64}[])

"""
function filterpeaks!(pks::NamedTuple, mask::Union{BitVector, Vector{Bool}})


    # Check lengths first to avoid a dimension mismatch 
    # after having filtered some features.
    # feature_mask = hasproperty.(pks, features_to_filter)
    check_known_fields_equal_length(pks)
    check_has_known_field(pks)

    # At this point we know that all feature_length are equal, and do not need to check it again
    # pks[1] returns the indices.
    if length(pks[1]) != length(mask)
        throw(DimensionMismatch(
        "Length of `mask` is $(length(mask)), but the length of each of the known fields of `pks` is $(length(pks[1])). 
        This means that the given mask can not be used to filter the given named tuple `pks`."
        ))
    end

    for field in known_fields  # Only risk mutating fields added by this package
        hasproperty(pks, field) || continue  # Do nothing if field is not present
        deleteat!(pks[field], .!mask)
    end
    return pks
end

# This method gets no docstring, as it is intended for internal use.
"""
    filterpeaks!(pks, feature; min=nothing, max=nothing)

Calculate `mask` as `[min < x < max for x in pks.feature]`, 
and call `filterpeaks!(pks, mask)`.

This method is intended for internal use. If you want to filter 
by a minimal or maximal value of some peak feature, concider using 
the keyword arguments `min` and `max` in the function that calculated 
that feature, such as `peakproms!(pks, min=2)`.
"""
function filterpeaks!(pks::NamedTuple, feature::Symbol; min=nothing, max=nothing)
    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(pks.data)))
        up = something(max, typemax(Base.nonmissingtype(eltype(pks.data))))
        mask = map(x -> !ismissing(x) && (lo ≤ x ≤ up), pks[feature])
        filterpeaks!(pks, mask)
    end
    return nothing
end

"""
    filterpeaks!(pred, pks) -> NamedTuple

Apply a predicate function `pred` to named tuple slices to get a filter-mask.
An example of a "named tuple slice" is `(indices=5, heights=3, proms=2)`.
If `pred` returns `false` for peak number `i`, element `i` is removed from 
the known fields of `pks`. 

# Examples
julia> data = [1, 2, 3, 2, 3, 4, 0];

julia> pks = findmaxima(data)
(indices = [3, 6], heights = [3, 4], data = [1, 2, 3, 2, 3, 4, 0])

julia> pks = peakwidths!(pks);

julia> data = [1, 2, 3, 2, 3, 4, 0];

julia> pks = findmaxima(data)
(indices = [3, 6], heights = [3, 4], data = [1, 2, 3, 2, 3, 4, 0])

julia> pks = peakwidths!(pks)
(indices = [3, 6], heights = [3, 4], data = [1, 2, 3, 2, 3, 4, 0], proms = Union{Missing, Int64}[1, 3], widths = [1.0, 1.875], edges = [(2.5, 3.5), (4.5, 6.375)])

julia> filterpeaks!(pks) do nt_slice
           nt_slice.heights > 3.5  &&  nt_slice.widths > 1.8
       end
(indices = [6], heights = [4], data = [1, 2, 3, 2, 3, 4, 0], proms = Union{Missing, Int64}[3], widths = [1.875], edges = [(4.5, 6.375)])

"""
function filterpeaks!(pred::Function, pks::NamedTuple)
    check_has_known_field(pks)
    check_known_fields_equal_length(pks)

    mask = map(eachindex(pks[1])) do i
        # :data is included in the nt_slice, but that should not be a problem
        itr = zip(keys(pks), getindex.(values(pks), i))
        nt_slice = NamedTuple(itr)
        return pred(nt_slice)
    end

    for field in known_fields  # Only risk mutating fields added by this package
        hasproperty(pks, field) || continue  # Do nothing if field is not present
        deleteat!(pks[field], .!mask)
    end
    return pks
end


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