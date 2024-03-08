const known_fields = (:indices, :proms, :heights, :widths, :edges)

function check_known_fields_equal_length(pks::NamedTuple)
    features_to_filter = known_fields

    feature_lengths = [length(pks[feature])
        for feature in features_to_filter if hasproperty(pks, feature)]

    # We refrain from using `allequal` to support Julia < 1.8
    if !all(first(feature_lengths) == feature_lengths[i]
            for i in eachindex(feature_lengths))
        length_pairs = [feature=>length(pks[feature])
            for feature in features_to_filter if hasproperty(pks, feature)]
        throw(DimensionMismatch("Expected all known fields of `pks` to be of equal length. Instead found the following pairs of known field and length: $length_pairs"))
    end
    return nothing
end

function check_has_required_fields(pks::NamedTuple)
    !haskey(pks, :indices) && throw(ArgumentError(
        "`pks` is missing required field `:indices`"))
    return nothing
end

"""
    filterpeaks!(pks::NT, feature; [min, max]) where {NT<:NamedTuple} -> pks::NT
    filterpeaks!(pks::NT, mask) -> pks::NT

Filter the standard `pks` fields where peaks are removed if `pks.\$feature` is less than
`min` and/or greater than `max`. If a `mask` is given, a given peak `i` is filtered (removed) if `mask[i]` is `false`.

Standard Peaks.jl fields of `pks` are `:indices`, `:proms`, `:heights`, `:widths`, `:edges`.

# Examples
```jldoctest
julia> pks = findmaxima([0,5,2,3,3,1,4,0])
(indices = [2, 4, 7], heights = [5, 3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> filterpeaks!(pks, :heights; max=4)
(indices = [4, 7], heights = [3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> pks = findmaxima([0,5,2,3,3,1,4,0]) |> peakproms!();

julia> mask = [pks.heights[i] < 5 &&  pks.proms[i] > 2 for i in eachindex(pks.indices)]
3-element Vector{Bool}:
 0
 0
 1

julia> filterpeaks!(pks, mask)
(indices = [7], heights = [4], data = [0, 5, 2, 3, 3, 1, 4, 0], proms = Union{Missing, Int64}[3])
```
"""
function filterpeaks!(pks::NamedTuple, mask::Union{BitVector, Vector{Bool}})
    # Check lengths first to avoid a dimension mismatch
    # after having filtered some features.
    # feature_mask = hasproperty.(pks, features_to_filter)
    check_has_required_fields(pks)
    check_known_fields_equal_length(pks)

    if length(pks.indices) != length(mask)
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

function filterpeaks!(pks::NamedTuple, feature::Symbol; min=nothing, max=nothing)
    haskey(pks, feature) || throw(ArgumentError("`pks` does not have key `$feature`"))
    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(pks.data)))
        hi = something(max, typemax(Base.nonmissingtype(eltype(pks.data))))
        mask = map(x -> !ismissing(x) && (lo ≤ x ≤ hi), pks[feature])
        filterpeaks!(pks, mask)
    end
    return pks
end

"""
    filterpeaks!(pred, pks) -> NamedTuple

Apply a predicate function `pred` to NamedTuple slices (the scalar values related to each
peak, e.g. `(indices=5, heights=3, proms=2)`)to get a filter-mask . A peak is removed if
`pred` returns `false`.

# Examples
```jldoctest
julia> pks = findmaxima([0,5,2,3,3,1,4,0])
(indices = [2, 4, 7], heights = [5, 3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> filterpeaks!(pks) do nt
           return nt.heights ≥ 5 || nt.heights ≤ 3
       end
(indices = [2, 4], heights = [5, 3], data = [0, 5, 2, 3, 3, 1, 4, 0])
```
"""
function filterpeaks!(pred::Function, pks::NamedTuple)
    check_has_required_fields(pks)
    check_known_fields_equal_length(pks)

    # `pks` key's except for `data`, which isn't filtered
    pks_keys = filter(!(==(:data)), keys(pks))
    mask = map(eachindex(pks.indices)) do i
        nt_slice = NamedTuple{pks_keys}(ntuple(j -> getindex(pks[pks_keys[j]], i), length(pks_keys)))
        return !pred(nt_slice)
    end

    for field in known_fields  # Only risk mutating fields added by this package
        hasproperty(pks, field) || continue  # Do nothing if field is not present
        deleteat!(pks[field], mask)
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
