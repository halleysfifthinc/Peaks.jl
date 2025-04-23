"""
    ismaxima(i, x[, w=1; strict=true]) -> Bool

Test if `i` is a maxima in `x`, where the maxima `i` is either the maximum of `x[i-w:i+w]`
or the first index of a plateau.

A plateau is defined as a maxima with consecutive equal (`==`) maximal values which
are bounded by lesser values immediately before and after the consecutive maximal values.

See also: [`findnextmaxima`](@ref)
"""
ismaxima(i, x, w=1; strict=true)::Bool = findnextextrema(<, x, i, w, strict) === i

"""
    isminima(i, x[, w=1; strict=true]) -> Bool

Test if `i` is a minima in `x`, where the minima `i` is either the minimum of `x[i-w:i+w]`
or the first index of a plateau.

A plateau is defined as a minima with consecutive equal (`==`) minimal values which
are bounded by greater values immediately before and after the consecutive minimal values.

See also: [`findnextminima`](@ref)
"""
isminima(i, x, w=1; strict=true)::Bool = findnextextrema(>, x, i, w, strict) === i

"""
    isplateau(i, x[, w=1; strict=true]) -> Union{Missing,Bool}

Test if `i` is a plateau in `x`, where a plateau is defined as a maxima or minima with
consecutive equal (`==`) extreme values which are bounded by lesser values immediately
before and after the consecutive values. Returns `false` if `i` is the last index in `x`.

See also: [`ismaxima`](@ref), [`isminima`](@ref)
"""
function isplateau(i, x, w=1; strict=true)
    if ismaxima(i, x, w; strict) || isminima(i, x, w; strict)
        if i === lastindex(x)
            # default unstrict assumption that first/last element can be peak means that we
            # should not assume first/last element is (also) a plateau (too much assuming)
            return false
        else
            return x[i] == x[i+1]
        end
    else
        return false
    end
end

const known_fields = (:indices, :proms, :heights, :widths, :edges)

function check_known_fields_equal_length(pks::NamedTuple)
    if !all(==(length(pks.indices))∘length, pks[k]
            for k in keys(pks) if k in known_fields)
        length_pairs = [feature=>length(pks[feature])
            for feature in known_fields if hasproperty(pks, feature)]
        throw(DimensionMismatch(
        "Expected all standard fields of `pks` (exc. `:data`) to be of equal length. Instead,"*
        "found the following lengths: "*join(length_pairs, ", ", " and ")))
    end
    return nothing
end

function check_has_required_fields(pks::NamedTuple)
    !hasproperty(pks, :indices) &&
        throw(ArgumentError("`pks` is missing required field `:indices`"))
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
function filterpeaks!(pks::NamedTuple, mask::Union{BitVector, Vector{Bool}}; mutatemask=false)
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

    if mutatemask
        mask .= .!mask
    else
        mask = .!mask
    end
    for field in known_fields  # Only risk mutating fields added by this package
        hasproperty(pks, field) || continue  # Do nothing if field is not present
        deleteat!(pks[field], mask)
    end
    return pks
end

function filterpeaks!(pks::NamedTuple, feature::Symbol; min=nothing, max=nothing)
    hasproperty(pks, feature) || throw(ArgumentError("`pks` does not have key `$feature`"))
    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(pks.data)))
        hi = something(max, typemax(Base.nonmissingtype(eltype(pks.data))))
        mask = map(x -> !ismissing(x) && (lo ≤ x ≤ hi), pks[feature])
        filterpeaks!(pks, mask; mutatemask=true)
    end
    return pks
end

"""
    filterpeaks!(pred, pks) -> NamedTuple

Apply a predicate function `pred` to NamedTuple slices (the scalar values related to each
peak, e.g. `(indices=5, heights=3, proms=2)`) to and remove a peak if `pred` returns
`false`.

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

# plotting helper functions
"Get `x[idx::Float]` using linear interpolation."
function interp(x::AbstractVector{<:Real}, idx::Real)
    isinteger(idx) && return x[Int(idx)]
    prev = Int(floor(idx))
    next = Int(ceil(idx))
    return x[prev] + ((idx - prev) * (x[next] - x[prev])) / (next - prev)
end

function interp(x::AbstractVector{<:Real}, idx::AbstractVector{<:Real})
    return interp.(Ref(x), idx)
end

function drop_irrelevant_side(i, peak, y, maxima)
    if maxima
        return ifelse(isminima(round(Int, i), y; strict=false), i, peak)
    else
        return ifelse(ismaxima(round(Int, i), y; strict=false), i, peak)
    end
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
