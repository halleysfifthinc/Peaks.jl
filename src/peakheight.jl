"""
    peakheights(indices, heights; [min, max]) -> (indices, heights)
    peakheights(pks::NamedTuple; [min, max]) -> NamedTuple


Return a copy of `indices` and `heights` where peaks are removed if their height is less than
`min` and/or greater than `max`.

If a NamedTuple `pks` is given, a new NamedTuple is returned with filtered copies of fields
from `pks`. `pks` must have `:indices` and `:heights` fields. The fields `:proms`,
`:widths`, and `:edges` will be filtered if present, and any remaining fields will be
copied unmodified.

See also: [`peakproms`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> pks = findmaxima([0,5,2,3,3,1,4,0])
(indices = [2, 4, 7], heights = [5, 3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> peakheights(pks; max=4)
(indices = [4, 7], heights = [3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> inds, heights = peakheights(pks.indices, pks.heights; max=4)
([4, 7], [3, 4])
```
"""
function peakheights(
    indices::AbstractVector{Int}, heights::AbstractVector;
    min=nothing, max=nothing
)
    peakheights!(copy(indices), copy(heights); min=min, max=max)
end

peakheights(pks::NamedTuple; kwargs...) = peakheights!(deepcopy(pks); kwargs...)

"""
    peakheights(; [min, max]) -> Function

Create a function, `f(pks::NamedTuple)`, that copies and filters the peak heights of its
argument, `pks`, using any given keyword arguments.

# Examples
```jldoctest
julia> findmaxima([0, 5, 2, 3, 3, 1, 4, 0]) |> peakheights(; max=4)
(indices = [4, 7], heights = [3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])
```
"""
peakheights(; kwargs...) = function _curried_peakheights(pks)
    return peakheights(deepcopy(pks); kwargs...)
end

"""
    peakheights!(indices, heights; [min, max]) -> (indices, heights)
    peakheights!(pks::NamedTuple; [min, max]) -> NamedTuple

Filter (mutate) and return `indices` and `heights` by removing peaks that are less than `min`
and/or greater than `max`.

If a NamedTuple `pks` is given, a new NamedTuple is returned with the same fields
(references) from `pks`. `pks` must have `:indices` and `:heights` fields. The fields
`:proms`, `:widths`, and `:edges` will be filtered (mutated) if present, and any remaining
fields will be referenced unmodified.

See also: [`peakproms`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref)
[`filterpeaks!`](@ref)

# Examples
```jldoctest
julia> pks = findmaxima([0,5,2,3,3,1,4,0])
(indices = [2, 4, 7], heights = [5, 3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> peakheights!(pks; max=4)
(indices = [4, 7], heights = [3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> inds, heights = peakheights!(pks.indices, pks.heights; min=3.5)
([7], [4])
```
"""
function peakheights!(
    peaks::Vector{Int}, heights::AbstractVector{T};
    min=nothing, max=nothing
) where {T}
    length(peaks) == length(heights) || throw(DimensionMismatch("length of `peaks`, $(length(peaks)), does not match the length of `heights`, $(length(heights))"))
    if !isnothing(min) || !isnothing(max)
        lo = something(min, typemin(Base.nonmissingtype(T)))
        hi = something(max, typemax(Base.nonmissingtype(T)))
        matched = findall(x -> !(lo ≤ x ≤ hi), heights)
        deleteat!(peaks, matched)
        deleteat!(heights, matched)
    end

    return peaks, heights
end

function peakheights!(pks::NamedTuple; min=nothing, max=nothing)
    filterpeaks!(pks, :heights; min, max)
    return pks
end

"""
    peakheights!(; [min, max]) -> Function

Create a function, `f(pks::NamedTuple)`, that calculates peak heights and then filters
(mutates) the fields of its argument, `pks`, using any given keyword arguments.

# Examples
```jldoctest
julia> findmaxima([0, 5, 2, 3, 3, 1, 4, 0]) |> peakheights!(; max=4)
(indices = [4, 7], heights = [3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])
```
"""
peakheights!(; kwargs...) = function _curried_peakheights!(pks)
    return peakheights!(pks; kwargs...)
end

