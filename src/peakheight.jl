"""
    peakheights(indices, heights; [min, max]) -> (indices, heights)
    peakheights(pks::NamedTuple; [min, max]) -> NamedTuple


Return a copy of `indices` and `heights` where peaks are removed if their height is less than
`min` and/or greater than `max`.

If a NamedTuple `pks` is given, a new NamedTuple is returned with filtered copies of the
fields in `pks`. `pks` must have `:indices` and `:heights` fields, at a minimum.
The following fields will also be copied and filtered if they exist: `:proms`, `:widths`,
and `:edges`.

See also: [`peakprom`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref),
[`filterpeaks`](@ref)

# Examples

## NamedTuple API:
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> nt = findmaxima(x)
(indices = [2, 4, 7], heights = [5, 3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> peakheights(nt; max=4)
(indices = [4, 7], heights = [3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])
```

## Seperate vector API:
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> indices, heights = findmaxima(x)
(indices = [2, 4, 7], heights = [5, 3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> peakheights(indices, heights; max=4)
([4, 7], [3, 4])
```
"""
function peakheights(
    indices::AbstractVector{Int}, heights::AbstractVector;
    minheight=nothing, maxheight=nothing,
    min=minheight, max=maxheight
)
    if !isnothing(minheight)
        Base.depwarn("Keyword `minheight` has been renamed to `min`", :peakheights!)
    end
    if !isnothing(maxheight)
        Base.depwarn("Keyword `maxheight` has been renamed to `max`", :peakheights!)
    end
    peakheights!(copy(indices), copy(heights); min=min, max=max)
end

peakheights(pks::NamedTuple; kwargs...) = peakheights!(deepcopy(pks); kwargs...)

"""
    peakheights(; min=nothing, max=nothing) -> Function

Return a function that acts just like the NamedTuple method for 
`peakheights`, but with a not-yet-specified first argument (named tuple). 
This allows a convenient workflow for chaining operations.

## Examples
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> nt = findmaxima(x) |> peakheights(max=4)
(indices = [4, 7], heights = [3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])
```
"""
peakheights(; kwargs...) = function curried_peakheights(pks)
    return peakheights(deepcopy(pks); kwargs...)
end

"""
    peakheights!(indices, heights; [min, max]) -> (indices, heights)
    peakheights!(pks::NamedTuple; [min, max]) -> NamedTuple

Filter (mutate) and return `indices` and `heights` by removing peaks that are less than `min`
or greater than `max`.

If a NamedTuple `pks` is given, a new NamedTuple is returned with the same fields as in
`pks`. `pks` must have `:indices` and `:heights` fields, at a minimum. The
following fields will also be filtered (mutated), if they exist: `:proms`, `:widths`, and
`:edges`.

See also: [`peakprom`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref),
[`filterpeaks!`](@ref)

# Examples
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> indices, heights = findmaxima(x)
(indices = [2, 4, 7], heights = [5, 3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> peakheights!(indices, heights; max=4)
(indices = [2], heights = [5], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> peakheights!((;indices, heights, data=x); max=4)
(indices = [2], heights = [5], data = [0, 5, 2, 3, 3, 1, 4, 0])
```
"""
function peakheights!(
    peaks::Vector{Int}, heights::AbstractVector{T};
    minheight=nothing, maxheight=nothing,
    min=minheight, max=maxheight
) where {T}
    if !isnothing(minheight)
        Base.depwarn("Keyword `minheight` has been renamed to `min`", :peakheights!)
    end
    if !isnothing(maxheight)
        Base.depwarn("Keyword `maxheight` has been renamed to `max`", :peakheights!)
    end
    length(peaks) == length(heights) || throw(DimensionMismatch("length of `peaks`, $(length(peaks)), does not match the length of `heights`, $(length(heights))"))
    if !isnothing(min) || !isnothing(max)
        lo = something(min, typemin(Base.nonmissingtype(T)))
        up = something(max, typemax(Base.nonmissingtype(T)))
        matched = findall(x -> !(lo ≤ x ≤ up), heights)
        deleteat!(peaks, matched)
        deleteat!(heights, matched)
    end

    return peaks, heights
end

function peakheights!(pks::NamedTuple; minheight=nothing, maxheight=nothing, min=minheight, max=maxheight)
    if !isnothing(minheight)
        Base.depwarn("Keyword `minheight` has been renamed to `min`", :peakheights!)
    end
    if !isnothing(maxheight)
        Base.depwarn("Keyword `maxheight` has been renamed to `max`", :peakheights!)
    end
    filterpeaks!(pks, min, max, :heights)
    return pks
end

"""
    peakheights!(; min=nothing, max=nothing) -> Function

Create a function that filters (mutates) the fields of a singular NamedTuple argument `pks`.

This allows for piping several `Peaks.jl` functions together.

# Examples
```jldoctest
julia> pks |> peakheights!(; max=4)
(indices = [2], heights = [5], data = [0, 5, 2, 3, 3, 1, 4, 0])
```
"""
peakheights!(; kwargs...) = function curried_peakheights!(pks)
    peakheights!(pks; kwargs...)
end

