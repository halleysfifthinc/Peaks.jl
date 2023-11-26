"""
    peakheights!(peaks, heights;
        minheight=nothing,
        maxheight=nothing
    ) -> (peaks, heights)

Modify and return `peaks` and `heights` by removing peaks that are less than `minheight` or greater
than `maxheight`.

See also: [`peakprom`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> xpks, vals = findmaxima(x)
(indices = [2, 4, 7], heights = [5, 3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> peakheights!(xpks, vals; maxheight=4);

julia> xpks, vals
([4, 7], [3, 4])
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


"""
    peakheights(peaks, heights;
        minheight=nothing,
        maxheight=nothing
    ) -> (peaks, heights)

Return a copy of `peaks` and `heights` where peak heights are removed if less than
`minheight` and/or greater than `maxheight`.

See also: [`peakprom`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> xpks, vals = findmaxima(x)
(indices = [2, 4, 7], heights = [5, 3, 4], data = [0, 5, 2, 3, 3, 1, 4, 0])

julia> peakheights(xpks, vals; maxheight=4)
([4, 7], [3, 4])

julia> peakheights(xpks, vals; minheight=4.5)
([2], [5])
```
"""
function peakheights(
    peaks::AbstractVector{Int}, heights::AbstractVector;
    minheight=nothing, maxheight=nothing, 
    min=minheight, max=maxheight
)
    if !isnothing(minheight)
        Base.depwarn("Keyword `minheight` has been renamed to `min`", :peakheights!)
    end
    if !isnothing(maxheight)
        Base.depwarn("Keyword `maxheight` has been renamed to `max`", :peakheights!)
    end
    peakheights!(copy(peaks), copy(heights); min=min, max=max)
end

##!===============================================================================================!##
##!==========================================  New API  ==========================================!##
##!===============================================================================================!##

"""
    peakheights!(pks) -> NamedTuple
    peakheights!() -> Function

# Optional keyword arguments
- `min`: Filter out any peak with a height smaller than `min`.
- `max`: Filter out any peak with a height greater than `min`.

Find the heights of the peaks in `pks`, and filter out any peak 
with a heights smaller than `min` or greater than `max`.
Note that because the peaks returned by `findpeaks` already have 
the feature `heights` calculated, this function is mainly useful to 
filter peaks by a minimum and/or maximum height.

If the positional argument `pks` is omitted, a function is returned such
that `peakheights!(pks)` is equivalent to `pks |> peakheights!`.

Note: This function mutates the vectors stored in the NamedTuple `pks`, 
and not the named tuple itself.

See also: [`peakproms!`](@ref), [`peakwidths!`](@ref)

# Examples
```jldoctest
julia> data = [1, 5, 1, 3, 2];

julia> pks = findmaxima(data)
(indices = [2, 4], heights = [5, 3], data = [1, 5, 1, 3, 2])

julia> pks = peakheights!(pks, min=4)
(indices = [2], heights = [5], data = [1, 5, 1, 3, 2])

julia> data |> findmaxima |> peakheights!(min=4)
(indices = [2], heights = [5], data = [1, 5, 1, 3, 2])
```
"""
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
peakheights!(; kwargs...) = pks -> peakheights!(pks; kwargs...)

"""
    peakheights(pks) -> NamedTuple
    peakheights() -> Function

# Optional keyword arguments
- `min`: Filter out any peak with a height smaller than `min`.
- `max`: Filter out any peak with a height greater than `min`.

Non-mutation version of `peakheights!`. Note that 
this copies all vectors in `pks`, including the data. 
This means that it is less performant. See the docstring for 
`peakheights!` for more information.
"""
peakheights(pks::NamedTuple; kwargs...) = peakheights!(deepcopy(pks); kwargs...)