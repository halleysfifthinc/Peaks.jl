"""
    peakheights(peaks, heights;
        minheight=nothing,
        maxheight=nothing
    ) -> (peaks, heights)

Return a copy of `peaks` and `heights` where peak heights are removed if less than
`minheight` and/or greater than `maxheight`.

See also: [`peakproms`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> xpks, vals = findmaxima(x)
([2, 4, 7], [5, 3, 4])

julia> peakheights(xpks, vals; maxheight=4)
([4, 7], [3, 4])

julia> peakheights(xpks, vals; minheight=4.5)
([2], [5])
```
"""
function peakheights(
    peaks::AbstractVector{Int}, heights::AbstractVector;
    minheight=nothing, maxheight=nothing
)
    peakheights!(copy(peaks), copy(heights); minheight=minheight, maxheight=maxheight)
end

"""
    peakheights!(peaks, heights;
        minheight=nothing,
        maxheight=nothing
    ) -> (peaks, heights)

Modify and return `peaks` and `heights` by removing peaks that are less than `minheight` or greater
than `maxheight`.

See also: [`peakproms`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> xpks, vals = findmaxima(x)
([2, 4, 7], [5, 3, 4])

julia> peakheights!(xpks, vals; maxheight=4);

julia> xpks, vals
([4, 7], [3, 4])
```
"""
function peakheights!(
    peaks::Vector{Int}, heights::AbstractVector{T};
    minheight=nothing, maxheight=nothing
) where {T}
    length(peaks) == length(heights) || throw(DimensionMismatch("length of `peaks`, $(length(peaks)), does not match the length of `heights`, $(length(heights))"))
    if !isnothing(minheight) || !isnothing(maxheight)
        lo = something(minheight, typemin(Base.nonmissingtype(T)))
        up = something(maxheight, typemax(Base.nonmissingtype(T)))
        matched = findall(x -> !(lo ≤ x ≤ up), heights)
        deleteat!(peaks, matched)
        deleteat!(heights, matched)
    end

    return peaks, heights
end


