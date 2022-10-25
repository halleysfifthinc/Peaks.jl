"""
    peakheights(peaks, x;
        minheight=nothing,
        maxheight=nothing
    ) -> (peaks, heights)

Return the `peaks` and peak heights in `x` (e.g. `x[peaks]`) for peaks between `minheight`
and `maxheight`.

See also: [`peakprom`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = [0,5,2,3,3,1,4,0];

julia> xpks = argmaxima(x)
3-element Vector{Int64}:
 2
 4
 7

julia> peakheights(xpks, x)
([2, 4, 7], [5, 3, 4])

julia> peakheights(xpks, x; maxheight=4)
([4, 7], [3, 4])

julia> peakheights(xpks, x; minheight=4.5)
([2], [5])
```
"""
function peakheights(
    peaks::AbstractVector{Int}, x::AbstractVector;
    minheight=nothing, maxheight=nothing,
)
    if !isnothing(minheight) || !isnothing(maxheight)
        _peaks = copy(peaks)
    else
        # peaks will not be modified
        _peaks = peaks
    end
    peakheights!(_peaks, x; minheight=minheight, maxheight=maxheight)
end

"""
    peakheights(peaks, x;
        minheight=nothing,
        maxheight=nothing
    ) -> (peaks, heights)

Modify `peaks` by filtering peaks that are not between `minheight` and `maxheight`. Return
the modified `peaks` and the peak heights (e.g. `x[peaks]`).

See also: [`peakprom`](@ref), [`peakwidths`](@ref), [`findmaxima`](@ref)
"""
function peakheights!(
    peaks::Vector{Int}, x::AbstractVector{T};
    minheight=nothing, maxheight=nothing
) where T
    if !isnothing(minheight) || !isnothing(maxheight)
        lo = something(minheight, typemin(Base.nonmissingtype(T)))
        up = something(maxheight, typemax(Base.nonmissingtype(T)))
        matched = findall(x -> !(lo ≤ x ≤ up), @view x[peaks])
        deleteat!(peaks, matched)
    end

    return peaks, x[peaks]
end


