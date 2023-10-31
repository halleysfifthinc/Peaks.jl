# I am putting everything in here for now. The contents of this file should be moved around in the future.

"""
    peakproms!(pks) --> NamedTuple
    peakproms!() --> Function

# Optional keyword arguments
- `min`: Filter out any peak with a height smaller than `min`
- `max`: Filter out any peak with a height greater than `min`
- `relheight`: BLABLA. Defaults to `0.5`
- `strict`: BLABLA. Default to `true`

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
        _, proms = _peakproms(pks.indices, pks.data; strict)
        pks = merge(pks, (; proms))
    end
    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(pks.data)))
        up = something(max, typemax(Base.nonmissingtype(eltype(pks.data))))
        mask = map(x -> !ismissing(x) && !(lo ≤ x ≤ up), pks.proms)
        filterpeaks!(pks, mask)
    end
    return pks
end
peakproms!(; kwargs...) = pks -> peakproms!(pks; kwargs...)

"""
    peakwidths!(pks) --> NamedTuple
    peakwidths!() --> Function

# Optional keyword arguments
- `min`: Filter out any peak with a height smaller than `min`
- `max`: Filter out any peak with a height greater than `min`
- `relheight`: BLABLA. Defaults to `0.5`
- `strict`: BLABLA. Default to `true`

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
        pks = peakproms!(pks; strict)
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
        lo = something(min, zero(eltype(pks.data)))
        up = something(max, typemax(Base.nonmissingtype(eltype(pks.data))))
        mask = map(x -> !ismissing(x) && !(lo ≤ x ≤ up), pks.widths)
        filterpeaks!(pks, mask)
    end
    return pks
end
peakwidths!(; kwargs...) = pks -> peakwidths!(pks; kwargs...)

"""
    peakheights!(pks) --> NamedTuple
    peakheights!() --> Function

# Optional keyword arguments
- `min`: Filter out any peak with a height smaller than `min`
- `max`: Filter out any peak with a height greater than `min`

Find the heights of the peaks in `pks`, and filter out any peak 
with a heights smaller than `min` or greater than `max`.
Note that because the peaks returned by `findpeaks` already have 
the feature `heights` calculated, this function is mainly useful to 
filter peaks by a minimum and/or maximum height.

If the positional argument `pks` is omitted, an anonymus function is returned that performs the action (adding field `heights` and filtering) to its input.
"""
function peakheights!(pks::NamedTuple; minheight=nothing, maxheight=nothing, min=minheight, max=maxheight)
    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(pks.data)))
        up = something(max, typemax(Base.nonmissingtype(eltype(pks.data))))
        mask = map(x -> !ismissing(x) && !(lo ≤ x ≤ up), pks.heights)
        filterpeaks!(pks, mask)
    end
    return pks
end
peakheights!(; kwargs...) = pks -> peakheights!(pks; kwargs...)