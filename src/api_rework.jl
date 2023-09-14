# I am putting everything in here for now. The contents of this file should be moved around in the future.
"""
    findpeaks(x) -> NamedTuple

Find the peaks in a vector `x`. 
A `NamedTuple` is returned with the original vector 
in the field `data`, and the indices of the peaks 
in the field `inds`.

This function serves as the entry-point for other 
functions such as `peakproms!` and `peakwidths!`
"""
function findpeaks(x::AbstractVector)
    data = x
    inds, _ = findmaxima(x)
    return (data=x, inds=inds)
end

"""
    _filter_fields!(pks, new_inds)

Internal function to mutate the vectors in all fields of `pks` 
to remove any element with a corresponding `pks.inds` entry 
not present in `new_inds`.
"""
function _filter_fields!(pks, new_inds)
    # Add checks on `pks` to see if it has fields "data" and "inds"?
    mask = pks.inds .∉ Ref(new_inds)
    for field in filter(!=(:data), propertynames(pks)) # Do not want to mutate data vector
        v_to_be_mutated = getfield(pks, field)
        if length(v_to_be_mutated) != length(mask)
            error("Internal error: The length of the vector in field `$field` is $(length(v_to_be_mutated)), but was expected to be $(length(mask))")
        end
        deleteat!(getfield(pks, field), mask)
    end
    return nothing
end

"""
    peakproms!(pks::NamedTuple; min=0, max=Inf, strict=true)
    peakproms!(; min=0, max=Inf, strict=true)

Find the prominences of the peaks in `pks`, and filter out any peak 
with a prominence smaller than `min` or greater than `max`.
The prominences are returned in the field `:proms` of the returned named tuple.

If the positional argument `pks` is omitted, an anonymus function is returned 
that performs the action (adding field `:proms` and filtering) to its input.
"""
function peakproms!(pks::NamedTuple; min=0, max=Inf, strict=true)
    if !hasproperty(pks, :proms)
        # Avoid filtering by min/max/strict here, so that it always happens outside if-statement.
        # Pro: one less edge case. Con: More internal allocations
        _, proms = _peakproms(pks.inds, pks.data)
        pks = merge(pks, (; proms))
    end
    new_inds, _ = _peakproms(pks.inds, pks.data; minprom=min, maxprom=max, strict)
    _filter_fields!(pks, new_inds)
    return pks
end
peakproms!(; kwargs...) = pks -> peakproms!(pks; kwargs...)


"""
    peakwidths!(pks::NamedTuple; min=0, max=Inf, relheight=0.5, strict=true)
    peakwidths!(; min=0, max=Inf, relheight=0.5, strict=true)

Find the widths of the peaks in `pks`, and filter out any peak 
with a width smaller than `min` or greater than `max`.
The widths are returned in the field `:widths` of the returned named tuple.
The edges of the peaks are also added in the field `:edges`.

If the positional argument `pks` is omitted, an anonymus function is returned 
that performs the action (adding fields `:widths` and `:edges` and filtering) to its input.

Note: If `pks` does not have a field `proms`, it is added. This is 
because it is needed to calculate the peak width.
"""
function peakwidths!(pks::NamedTuple; min=0, max=Inf, relheight=0.5, strict=true)
    if !hasproperty(pks, :proms)  # Add proms if needed
        pks = peakproms!(pks)
    end
    if xor(hasproperty(pks, :widths), hasproperty(pks, :edges))
        throw(ArgumentError("The named tuple `pks` (first argument to `peakwidths!` is expected have both the fields `:widths` and `:edges`, or to have neither of them. The provided `pks` only has one of them."))
    end
    if !hasproperty(pks, :widths)
        # Avoid filtering by min/max/strict here, so that it always happens outside if-statement.
        # Pro: one less edge case. Con: More internal allocations
        _, widths, leftedges, rightedges = _peakwidths(pks.inds, pks.data, pks.proms)
        pks = merge(pks, (; widths, edges=collect(zip(leftedges, rightedges))))
    end
    new_inds, _ = _peakwidths(pks.inds, pks.data, pks.proms; minwidth=min, maxwidth=max, relheight, strict)
    _filter_fields!(pks, new_inds)
    return pks
end
peakwidths!(; kwargs...) = pks -> peakwidths!(pks; kwargs...)

"""
    peakheights!(pks::NamedTuple; min=0, max=Inf)
    peakheights!(; min=0, max=Inf)

Find the heights of the peaks in `pks`, and filter out any peak 
with a heights smaller than `min` or greater than `max`.
The heights are returned in the field `:heights` of the returned named tuple.

If the positional argument `pks` is omitted, an anonymus function is returned 
that performs the action (adding field `heights` and filtering) to its input.
"""
function peakheights!(pks::NamedTuple; min=0, max=Inf)
    if !hasproperty(pks, :heights)
        # Avoid filtering by min/max here, so that it always happens outside if-statement.
        # Pro: one less edge case. Con: More internal allocations
        heights = pks.data[pks.inds]
        pks = merge(pks, (; heights))
    end
    new_inds, _ = _peakheights(pks.inds, pks.heights; minheight=min, maxheight=max)
    _filter_fields!(pks, new_inds)
    return pks
end
peakheights!(; kwargs...) = pks -> peakheights!(pks; kwargs...)