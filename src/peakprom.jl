using UnsafeArrays

struct Maxima; end
struct Minima; end

for (comps, Extrema) in (((:>, minimum, max), Maxima),
                       ((:<, maximum, min), Minima))

    comp, argcomp1, argcomp2 = comps

    @eval begin
        function peakprom(x::AbstractVector{T}, ::$Extrema, w=1, minprom::T=zero(T)) where T
            if ($Extrema) === Maxima
                m = argmaxima(x, w)
            else
                m = argminima(x, w)
            end

            M = lastindex(m)
            lbegin = firstindex(x)
            lend = lastindex(x)
            proms = similar(x,M)

            @inbounds for i in eachindex(m)
                lb, rb = lbegin, lend
                lcan, rcan = zero(T), zero(T)

                # Find left bound
                for ii in (i-1):-1:1
                    if ($comp)(x[m[ii]], x[m[i]])
                        lb = m[ii]
                        break
                    end
                end

                # Find left bound
                for ii in (i+1):M
                    if ($comp)(x[m[ii]], x[m[i]])
                        rb = m[ii]
                        break
                    end
                end

                lcan = ($argcomp1)(skipmissing(uview(x, lb:(m[i] - 1)))) # Slice from lower bound to the index prior to the current extrema
                rcan = ($argcomp1)(skipmissing(uview(x, (m[i] + 1):rb))) # Slice corollary upper side

                proms[i] = abs(x[m[i]] - ($argcomp2)(lcan, rcan))
            end

            if minprom != zero(T)
                matched = findall(x -> x >= minprom, proms)
                return (m[matched], proms[matched])
            else
                return (m, proms)
            end
        end

    end
end

@doc """
    peakprom(x, ::Maxima[, w=1, minprom=0])

Find the indices of all local maxima and their prominences in `x` matching the conditions `w` and `minprom`.
`w` sets the minimum allowed distance between maxima. `minprom` sets the minimum prominence
(inclusive) of returned maxima.

Peak prominence is calculated as the distance between the current maxima and the highest of
the minimums of the lower and upper bounds. Bounds extend from the next index from the
current maxima to the next maxima higher than the current maxima, or the end of the signal,
which ever comes first.

# Examples
```jldoctest
julia> x = rand(1000);

julia> ma, pa = peakprom(x, Maxima());

julia> mi, pi = peakprom(-x, Minima());

julia> @assert (mi == ma) && (pa == pi)

```

See also: [`argmaxima`](@ref), [`maxima`](@ref)
"""
peakprom(x, ::Maxima)

@doc """
    peakprom(x, ::Minima[, w=1, minprom=0])

Find the indices of all local minima and their prominences in `x` matching the conditions `w` and `minprom`.
`w` sets the minimum allowed distance between minima. `minprom` sets the minimum prominence
(inclusive) of returned minima.

Peak prominence is calculated as the distance between the current minima and the lowest of
the maximums of the lower and upper bounds. Bounds extend from the next index from the
current minima to the next minima lower than the current minima, or the end of the signal,
which ever comes first.

See also: [`argminima`](@ref), [`minima`](@ref)
"""
peakprom(x, ::Minima)

