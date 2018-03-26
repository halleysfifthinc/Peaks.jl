__precompile__()

module Peaks

export PeakAlg

export findpeaks,
       maxima,
       minima

function findpeaks(x::AbstractVector)
    return maxima(x)
end

for (funcname, comp) in ((:maxima, :<),
                         (:minima, :>))

    @eval begin
        function ($funcname)(x::AbstractVector,
                        w::Int=1,
                        includebounds::Bool=false)
            w > 0 || throw(ArgumentError("window cannot be negative"))
            xlen = length(x)
            idxs = Int[]

            # There can't be more than one peak every `w` elements
            maxN = includebounds ? cld(xlen,2*w) : fld(xlen,2*w)
            N = 0

            ww = [ collect(w:-1:1); collect(-w:-1) ]

            sizehint!(idxs,maxN)

            i = 1

            if includebounds
                @inbounds while i <= w
                    peak = true
                    for j in max(i-w, 1):(i+w)
                        if ($comp)(x[i], x[j])
                            peak &= false
                            break # No reason to continue checking if the rest of the elements qualify
                        end
                    end

                    if peak
                        push!(idxs,i)
                        N += 1
                        i += w # There can't be another peak for at least `w` more elements
                    end
                    i += 1
                end
            end

            i = 1+w
            maxI = xlen-w

            @inbounds while i <= maxI
                peak = true
                for j in ww # For all elements within the window
                    if ($comp)(x[i], x[i+j])
                        peak &= false
                        break # No reason to continue checking if the rest of the elements qualify
                    end
                end

                if peak
                    push!(idxs,i)
                    N += 1
                    i += w # There can't be another peak for at least `w` more elements
                end
                i += 1
            end

            if includebounds
                @inbounds while i <= xlen
                    peak = true
                    for j in (i-w):min((i+w),xlen)
                        if ($comp)(x[i], x[j])
                            peak &= false
                            break # No reason to continue checking if the rest of the elements qualify
                        end
                    end

                    if peak
                        push!(idxs,i)
                        N += 1
                        i += w # There can't be another peak for at least `w` more elements
                    end
                    i += 1
                end
            end

            resize!(idxs,N)
            sizehint!(idxs,N)

            return idxs
        end
    end
end

end # module
