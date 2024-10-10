var documenterSearchIndex = {"docs":
[{"location":"#Peaks.jl","page":"Peaks.jl","title":"Peaks.jl","text":"","category":"section"},{"location":"","page":"Peaks.jl","title":"Peaks.jl","text":"argmaxima\nargminima\nfindmaxima\nfindminima\npeakproms\npeakproms!\npeakwidths\npeakwidths!\nfindnextmaxima\nfindnextminima","category":"page"},{"location":"#Peaks.argmaxima","page":"Peaks.jl","title":"Peaks.argmaxima","text":"argmaxima(x[, w=1; strict=true])\n\nFind the indices of the local maxima of x where each maxima i is either the maximum of x[i-w:i+w] or the first index of a plateau.\n\nIf strict is true, all elements of x[i-w:i+w] or the bounds of a plateau must exist and must not be missing or NaN. For strict == false, a maxima is the maximum or first element of all existing and non-NaN or missing elements in x[i-w:i+w] or the bounds of a plateau.\n\nSee also: findmaxima, findnextmaxima\n\nExamples\n\njulia> argmaxima([0,2,0,1,1,0])\n2-element Vector{Int64}:\n 2\n 4\n\njulia> argmaxima([2,0,1,1])\nInt64[]\n\njulia> argmaxima([2,0,1,1]; strict=false)\n2-element Vector{Int64}:\n 1\n 3\n\n\n\n\n\n","category":"function"},{"location":"#Peaks.argminima","page":"Peaks.jl","title":"Peaks.argminima","text":"argminima(x[, w=1; strict=false])\n\nFind the indices of the local minima of x where each minima i is either the minimum of x[i-w:i+w] or the first index of a plateau.\n\nIf strict is true, all elements of x[i-w:i+w] or the bounds of a plateau must exist and must not be missing or NaN. For strict == false, a minima is the minimum or first element of all existing and non-NaN or missing elements in x[i-w:i+w] or the bounds of a plateau.\n\nSee also: findminima, findnextminima\n\nExamples\n\njulia> argminima([3,2,3,1,1,3])\n2-element Vector{Int64}:\n 2\n 4\n\njulia> argminima([2,3,1,1])\nInt64[]\n\njulia> argminima([2,3,1,1]; strict=false)\n2-element Vector{Int64}:\n 1\n 3\n\n\n\n\n\n","category":"function"},{"location":"#Peaks.findmaxima","page":"Peaks.jl","title":"Peaks.findmaxima","text":"findmaxima(x[, w=1; strict=true]) -> (idxs, vals)\n\nFind the indices and values of the local maxima of x where each maxima i is either the maximum of x[i-w:i+w] or the first index of a plateau.\n\nSee also: argmaxima, findnextmaxima\n\n\n\n\n\n","category":"function"},{"location":"#Peaks.findminima","page":"Peaks.jl","title":"Peaks.findminima","text":"findminima(x[, w=1; strict=true]) -> (idxs, vals)\n\nFind the indices and values of the local minima of x where each minima i is either the minimum of x[i-w:i+w] or the first index of a plateau.\n\nSee also: argminima, findnextminima\n\n\n\n\n\n","category":"function"},{"location":"#Peaks.peakproms","page":"Peaks.jl","title":"Peaks.peakproms","text":"peakproms(peaks, x;\n    strict=true,\n    minprom=nothing,\n    maxprom=nothing\n) -> (peaks, proms)\n\nCalculate the prominences of peaks in x, filtering peaks with prominences less than minprom and greater than maxprom, if either are given.\n\nPeak prominence is the absolute height difference between the current peak and the larger of the two adjacent smallest magnitude points between the current peak and adjacent larger peaks or signal ends.\n\nThe prominence for a peak with a NaN or missing between the current peak and either adjacent larger peaks will be NaN or missing if strict == true, or it will be the larger of the smallest non-NaN or missing values between the current peak and adjacent larger peaks for strict == false.\n\nSee also: findminima, findmaxima\n\nExamples\n\njulia> x = [0,5,2,3,3,1,4,0];\n\njulia> xpks = argmaxima(x)\n3-element Vector{Int64}:\n 2\n 4\n 7\n\njulia> peakproms(xpks, x)\n([2, 4, 7], Union{Missing, Int64}[5, 1, 3])\n\njulia> x = [missing,5,2,3,3,1,4,0];\n\njulia> peakproms(xpks, x)\n([2, 4, 7], Union{Missing, Int64}[missing, 1, 3])\n\njulia> peakproms(xpks, x; strict=false)\n([2, 4, 7], Union{Missing, Int64}[5, 1, 3])\n\n\n\n\n\n","category":"function"},{"location":"#Peaks.peakproms!","page":"Peaks.jl","title":"Peaks.peakproms!","text":"peakproms!(peaks, x;\n    strict=true,\n    minprom=nothing,\n    maxprom=nothing\n) -> (peaks, proms)\n\nCalculate the prominences of peaks in x, removing peaks with prominences less than minprom or greater than maxprom, if either are given. Returns the modified peaks and their prominences.\n\nSee also: peakproms, findminima, findmaxima\n\n\n\n\n\n","category":"function"},{"location":"#Peaks.peakwidths","page":"Peaks.jl","title":"Peaks.peakwidths","text":"peakwidths(peaks, x, proms;\n    strict=true,\n    relheight=0.5,\n    minwidth=nothing,\n    maxwidth=nothing\n) -> (peaks, widths, leftedge, rightedge)\n\nCalculate the widths of peaks in x at a reference level based on proms and relheight, removing peaks with widths less than minwidth or greater than maxwidth, if either are given. Returns the modified peaks, widths, and the left and right edges at the reference level.\n\nPeak width is the distance between the signal crossing a reference level before and after the peak. Signal crossings are linearly interpolated between indices. The reference level is the difference between the peak height and relheight times the peak prominence. Width cannot be calculated for a NaN or missing prominence.\n\nThe width for a peak with a gap in the signal (e.g. NaN, missing) at the reference level will match the value/type of the signal gap if strict == true. For strict == false, the signal crossing will be linearly interpolated between the edges of the gap.\n\nSee also: peakprom, findminima, findmaxima\n\nExamples\n\njulia> x = [0,1,0,-1.];\n\njulia> xpks = argmaxima(x)\n1-element Vector{Int64}:\n 2\n\njulia> peakwidths(xpks, x, [1])\n([2], [1.0], [1.5], [2.5])\n\njulia> x[3] = NaN;\n\njulia> peakwidths(xpks, x, [1])\n([2], [NaN], [1.5], [NaN])\n\njulia> peakwidths(xpks, x, [1]; strict=false)\n([2], [1.0], [1.5], [2.5])\n\n\n\n\n\n","category":"function"},{"location":"#Peaks.peakwidths!","page":"Peaks.jl","title":"Peaks.peakwidths!","text":"peakwidths!(peaks, x, proms;\n    strict=true,\n    relheight=0.5,\n    minwidth=nothing,\n    maxwidth=nothing\n) -> (peaks, widths, leftedge, rightedge)\n\nCalculate the widths of peaks in x at a reference level based on proms and relheight, removing peaks with widths less than minwidth or greater than maxwidth, if either are given. Returns the modified peaks, widths, and the left and right edges at the reference level.\n\nSee also: peakwidths, peakproms, findminima, findmaxima\n\n\n\n\n\n","category":"function"},{"location":"#Peaks.findnextmaxima","page":"Peaks.jl","title":"Peaks.findnextmaxima","text":"findnextmaxima(x, i[, w=1; strict=true])\n\nFind the index of the next maxima in x after or including i, where the maxima i is either the maximum of x[i-w:i+w] or the first index of a plateau. Returns lastindex(x) + 1 if no maxima occur after i.\n\nIf strict is true, all elements in x[i-w:i+w] or the bounds of a plateau must exist and must not be missing or NaN. For strict == false, a maxima is the maximum or first element of all existing and non-NaN or missing elements in x[i-w:i+w] or the bounds of a plateau.\n\nSee also: argmaxima\n\nExamples\n\njulia> findnextmaxima([0,2,0,1,1,0], 2)\n2\n\njulia> findnextmaxima([0,2,0,1,1,0], 3)\n4\n\n\n\n\n\n\n","category":"function"},{"location":"#Peaks.findnextminima","page":"Peaks.jl","title":"Peaks.findnextminima","text":"findnextminima(x, i[, w=1, strict=true])\n\nFind the index of the next minima in x after or including i, where the minima i is either the minimum of x[i-w:i+w] or the first index of a plateau. Returns lastindex(x) + 1 if no minima occur after i.\n\nIf strict is true, all elements in x[i-w:i+w] or the bounds of a plateau must exist and must not be missing or NaN. For strict == false, a minima is the minimum or first element of all existing and non-NaN or missing elements in x[i-w:i+w] or the bounds of a plateau.\n\nSee also: argminima\n\nExamples\n\njulia> findnextminima([3,2,3,1,1,3], 2)\n2\n\njulia> findnextminima([3,2,3,1,1,3], 3)\n4\n\n\n\n\n\n\n","category":"function"}]
}
