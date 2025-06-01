function exercise_Peaks()
    pks = findmaxima([0,5,2,3,3,1,4,0])
    pks = peakproms!(pks)
    pks = peakwidths!(pks)

    pks = findmaxima([0,5,2,3,missing,1,4,0])
    pks = peakproms!(pks)
    pks = peakwidths!(pks)

    pks = findmaxima([0,5,2,3,missing,1,4,0]; strict=false)
    pks = peakproms!(pks; strict=false)
    pks = peakwidths!(pks; strict=false)

    return nothing
end

