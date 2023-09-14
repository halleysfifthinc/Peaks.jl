# I am putting everything in here for now. The contents of this file should be moved around in the future.
"""
    Docs here
"""
function findpeaks(x)
    data = x
    inds, _ = findmaxima(x)
    return (data=x, inds=inds)
end