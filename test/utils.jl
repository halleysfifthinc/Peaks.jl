using Peaks: check_known_fields_equal_length, check_has_required_fields
@testset "Utility functions" begin
    @test_throws DimensionMismatch check_known_fields_equal_length((;indices=rand(4), heights=rand(3)))
    @test_throws ArgumentError check_has_required_fields((;nonstandard=1))

    @test_throws DimensionMismatch filterpeaks!((;indices=rand(4), heights=rand(4)), rand(Bool, 5))
    @test_throws ArgumentError filterpeaks!((;indices=rand(4), heights=rand(4)), :proms; max=1)

    let
        x = [0,5,2,3,3,1,4,0]
        pks = findmaxima(x)
        inds, heights = peakheights(pks.indices, pks.heights; max=4)
        pksh = peakheights(pks; max=4)
        pkshfilt = filterpeaks!(peakheights(pks), :heights; max=4)
        pkshfilt2 = filterpeaks!(peakheights(pks), Bool[0, 1, 1])
        pkshfilt3 = filterpeaks!(peakheights(pks)) do pk
            pk.heights ≤ 4
        end

        @test heights == pksh.heights
        @test heights == pkshfilt.heights
        @test heights == pkshfilt2.heights
        @test heights == pkshfilt3.heights
    end

end
