using Peaks: check_known_fields_equal_length, check_has_required_fields
@testset "Utility functions" begin
    @test_throws DimensionMismatch check_known_fields_equal_length((;indices=rand(4), heights=rand(3)))
    @test_throws ArgumentError check_has_required_fields((;nonstandard=1))

    @test_throws DimensionMismatch filterpeaks!((;indices=rand(4), heights=rand(4)), rand(Bool, 5))
    @test_throws ArgumentError filterpeaks!((;indices=rand(4), heights=rand(4)), :proms; max=1)

end
