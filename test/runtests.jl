using Peaks
using Test, Aqua, JET, ExplicitImports, ReferenceTests, ImageIO
using OffsetArrays, Plots, CairoMakie

@testset verbose=true "Peaks" begin
    @testset "Aqua tests" begin
        Aqua.test_all(Peaks)
    end
    @testset "Explicit and public imports" begin
        @test check_no_implicit_imports(Peaks; ignore=(
            Symbol("@recipe"), Symbol("@series"), Symbol("@shorthands")
            )) === nothing
        @test check_no_stale_explicit_imports(Peaks) === nothing
        @test check_all_explicit_imports_via_owners(Peaks) === nothing


        @static if VERSION < v"1.11"
            @test check_all_qualified_accesses_are_public(Peaks; ignore=(
                :_overflowind, :_blsr, :_toind, # for findall_offset function
                :Experimental, :register_error_hint, :argument_names, # should be public
                :VecTypes, # from SIMD
                :Fix2, :depwarn
                )) === nothing
        else
            @test check_all_qualified_accesses_are_public(Peaks; ignore=(
                :_overflowind, :_blsr, :_toind, # for findall_offset function
                :Experimental, :register_error_hint, :argument_names, # should be public
                :VecTypes, # from SIMD
                )) === nothing
        end
    end
    @static if isempty(VERSION.prerelease)
        @testset "JET inference tests" begin
            JET.test_package(Peaks; target_modules=(Peaks,))
        end
    end

    include("simd.jl")
    include("minmax.jl")
    include("peakprom.jl")
    include("peakwidth.jl")
    include("peakheight.jl")
    include("plotting.jl")
    include("utils.jl")
end
