# Unit tests for MembraneAleFem.jl
#
# References:
# https://www.matecdev.com/posts/julia-testing.html
# https://docs.julialang.org/en/v1.10/stdlib/Test/#Testing-Our-Package
#
# Written by Amaresh Sahu (asahu@che.utexas.edu) on 17 January 2024

using Test, MembraneAleFem
using LinearAlgebra, DelimitedFiles, StaticArrays, Parameters


@testset "MembraneAleFem.jl tests" begin

  @testset "Input tests" begin
    include("../src/Input.jl");

    @testset "Spline.jl tests" begin
      include("input/spline_tests.jl");
    end
    @testset "GaussPoint.jl tests" begin
      include("input/gauss_pt_tests.jl");
    end
    @testset "GpBasisFn.jl tests" begin
      include("input/gp_basis_fn_tests.jl");
    end
  end

  @testset "Analysis tests" begin
    include("analysis/geo_dyn_stress_tests.jl");
  end

end
