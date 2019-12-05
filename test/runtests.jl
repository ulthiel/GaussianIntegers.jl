using Test
using GaussianIntegers
using AbstractAlgebra

@testset "GaussianIntegers" begin
  R = GaussianIntegerRing()
  A=matrix(R,2,2,[R(2,-1), R(2,0), R(7,-1), R(3,1)])
  @test hnf(A)==matrix(R,2,2,[R(1,2), R(1,-1), R(0,0), R(3,1)])
  @test snf(A)==matrix(R,2,2,[R(0,1), R(0,0), R(0,0), R(1,7)])
end
