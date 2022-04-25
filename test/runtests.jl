using GaussianIntegers
using AbstractAlgebra
using Test

@testset "GaussianIntegers" begin
  R = GaussianIntegerRing()

  # https://math.stackexchange.com/questions/82350/how-to-calculate-gcd-of-gaussian-integers
  x = R(11,7)
  y = R(18,-1)
  @test gcd(x,y) == one(R)
  (g,a,b) = gcdx(x,y)
  @test g == one(R)
  @test a*x + b*y == g

  # My example from class
  A=matrix(R,2,2,[(2,-1), (2,0), (7,-1), (3,1)])
  @test hnf(A)==matrix(R,2,2,[(1,2), (1,-1), (0,0), (3,1)])
  @test snf(A)==matrix(R,2,2,[(1,0), (0,0), (0,0), (1,7)])

  # From https://www.maplesoft.com/support/help/Maple/view.aspx?path=GaussInt/GIhermite&cid=411
  A = matrix(R, 2, 3, [(3,-7), (7, 11), (0,11), (13,-4), (17,12), (19,0) ])
  @test hnf(A) == matrix(R,2,3,[(1,0),(-59,-2),(-82,-8),(0,0),(198,0),(276,13)])

  # From https://www.maplesoft.com/support/help/Maple/view.aspx?path=GaussInt%2fGIsmith
  A = matrix(R,3,3,[(-4,7),(8,10),(-6,-8),(-5,7),(6,-6),(0,5),(-10,1),(1,-3),(-10,5)])
  @test snf(A) == matrix(R,3,3,[(1,0),(0,0),(0,0),(0,0),(1,0),(0,0),(0,0),(0,0),(1797,791)])
end
