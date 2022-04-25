# GaussianIntegers.jl

This simple [Julia](https://julialang.org) module illustrates how to implement the ring of [Gaussian integers](https://en.wikipedia.org/wiki/Gaussian_integer) using the generic (Euclidean) ring interface of [AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl). This allows for example to compute Hermite and Smith normal forms of matrices over the Gaussian integers without further ado.

By [Ulrich Thiel](https://ulthiel.com/math) (University of Kaiserslautern), 2019 (updated 2022)

## Background

In my algebraic number theory course I (un)fortunately gave as an exercise to compute Hermite and Smith normal forms of matrices over the Gaussian integers to illustrate that the theory works over Euclidean rings in general and not just over ℤ and K[X].[^1] I noticed that by hand it is easy to make mistakes and therefore I wanted to verify my results with the computer. But if one creates the Gaussian integers as a maximal order in a computer algebra system like Magma, the system doesn't know about the Euclidean ring structure, so one can't compute division with remainder and normal forms of matrices. AbstractAlgebra.jl supports generic (Euclidean) rings and algorithms, and once all the basic ring functions are implemented, everything else is handled without further ado.

**Remark.** The Hermite normal form is only unique after fixing a choice of a system of representatives of non-associates and of residues. In the implementation I have chosen the points in the first quadrant as representatives of non-associates—this is handled by the function ```canonical_unit```. The division with remainder algorithm (which I do geometrically) *by itself* defines a system of representatives of residues, so we don't (have to) deal with this any further.

## Installation
```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/ulthiel/GaussianIntegers.jl")
```

## Usage

```julia
julia> using GaussianIntegers, AbstractAlgebra

julia> R = GaussianIntegerRing() #Create the ring of Gaussian integers
Ring of Gaussian integers

julia> x=R(2,1) #Create the Gaussian integer 2+i*1
(2, 1)

julia> x+x #Addition
(4, 2)

julia> x*x #Multiplication
(3, 4)

julia> A=matrix(R,2,2,[(2,-1), (2,0), (7,-1), (3,1)]) #Creating a 2x2-matrix
[(2, -1)  (2, 0)]
[(7, -1)  (3, 1)]

julia> hnf(A) #The Hermite normal form of A (hnf_with_transform will also return the transformation matrix)
[(1, 2)  (1, -1)]
[(0, 0)   (3, 1)]

julia> snf(A) #The Smith normal form of A
[(1, 0)  (0, 0)]
[(0, 0)  (1, 7)]
```

## References

[^1]: Adkins, W. A. & Weintraub, S. H. (1992). *Algebra*). Chapter 5.
