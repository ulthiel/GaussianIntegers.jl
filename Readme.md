# GaussianIntegers.jl

This simple [Julia](https://julialang.org) module illustrates how to implement an own ring using the
[AbstractAlgebra](https://github.com/Nemocas/AbstractAlgebra.jl)/[Nemo](https://github.com/Nemocas/Nemo.jl) computer algebra packages.

I wanted to compute Hermite and Smith normal forms of matrices
over the Gaussian integers (I (un-)fortunately gave an exercise in class about this). The problem is that if you create the Gaussian integers as a maximal order in a computer algebra system, it doesn't know about the Euclidean ring structure, so you can't compute division with remainder and normal forms. AbstractAlgebra supports generic ring types and algorithms, and once you implement all the necessary functions, everything else is handled by AbstractAlgebra.

As a system of representatives of non-associates I have chosen the points in the first quadrant (this is handled by the function ```canonical_unit``` of AbstractAlgebra). The division with remainder (which I do geometrically) automatically yields a system of representatives of residues and I (don't have to) deal with this any further.

## Installation
```julia
julia> using Pkg

julia> Pkg.add(PackageSpec(url="https://github.com/ulthiel/GaussianIntegers.jl", rev="master" ))
```

## Usage

```julia
julia> using AbstractAlgebra
julia> using GaussianIntegers

julia> R = GaussianIntegerRing() #create the ring of Gaussian integers
Ring of Gaussian integers

julia> x=R(2,1) #Create the Gaussian integer 2+i*1
(2, 1)

julia> x+x #Addition
(4, 2)

julia> x*x #Multiplication
(3, 4)

julia> A=matrix(R,2,2,[R(2,-1), R(2,0), R(7,-1), R(3,1)]) #creating a 2x2-matrix
[(2, -1)  (2, 0)]
[(7, -1)  (3, 1)]

julia> hnf(A) #the Hermite normal form of A (hnf_with_transform will also return the transformation matrix)
[(1, 2)  (1, -1)]
[(0, 0)   (3, 1)]

julia> snf(A) #the Smith normal form of A
[(0, 1)  (0, 0)]
[(0, 0)  (1, 7)]
```
