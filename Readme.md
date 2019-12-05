# GaussianIntegers

This simple [Julia](https://julialang.org) module illustrates how to implement an own ring using the
[AbstractAlgebra](https://github.com/Nemocas/AbstractAlgebra.jl)/[Nemo](https://github.com/Nemocas/Nemo.jl) computer algebra packages.

I wanted to compute Hermite and Smith normal forms of matrices
over the Gaussian integers (I unfortunately gave an exercise in class about this). The problem is that if you create the Gaussian integers as a maximal order in a computer algebra system, it doesn't know about the Euclidean ring structure, so you can't compute normal forms. AbstractAlgebra supports generic ring types and algorithms, and once you implement all the necessary functions, everything else is handled by the generic algorithms in AbstractAlgebra.

## Usage

```julia
julia> using AbstractAlgebra
julia> using GaussianIntegers

julia> R = GaussianIntegerRing() #create the ring of Gaussian integers
Ring of Gaussian integers

julia> x=R(2,1) #Create the element 2+i*1
(2, 1)

julia> x+x #Addition of elements
(4, 2)

julia> x*x #Multiplication of elements
(3, 4)

julia> A=matrix(R,2,2,[R(2,-1), R(2,0), R(7,-1), R(3,1)]) #a 2x2-matrix
[(2, -1)  (2, 0)]
[(7, -1)  (3, 1)]

julia> hnf(A) #the Hermite normal form of A
[(1, 2)  (1, -1)]
[(0, 0)   (3, 1)]

julia> snf(A) #the Smith normal form of A
[(0, 1)  (0, 0)]
[(0, 0)  (1, 7)]
```