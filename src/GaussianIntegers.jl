################################################################################
# This simple Julia module illustrates how to implement the ring of Gaussian
# integers using the generic (Euclidean) ring interface of AbstractAlgebra.jl.
# This allows for example to compute Hermite and Smith normal forms of matrices
# over the Gaussian integers without further ado.
#
# By Ulrich Thiel (University of Kaiserslautern), 2019 (updated 2022)
################################################################################

module GaussianIntegers

################################################################################
# Imports
################################################################################
import Base:
    show, hash, +, -, *, //, ==, ^, isequal, zero, one, iszero, isone,
    inv, deepcopy_internal, divrem

import AbstractAlgebra:
    ZZ, Ring, RingElem, elem_type, parent_type, parent, base_ring,
    isdomain_type, isexact_type, rand, mul!, add!, addeq!, isunit,
    canonical_unit, characteristic, divexact

# Note: In AbstractAlgebra ZZ is shortcut for BigInt.

################################################################################
# Exports
################################################################################
export GaussianIntegerRing

################################################################################
# Now, the actual implementation starts.
# For the description of the ring interface of AbstractAlgebra.jl see
# https://nemocas.github.io/AbstractAlgebra.jl/stable/ring_interface/.
# Note: We do not have to implement every function mentioned there since some
# follow from others generically, e.g. once we have divrem, we get mod,
# gcd, etc. generically (one could implement those if there are faster ways to
# determine the result but this is not necessary here).
################################################################################

################################################################################
# Ring and ring element type
################################################################################
mutable struct GaussianIntegerRing <: Ring end

mutable struct GaussianIntegerRingElem <: RingElem
    re::BigInt #real part
    im::BigInt #imaginary part
    parent::GaussianIntegerRing #parent
end

################################################################################
# Data type and parent object methods
################################################################################

# Parent type
parent_type(::Type{GaussianIntegerRingElem}) = GaussianIntegerRing

# Element type
elem_type(::Type{GaussianIntegerRing}) = GaussianIntegerRingElem

# The base ring (we consider Z[i] as a Z-algebra, so the base ring is Z)
base_ring(R::GaussianIntegerRing) = return ZZ

# Parent of an element
parent(x::GaussianIntegerRingElem) = return x.parent

# Yes, the parent is an integral domain
isdomain_type(x::GaussianIntegerRingElem) = return true

# Yes, every element is represented by exact means
isexact_type(x::GaussianIntegerRingElem) = return true

# A simple hash function for element (hash of real and imaginary part mixed)
hash(x::GaussianIntegerRingElem, h::UInt) = return hash(x.im, hash(x.re, h))

# Deepcopy
function deepcopy_internal(x::GaussianIntegerRingElem)
    R = x.parent
    return R(deepcopy_internal(x.re),deepcopy_internal(x.im))
end

################################################################################
# Constructors
################################################################################
(R::GaussianIntegerRing)(a::Integer, b::Integer) = GaussianIntegerRingElem(ZZ(a),ZZ(b), R)
(R::GaussianIntegerRing)() = R(0,0)
(R::GaussianIntegerRing)(a::Integer) = R(ZZ(a),ZZ(0))
(R::GaussianIntegerRing)(x::Tuple{Integer, Integer}) = R(ZZ(x[1]), ZZ(x[2]))

function (R::GaussianIntegerRing)(x::GaussianIntegerRingElem)
    return R(x.re, x.im)
end

################################################################################
# Basic manipulation
################################################################################
zero(R::GaussianIntegerRing) = R(0,0)
one(R::GaussianIntegerRing) = R(1,0)

iszero(x::GaussianIntegerRingElem) = iszero(x.re) && iszero(x.im)
isone(x::GaussianIntegerRingElem) = isone(x.re) && iszero(x.im)

################################################################################
# Canonicalisation
################################################################################

# This stuff is about choosing a system of representatives of ring elements
# under the associates relation (i.e. x and y are equivalent if y = xu^{-1} for
# a unit u). A system of representatives for Z[i] is the first quadrant (more
# precisely {0} union {re > 0 and im >= 0}).
# Given an element x, its *canonical unit* (in AbstractAlgebra.jl) is the unit u
# such that xu^{-1} is in the fixed system of representatives under the
# associates relation.
# We move from one quadrant to the quadrant before by multiplication by -i
# (rotation by 90 degrees clockwise).
function canonical_unit(x::GaussianIntegerRingElem)
    a = x.re
    b = x.im
    R = parent(x)
    if a == 0 && b == 0 #origin
        return one(R)
    elseif a > 0 && b >= 0 #I
        return one(R)
    elseif a <=0 && b > 0 #II
        return inv(R(0,-1))
    elseif a < 0 && b <= 0 #III
        return inv(R(-1,0))
    else
        return inv(R(0,1))
    end
end

################################################################################
# String I/O
################################################################################
function show(io::IO, R::GaussianIntegerRing)
    print(io, "Ring of Gaussian integers" )
end

function show(io::IO, x::GaussianIntegerRingElem)
    print(io, (x.re, x.im) )
end

################################################################################
# Expressions
################################################################################

# This is just for printing. I'm too lazy and this is also not really necessary.

################################################################################
# Unary operations
################################################################################
(-)(x::GaussianIntegerRingElem) = GaussianIntegerRingElem(-x.re,-x.im, x.parent)

################################################################################
# Binary operations
################################################################################
(+)(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = GaussianIntegerRingElem(x.re+y.re, x.im+y.im, x.parent)

(-)(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = GaussianIntegerRingElem(x.re-y.re, x.im-y.im, x.parent)

(*)(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = GaussianIntegerRingElem(x.re*y.re-x.im*y.im, x.re*y.im+x.im*y.re, x.parent)

################################################################################
# Comparison
################################################################################
(==)(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = x.re == y.re && x.im == y.im

isequal(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = x==y

################################################################################
# Exact division
################################################################################

# This uses the function divrem implemented below
# (There's no generic code yet)
function divexact(a::GaussianIntegerRingElem, b::GaussianIntegerRingElem; check::Bool=true)
    if iszero(b)
        throw(DivideError())
    end
    (q,r) = divrem(a,b)
    if check && !iszero(r)
        throw(DivideError())
    end
    return q
end

################################################################################
# Inverse
################################################################################
# There are just 4 units, so this is easy.
function inv(x::GaussianIntegerRingElem)
    R = parent(x)
    u1 = R(1,0)
    u2 = R(0,1)
    if x==u1
        return u1
    elseif x==-u1
        return -u1
    elseif x==u2
        return -u2
    elseif x==-u2
        return u2
    else
        throw(DivideError())
    end
end

################################################################################
# Unsafe operators
################################################################################
function zero!(x::GaussianIntegerRingElem)
    x.re = ZZ(0)
    x.im = ZZ(0)
    return x
end

function mul!(c::GaussianIntegerRingElem, a::GaussianIntegerRingElem, b::GaussianIntegerRingElem)
    c.re = a.re*b.re-a.im*b.im
    c.im = a.re*b.im+a.im*b.re
    return c
end

function add!(c::GaussianIntegerRingElem, a::GaussianIntegerRingElem, b::GaussianIntegerRingElem)
    c.re = a.re+b.re
    c.im = a.im+b.im
    return c
end

function addeq!(a::GaussianIntegerRingElem, b::GaussianIntegerRingElem)
    a.re = a.re+b.re
    a.im = a.im+b.im
    return a
end

################################################################################
# Random generation
################################################################################
rand(R::GaussianIntegerRing) = return R(rand(Int), rand(Int))

################################################################################
# Optional basic manipulation functionality
################################################################################
function isunit(x::GaussianIntegerRingElem)
    R = parent(x)
    u1 = R(1,0)
    u2 = R(0,1)
    return x==u1 || x==-u1 || x==u2 || x==-u2
end

characteristic(R::GaussianIntegerRing) = return 0

################################################################################
# Now, we come to the euclidean ring structure.
# https://nemocas.github.io/AbstractAlgebra.jl/stable/euclidean_interface/.
################################################################################

# The key function is division with remainder.
# Given x, y, we find q,r such that x = y*q + r and N(r) < N(y) or r=0
# This can be done using geometry (closest lattice point to x/y \in Q(i)).
# The fun fact is that this also fixes a system of residues, so we do not
# have to implement any of this to make HNF unique.
function divrem(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem)
    if iszero(y)
        throw(DivideError())
    end

    #compute x/y=a+bi (rational a and b)
    a = (x.re*y.re+x.im*y.im)//( (y.re)^2 + (y.im)^2 )
    b = (x.im*y.re - x.re*y.im)//( (y.re)^2 + (y.im)^2 )

    #now, take closest integral point to (a,b)
    R = x.parent
    q = R(ZZ(round(a)), ZZ(round(b)))

    #remainder
    r = x-q*y

    return (q,r)
end

end #module
