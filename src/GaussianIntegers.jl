################################################################################
# This simple module illustrates how to implement an own ring using the
# AbstractAlgebra/Nemo computer algebra packages.
# The goal was to be able to compute Hermite and Smith normal forms of matrices
# over the Gaussian integers. To this end, we need to implement the ring and
# its Euclidean division with remainder. The rest is then handled by
# AbstractAlgebra.
#
# By Ulrich Thiel, 2019
# math@ulthiel.com
################################################################################

module GaussianIntegers

import Base
import AbstractAlgebra
import Nemo
import Nemo: fmpz, QQ, ZZ
import Hecke: round #needed for round(fmpq); not in Nemo

export GaussianIntegerRing

################################################################################
# Type definition of the Gaussian integer ring and its element.
# Both types are derived from the generic types in AbstractAlgebra.
################################################################################
mutable struct GaussianIntegerRing <: AbstractAlgebra.Ring end

mutable struct GaussianIntegerRingElem <: AbstractAlgebra.RingElem
  re::fmpz #real part
  im::fmpz #imaginary part
  P::GaussianIntegerRing #parent
end

# To connect these two structures, the following two functions need to be
# implemented
AbstractAlgebra.elem_type(::Type{GaussianIntegerRing}) = GaussianIntegerRingElem
AbstractAlgebra.parent_type(::Type{GaussianIntegerRingElem}) = GaussianIntegerRing

################################################################################
# Basic functions
################################################################################

# Element constructors
(R::GaussianIntegerRing)() = zero(R)
(R::GaussianIntegerRing)(a::Integer, b::Integer) = GaussianIntegerRingElem(ZZ(a),ZZ(b), R)
(R::GaussianIntegerRing)(a::fmpz, b::fmpz) = GaussianIntegerRingElem(a,b, R)
(R::GaussianIntegerRing)(a::Integer) = GaussianIntegerRingElem(ZZ(a),ZZ(0), R)
(R::GaussianIntegerRing)(x::GaussianIntegerRingElem) = x

# Printing functions
function Base.show(io::IO, R::GaussianIntegerRing)
  print(io, "Ring of Gaussian integers" )
end

function Base.show(io::IO, x::GaussianIntegerRingElem)
  print(io, (x.re, x.im) )
end

# This is just someting for printing
AbstractAlgebra.needs_parentheses(x::GaussianIntegerRingElem) = return false

# Parent of an element
AbstractAlgebra.parent(x::GaussianIntegerRingElem) = return x.P

# The base ring (I consider Z[i] as a Z-algebra, so the base ring is Z)
AbstractAlgebra.base_ring(R::GaussianIntegerRing) = return ZZ

# Yes, the parent is an integral domain
AbstractAlgebra.isdomain_type(x::GaussianIntegerRingElem) = return true

# Yes, every element is represented by exact means
AbstractAlgebra.isexact_type(x::GaussianIntegerRingElem) = return true

# A simple hash function for element (hash of real and imaginary part mixed)
Base.hash(x::GaussianIntegerRingElem, h::UInt) = return Base.hash(x.im, Base.hash(x.re, h))

# This, I don't know
#Base.deepcopy_internal = ??

# Random element creation
function AbstractAlgebra.rand(R::GaussianIntegerRing)
  return R(rand(Int), rand(Int))
end

################################################################################
# Arithmetic
################################################################################

Base.:(+)(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = GaussianIntegerRingElem(x.re+y.re, x.im+y.im, x.P)

Base.:(-)(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = GaussianIntegerRingElem(x.re-y.re, x.im-y.im, x.P)

Base.:(-)(x::GaussianIntegerRingElem) = GaussianIntegerRingElem(-x.re,-x.im, x.P)

Base.:(*)(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = GaussianIntegerRingElem(x.re*y.re-x.im*y.im, x.re*y.im+x.im*y.re, x.P)

AbstractAlgebra.mul!(a::GaussianIntegerRingElem, b::GaussianIntegerRingElem, c::GaussianIntegerRingElem) = b*c

AbstractAlgebra.add!(a::GaussianIntegerRingElem, b::GaussianIntegerRingElem, c::GaussianIntegerRingElem) = b+c

AbstractAlgebra.addeq!(a::GaussianIntegerRingElem, b::GaussianIntegerRingElem) = a+b

Base.:(==)(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = x.re == y.re && x.im == y.im

Base.isequal(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem) = x==y

Base.zero(R::GaussianIntegerRing) = GaussianIntegerRingElem(0,0,R)
Base.one(R::GaussianIntegerRing) = GaussianIntegerRingElem(1,0,R)

Base.iszero(x::GaussianIntegerRingElem) = iszero(x.re) && iszero(x.im)
Base.isone(x::GaussianIntegerRingElem) = isone(x.re) && iszero(x.im)

function AbstractAlgebra.isunit(x::GaussianIntegerRingElem)
  R = parent(x)
  u1 = R(1,0)
  u2 = R(0,1)
  return x==u1 || x==-u1 || x==u2 || x==-u2
end

function Base.inv(x::GaussianIntegerRingElem)
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

function Base.:(^)(x::GaussianIntegerRingElem, n::Integer)
  if n == 0
    return one(parent(x))
  elseif n == 1
    return x
  elseif n > 1
    y = x
    for i=2:n
      y = y*x
    end
    return y
  elseif n < 0
    y = inv(x)^abs(n)
  end
end

################################################################################
# Euclidean ring structure
################################################################################

# The key function is division with remainder
# Given x, y, we find q,r such that x = y*q + r and N(r) < N(y) or r=0
# This can be done using geometry (closest lattice point to x/y \in Q(i)).
# The fun fact is that this also fixes a system of residues, so I do not
# have to implement this.
function divrem(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem)
  if iszero(y)
    throw(DivideError())
  end

  #compute x/y=a+bi (rational a and b)
  a = QQ(x.re*y.re+x.im*y.im)//QQ( (y.re)^2 + (y.im)^2 )
  b = QQ(x.im*y.re - x.re*y.im)//QQ( (y.re)^2 + (y.im)^2 )

  #now, take closest integral point to (a,b)
  R = x.P
  q = R(round(a), round(b))

  #remainder
  r = x-q*y

  return (q,r)
end

# From divrem, we get all the other division functions
function Base.mod(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem)
  q,r = divrem(x,y)
  return r
end

function Base.div(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem)
  q,r = divrem(x,y)
  return q
end

# Wikipedia pseudo code for gcd
function Base.gcd(x::GaussianIntegerRingElem, y::GaussianIntegerRingElem)
  while !iszero(y)
    t = y
    y = mod(x,y)
    x = t
  end
  return x
end

# Wikipedia pseudo code for gcdx
function Base.gcdx(a::GaussianIntegerRingElem, b::GaussianIntegerRingElem)
  R = parent(a)
  s = zero(R)
  t = one(R)
  r = b
  old_s = one(R)
  old_t = zero(R)
  old_r = a

  while !iszero(r)
    quotient = div(old_r, r)

    prov = r
    r = old_r - quotient*r
    old_r = prov

    prov = s
    s = old_s - quotient*s
    old_s = prov

    prov = t
    t = old_t - quotient*t
    old_t = prov

  end

  return (old_r, old_s, old_t)

end

function AbstractAlgebra.divides(a::GaussianIntegerRingElem, b::GaussianIntegerRingElem)
  return iszero(mod(a,b))
end

function AbstractAlgebra.divexact(a::GaussianIntegerRingElem, b::GaussianIntegerRingElem)
  if iszero(b)
    throw(DivideError())
  end
  (q,r) = divrem(a,b)
  if !iszero(r)
    throw(DivideError())
  end
  return q
end

# As a last action: fix canonical forms of elements which can be reached by
# multiplication with a unit.
# I choose the first quadrant, so x=a+ib with a,b>=0.
# The canonical unit of x is the unit u such that xu^-1 is in the first quadrant.
function AbstractAlgebra.canonical_unit(x::GaussianIntegerRingElem)
  a = x.re
  b = x.im
  R = parent(x)
  if a >= 0 && b >= 0
    return R(1,0)
  elseif a >= 0 && b < 0
    return R(0,-1)
  elseif a < 0 && b >= 0
    return R(0,1)
  else
    return R(-1,0)
  end
end

end # module
