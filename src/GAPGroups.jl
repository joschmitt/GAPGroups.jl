# further possible functions: similar, iterate, literal_pow, parent_type

module GAPGroups

using GAP

import Base.==
import Base.rand
import Base.conj
import Base.conj!
import Base.parent
import Base.eltype

export symmetric_group, order, perm, cperm, hasorder, hasgens, gens, ngens, comm, comm!, inv!, rand_pseudo, one!, div_right, div_left, div_right!, div_left!, elem_type, mul, mul!, deg
     #conj!, conj

struct GAPGroup
   X::GapObj
end

struct GAPGroupElem
   X::GapObj
   par::Int64       # Sym(par) = parent of the object
end

function symmetric_group(n::Int64)
   if n<1
     throw(ArgumentError("it must be a positive integer"))
   else
      return GAPGroup(GAP.Globals.SymmetricGroup(n))
   end
end

# to be fixed later
Base.:hash(x::GAPGroupElem) = 0

#shorter writing
#Base.:parent(x::GAPGroupElem) = x==one(symmetric_group(2)) ? symmetric_group(1) : symmetric_group(GAP.Globals.LargestMovedPointPerm(x.X))

function parent(x::GAPGroupElem)
#   if x==one(symmetric_group(1))
#      n = 1
#   else
#      n = GAP.Globals.LargestMovedPointPerm(x.X)
#   end
   return symmetric_group(x.par)
end

Base.:eltype(::Type{G}) where G<:GAPGroup = GAPGroup
Base.:eltype(::Type{G}) where G<:GAPGroupElem = GAPGroupElem

function order(x::GAPGroup)
   return GAP.Globals.Size(x.X)
end

function order(x::GAPGroupElem)
   return GAP.Globals.Order(x.X)
end

function order(::Type{T}, x::Union{GAPGroupElem, GAPGroup}) where T<:Number
   return T(order(x))
end

# degree of the Permutation Group
deg(x::GAPGroup) = GAP.Globals.LargestMovedPoint(x.X)

function rand(x::GAPGroup)
   s=GAP.Globals.Random(x.X)
   return GAPGroupElem(s,deg(x))
end

# one of the following should be non-parametric
elem_type(G::GAPGroup) = GAPGroupElem
rand_pseudo(G::GAPGroup) = rand(G)

Base.:*(x::GAPGroupElem, y::GAPGroupElem) = GAPGroupElem(x.X * y.X, max(x.par, y.par))

function ==(x::GAPGroup, y::GAPGroup)
   return x.X == y.X
end

function ==(x::GAPGroupElem, y::GAPGroupElem)
   return x.X == y.X
end

Base.:one(x::GAPGroup) = GAPGroupElem(GAP.Globals.Identity(x.X),deg(x))
Base.:one(x::GAPGroupElem) = one(parent(x),x.par)
one!(x::GAPGroupElem) = one(parent(x),x.par)

Base.:isone(x::GAPGroupElem) = x == one(parent(x))

Base.:inv(x::GAPGroupElem) = GAPGroupElem(GAP.Globals.Inverse(x.X),x.par)

inv!(out::GAPGroupElem, x::GAPGroupElem) = inv(x)  #if needed later

Base.:^(x::GAPGroupElem, y::Integer) = GAPGroupElem(x.X ^ y,x.par)

Base.:^(x::GAPGroupElem, y::GAPGroupElem) = GAPGroupElem(x.X ^ y.X, max(x.par,y.par))

Base.:<(x::GAPGroupElem, y::GAPGroupElem) = x.X < y.X

Base.:/(x::GAPGroupElem, y::GAPGroupElem) = x*y^-1

mul(x::GAPGroupElem, y::GAPGroupElem) = x*y
mul!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x*y

div_right(x::GAPGroupElem, y::GAPGroupElem) = x*inv(y)
div_left(x::GAPGroupElem, y::GAPGroupElem) = inv(y)*x
div_right!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x*inv(y)
div_left!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = inv(y)*x

conj(x::GAPGroupElem, y::GAPGroupElem) = x^y
conj!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x^y

comm(x::GAPGroupElem, y::GAPGroupElem) = x^-1*x^y
comm!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x^-1*x^y

#maybe in future add more checks
hasorder(x::GAPGroup) = true
hasorder(x::GAPGroupElem) = true
hasgens(x::GAPGroup) = true

function perm(L::Array{Int64,1})
   return GAPGroupElem(GAP.Globals.PermList(GAP.julia_to_gap(L)), length(L))
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of arrays (not necessarly disjoint)
function cperm(L::Union{Array{Int64,1},UnitRange{Int64}}...)
   if length(L)==0
      return one(symmetric_group(1))
   else
      return prod([GAPGroupElem(GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y))), maximum(y)) for y in L])
   end
end

function gens(G::GAPGroup)
   L=GAP.Globals.GeneratorsOfGroup(G.X)
   l=length(L)
   return [GAPGroupElem(L[i],deg(G)) for i in 1:l]
end

function gens(G::GAPGroup, i::Integer)
   L=GAP.Globals.GeneratorsOfGroup(G.X)
   return GAPGroupElem(L[i],deg(G))
end

ngens(G::GAPGroup) = length(gens(G))

Base.:sign(x::GAPGroupElem) = GAP.Globals.SignPerm(x.X)

Base.:isless(x::GAPGroupElem, y::GAPGroupElem) = x<y

#evaluation function
function (x::GAPGroupElem)(n)
   return GAP.Globals.OnPoints(n,x.X)
end

end
