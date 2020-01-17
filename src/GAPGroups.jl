module GAPGroups

using GAP
using GAPTypes

import Base.==
import Base.rand

export symmetric_group, order, perm

struct GAPGroup
   X::GapObj
end

struct GAPGroupElem
   X::GapObj
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

function order(x::GAPGroup)
   return GAP.Globals.Size(x.X)
end

function order(x::GAPGroupElem)
   return GAP.Globals.Order(x.X)
end

function rand(x::GAPGroup)
   s=GAP.Globals.Random(x.X)
   return GAPGroupElem(s)
end

Base.:*(x::GAPGroupElem, y::GAPGroupElem) = GAPGroupElem(x.X * y.X)

function ==(x::GAPGroupElem, y::GAPGroupElem)
   return x.X == y.X
end

Base.:one(x::GAPGroup) = GAPGroupElem(GAP.Globals.Identity(x.X))

Base.:inv(x::GAPGroupElem) = GAPGroupElem(GAP.Globals.Inverse(x.X))

Base.:^(x::GAPGroupElem, y::Int64) = GAPGroupElem(x.X ^ y)

Base.:^(x::GAPGroupElem, y::GAPGroupElem) = GAPGroupElem(x.X ^ y.X)

Base.:<(x::GAPGroupElem, y::GAPGroupElem) = x.X < y.X

Base.:>(x::GAPGroupElem, y::GAPGroupElem) = x.X > y.X

Base.:/(x::GAPGroupElem, y::GAPGroupElem) = x*y^-1

# takes as input a list of arrays (not necessarly disjoint)
function perm(L::Array{Int64,1}...)
   if length(L)==0
      return one(symmetric_group(1))
   else
      return prod([GAPGroupElem(GAP.Globals.CycleFromList(GAP.julia_to_gap(y))) for y in L])
   end
end

Base.:sign(x::GAPGroupElem) = GAP.Globals.SignPerm(x.X)

Base.:isless(x::GAPGroupElem, y::GAPGroupElem) = x<y

#evaluation function
function (x::GAPGroupElem)(n)
   return GAP.Globals.OnPoints(n,x.X)
end

end
