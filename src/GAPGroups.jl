using GAP
import Base.*
import Base.^
import Base.inv
import Base.==
import Base.hash
import Base.<
import Base.>
import Base.isless
# export symmetric_group, hash, order, rand, *, ==, identity, inv, ^, >, <, perm

module GAPGroups

Tg=Main.ForeignGAP.MPtr

struct GAPGroup
   X::Tg
end

struct GAPGroupElem
   X::Tg
end

end


function symmetric_group(n::Int64)
   if n<1
     throw(ArgumentError("it must be a positive integer"))
   else
      return GAPGroups.GAPGroup(GAP.Globals.SymmetricGroup(n))
   end
end

# to be fixed later
function hash(x::GAPGroups.GAPGroupElem)
   return 0
end

function order(x::GAPGroups.GAPGroup)
   return GAP.Globals.Size(x.X)
end

function order(x::GAPGroups.GAPGroupElem)
   return GAP.Globals.Order(x.X)
end

function rand(x::GAPGroups.GAPGroup)
   s=GAP.Globals.Random(x.X)
   return GAPGroups.GAPGroupElem(s)
end

function *(x::GAPGroups.GAPGroupElem, y::GAPGroups.GAPGroupElem)
   return GAPGroups.GAPGroupElem(x.X * y.X)
end

function ==(x::GAPGroups.GAPGroupElem, y::GAPGroups.GAPGroupElem)
   return x.X == y.X
end

function identity(x::GAPGroups.GAPGroup)
   return GAPGroups.GAPGroupElem(GAP.Globals.Identity(x.X))
end

function inv(x::GAPGroups.GAPGroupElem)
   return GAPGroups.GAPGroupElem(GAP.Globals.Inverse(x.X))
end


function ^(x::GAPGroups.GAPGroupElem, y::Int64)
   return GAPGroups.GAPGroupElem(x.X ^ y)
end

function ^(x::GAPGroups.GAPGroupElem, y::GAPGroups.GAPGroupElem)
   return GAPGroups.GAPGroupElem(x.X ^ y.X)
end

function <(x::GAPGroups.GAPGroupElem, y::GAPGroups.GAPGroupElem)
   return x.X < y.X
end

function >(x::GAPGroups.GAPGroupElem, y::GAPGroups.GAPGroupElem)
   return x.X > y.X
end

function perm(L::Array{Int64,1})
   z=GAP.Globals.CycleFromList(GAP.julia_to_gap(L))
   return GAPGroups.GAPGroupElem(z)
end

function isless(x::GAPGroups.GAPGroupElem, y::GAPGroups.GAPGroupElem)
   return x<y
end
