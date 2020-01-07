import Base.*
import Base.^
module PermModule

using GAP
export symmetric_group, MyRandom

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
      return PermModule.GAPGroup(GAP.Globals.SymmetricGroup(n))
   end
end

function order(x::PermModule.GAPGroup)
   return GAP.Globals.Size(x.X)
end

function rand(x::PermModule.GAPGroup)
   s=GAP.Globals.Random(x.X)
   return PermModule.GAPGroupElem(s)
end

function *(x::PermModule.GAPGroupElem, y::PermModule.GAPGroupElem)
   return PermModule.GAPGroupElem(x.X * y.X)
end

function ^(x::PermModule.GAPGroupElem, y::PermModule.GAPGroupElem)
   return PermModule.GAPGroupElem(x.X ^ y.X)
end

function perm(L::Array{Int64,1})
   z=GAP.Globals.CycleFromList(GAP.julia_to_gap(L))
   return PermModule.GAPGroupElem(z)
end
