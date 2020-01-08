using GAP
import Base.*
import Base.^
import Base.inv


module PermModule

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

function order(x::PermModule.GAPGroupElem)
   return GAP.Globals.Order(x.X)
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

function identity(x::PermModule.GAPGroup)
   return PermModule.GAPGroupElem(GAP.Globals.Identity(x.X))
end

function inv(x::PermModule.GAPGroupElem)
   return PermModule.GAPGroupElem(GAP.Globals.Inverse(x.X))
end

function ^(x::PermModule.GAPGroupElem, y::Int64)
   if y<0
      return PermModule.GAPGroupElem(inv(x).X^(-y))
   else
      return PermModule.GAPGroupElem(x.X ^ y)
   end
end

function perm(L::Array{Int64,1})
   z=GAP.Globals.CycleFromList(GAP.julia_to_gap(L))
   return PermModule.GAPGroupElem(z)
end
