using GAP
import Base.*
import Base.^
import Base.inv
import Base.==
import Base.hash
import Base.<
import Base.>
import Base.isless

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

# to be fixed later
function hash(x::PermModule.GAPGroupElem)
   return 0
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

function ==(x::PermModule.GAPGroupElem, y::PermModule.GAPGroupElem)
   return x.X == y.X
end

function identity(x::PermModule.GAPGroup)
   return PermModule.GAPGroupElem(GAP.Globals.Identity(x.X))
end

function inv(x::PermModule.GAPGroupElem)
   return PermModule.GAPGroupElem(GAP.Globals.Inverse(x.X))
end


function ^(x::PermModule.GAPGroupElem, y::Int64)
   return PermModule.GAPGroupElem(x.X ^ y)
end

function ^(x::PermModule.GAPGroupElem, y::PermModule.GAPGroupElem)
   return PermModule.GAPGroupElem(x.X ^ y.X)
end

function <(x::PermModule.GAPGroupElem, y::PermModule.GAPGroupElem)
   return x.X < y.X
end

function >(x::PermModule.GAPGroupElem, y::PermModule.GAPGroupElem)
   return x.X > y.X
end

function perm(L::Array{Int64,1})
   z=GAP.Globals.CycleFromList(GAP.julia_to_gap(L))
   return PermModule.GAPGroupElem(z)
end

function isless(x::PermModule.GAPGroupElem, y::PermModule.GAPGroupElem)
   return x<y
end
