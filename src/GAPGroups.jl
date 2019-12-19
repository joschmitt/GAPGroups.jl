module PermModule

using GAP
export symmetric_group, MyRandom

Tg=Main.ForeignGAP.MPtr

struct GAPGroup
   X::Tg
end

struct GAPGroupElem
   x::Tg
end

function MyOrder(x::GAPGroup)
   return GAP.Globals.Size(x.X)
end

function MyRandom(x::GAPGroup)
   s=GAP.Globals.Random(x.X)
   y=GAPGroupElem(s)
   return y
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
   return PermModule.MyOrder(x)
end

function random(x::PermModule.GAPGroup)
   return PermModule.MyRandom(x)
end
