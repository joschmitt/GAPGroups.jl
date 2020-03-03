abstract type Group end

abstract type GroupElem end

struct PermGroup <: Group
   X::GapObj
   deg::Int64       # G < Sym(deg)
   
   function PermGroup(G::GapObj)
     n = GAP.gap_to_julia(Int64, GAP.Globals.LargestMovedPoint(G))
     z = new(G, n)
     return z
   end
   
   function PermGroup(G::GapObj, deg::Int)
     z = new(G, deg)
     return z
   end
end

mutable struct PermGroupElem <: GroupElem
   parent::PermGroup
   X::GapObj   
end

struct MatrixGroup <: Group
  X::GapObj
end

mutable struct MatrixGroupElem <: GroupElem
   parent::MatrixGroup
   X::GapObj       
end

struct PolycyclicGroup <: Group
  X::GapObj
end

mutable struct PolycyclicGroupElem <: GroupElem
   parent::PolycyclicGroup
   X::GapObj
end

struct FPGroup <: Group
  X::GapObj
end

mutable struct FPGroupElem <: GroupElem
   parent::FPGroup
   X::GapObj
end


const gap_group_types = 
[(GAP.Globals.IsPermGroup, PermGroup), (GAP.Globals.IsPcGroup, PolycyclicGroup), 
 (GAP.Globals.IsMatrixGroup, MatrixGroup), (GAP.Globals.IsFpGroup, FPGroup)               
]


################################################################################
#
#  Group Homomorphism
#
################################################################################

mutable struct GAPGroupHomomorphism{S, T}
   domain::S
   codomain::T
   image::GapObj
   
   function GAPGroupHomomorphism{S, T}(domain::S, codomain::T, image::GapObj) where S <: Group where T <: Group
     z = new{S, T}()
     z.domain = domain
     z.codomain = codomain
     z.image = image
     return z
   end
end
