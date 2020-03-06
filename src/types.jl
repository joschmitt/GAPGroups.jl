abstract type Group end

abstract type GroupElem end

export PermGroup, PermGroupElem, MatrixGroup, MatrixGroupElem, PcGroup, PcGroupElem, 
       FPGroup, FPGroupElem, AutomorphismGroup, AutomorphismGroupElem

struct PermGroup <: Group
   X::GapObj
   deg::Int64       # G < Sym(deg)
   
   function PermGroup(G::GapObj)
     @assert GAP.Globals.IsPermGroup(G)
     n = GAP.gap_to_julia(Int64, GAP.Globals.LargestMovedPoint(G))
     z = new(G, n)
     return z
   end
   
   function PermGroup(G::GapObj, deg::Int)
     @assert GAP.Globals.IsPermGroup(G)
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
  function MatrixGroup(G::GapObj)
    @assert GAP.Globals.IsMatrixGroup(G)
    z = new(G)
    return z
  end
end

mutable struct MatrixGroupElem <: GroupElem
   parent::MatrixGroup
   X::GapObj       
end

struct PcGroup <: Group
  X::GapObj
  function PcGroup(G::GapObj)
    @assert GAP.Globals.IsPcGroup(G)
    z = new(G)
    return z
  end
end

mutable struct PcGroupElem <: GroupElem
   parent::PcGroup
   X::GapObj
end

struct FPGroup <: Group
  X::GapObj
  
  function FPGroup(G::GapObj)
    @assert GAP.Globals.IsFpGroup(G)
    z = new(G)
    return z
  end
end

mutable struct FPGroupElem <: GroupElem
   parent::FPGroup
   X::GapObj
end

mutable struct AutomorphismGroup{T} <: Group
  X::GapObj
  G::T
  function AutomorphismGroup{T}(G::GapObj, H::T) where T
    @assert GAP.Globals.IsAutomorphismGroup(G)
    z = new{T}(G, H)
    return z
  end
end

mutable struct AutomorphismGroupElem{T} <: GroupElem
   parent::AutomorphismGroup{T}
   X::GapObj
end


const gap_group_types = []

function _get_type(G::GapObj)
  for i = 1:length(gap_group_types)
    if gap_group_types[i][1](G)
      return gap_group_types[i][2]
    end
  end
  error("Not a known type of group")
end

function _get_gap_function(T)
  for i = 1:length(gap_group_types)
    if gap_group_types[i][2] == T
      return gap_group_types[i][1]
    end
  end
  error("Not a known type of group")
end


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
