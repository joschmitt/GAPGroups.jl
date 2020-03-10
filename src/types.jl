"""
TODO: document this
"""
abstract type Group end

"""
TODO: document this
"""
abstract type GroupElem end

export PermGroup, PermGroupElem, MatrixGroup, MatrixGroupElem, PcGroup, PcGroupElem, 
       FPGroup, FPGroupElem, AutomorphismGroup, AutomorphismGroupElem

"""
TODO: document this
"""
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

"""
TODO: document this
"""
mutable struct PermGroupElem <: GroupElem
   parent::PermGroup
   X::GapObj   
end

"""
TODO: document this
"""
struct MatrixGroup <: Group
  X::GapObj
  function MatrixGroup(G::GapObj)
    @assert GAP.Globals.IsMatrixGroup(G)
    z = new(G)
    return z
  end
end

"""
TODO: document this
"""
mutable struct MatrixGroupElem <: GroupElem
   parent::MatrixGroup
   X::GapObj       
end

"""
TODO: document this
"""
struct PcGroup <: Group
  X::GapObj
  function PcGroup(G::GapObj)
    @assert GAP.Globals.IsPcGroup(G)
    z = new(G)
    return z
  end
end

"""
TODO: document this
"""
mutable struct PcGroupElem <: GroupElem
   parent::PcGroup
   X::GapObj
end

"""
TODO: document this
"""
struct FPGroup <: Group
  X::GapObj
  
  function FPGroup(G::GapObj)
    @assert GAP.Globals.IsFpGroup(G)
    z = new(G)
    return z
  end
end

"""
TODO: document this
"""
mutable struct FPGroupElem <: GroupElem
   parent::FPGroup
   X::GapObj
end

"""
TODO: document this
"""
mutable struct AutomorphismGroup{T} <: Group
  X::GapObj
  G::T
  function AutomorphismGroup{T}(G::GapObj, H::T) where T
    @assert GAP.Globals.IsAutomorphismGroup(G)
    z = new{T}(G, H)
    return z
  end
end

"""
TODO: document this
"""
mutable struct AutomorphismGroupElem{T} <: GroupElem
   parent::AutomorphismGroup{T}
   X::GapObj
end

#
# The array _gap_group_types contains pairs (X,Y) where
# X is a GAP filter such as IsPermGroup, and Y is a corresponding
# Julia type such as `PermGroup`.
#
const _gap_group_types = []

function _get_type(G::GapObj)
  for pair in _gap_group_types
    if pair[1](G)
      return pair[2]
    end
  end
  error("Not a known type of group")
end

function _get_gap_function(T)
  for pair in _gap_group_types
    if pair[2] == T
      return pair[1]
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
