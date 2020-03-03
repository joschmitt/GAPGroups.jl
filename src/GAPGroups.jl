# further possible functions: similar, iterate, literal_pow, parent_type


module GAPGroups 

using GAP

import Base.==
import Base.rand
import Base.conj
import Base.conj!
import Base.parent
import Base.eltype
import Base.iterate
import Base.collect

export symmetric_group, order, perm, cperm, hasorder, hasgens, gens, ngens, comm, comm!, inv!, rand_pseudo, one!, div_right, div_left, div_right!, div_left!, elem_type, deg, mul, mul!, listperm
     #conj!, conj

include("./types.jl")

elem_type(::PermGroup) = PermGroupElem
elem_type(::MatrixGroup) = MatrixGroupElem
elem_type(::FPGroup) = FPGroupElem
elem_type(::PolycyclicGroup) = PolycyclicGroupElem


function group_element(G::PermGroup, x::GapObj)
  return PermGroupElem(G, x)
end

function group_element(G::MatrixGroup, x::GapObj)
  return MatrixGroupElem(G, x)
end

################################################################################
#
#  Some basic constructors
#  
################################################################################

function symmetric_group(n::Int64)
  if n < 1
    throw(ArgumentError("it must be a positive integer"))
  end
  return PermGroup(GAP.Globals.SymmetricGroup(n), n)
end

function symmetric_group(::Type{T}, n::Int) where T <: Group
  if n < 1
    throw(ArgumentError("n must be a positive integer"))
  end
  return T(GAP.Globals.SymmetricGroup(_get_gap_function(T), n))
end

function alternating_group(n::Int64)
  if n < 1
    throw(ArgumentError("it must be a positive integer"))
  end
  return PermGroup(GAP.Globals.AlternatingGroup(n), n)
end

function alternating_group(::Type{T}, n::Int) where T <: Group
  if n < 1
    throw(ArgumentError("it must be a positive integer"))
  end
  return T(GAP.Globals.AlternatingGroup(_get_gap_function(T), n))
end

function small_group(n::Int, m::Int)
  return PolycyclicGroup(GAP.Globals.SmallGroup(n, m))
end

function cyclic_group(n::Int)
  return PolycyclicGroup(GAP.Globals.CyclicGroup(n))
end

function cyclic_group(::Type{T}, n::Int) where T <: Group
  return T(GAP.Globals.CyclicGroup(_get_gap_function(T), n))
end

function abelian_group(v::Vector{Int})
  for i = 1:length(v)
    iszero(v[i]) && error("Cannot represent an infinite group as a polycyclic group")
  end
  v1 = GAP.julia_to_gap(v)
  return PolycyclicGroup(GAP.Globals.AbelianGroup(v1))
end

function abelian_group(::Type{T}, v::Vector{Int]) where T <: Group
  v1 = GAP.julia_to_gap(v)
  return T(GAP.Globals.AbelianGroup(_get_gap_function(T), v1))
end

function mathieu_group(n::Int)
  @assert n in Int[9, 10, 11, 12, 21, 22, 23, 24]
  return PermGroup(GAP.Globals.MathieuGroup(n), n)
end

function free_abelian_group(n::Int)
  return FPGroup(GAP.Globals.FreeAbelianGroup(n))
end

function dihedral_group(n::Int)
  @assert iseven(n)
  return PolycyclicGroup(GAP.Globals.DihedralGroup(n))
end

function dihedral_group(::Type{T}, n::Int) where T <: Group
  @assert iseven(n)
  return T(GAP.Globals.DihedralGroup(_get_gap_function(T), n))
end

function quaternion_group(n::Int)
  return PolycyclicGroup(GAP.Globals.QuaternionGroup(n))
end

function quaternion_group(::Type{T}, n::Int) where T <: Group 
  @assert divisible(n, 4)
  return T(GAP.Globals.QuaternionGroup(_get_gap_function(T)), n)
end


function GL(n::Int, q::Int)
  return MatrixGroup(GAP.Globals.GL(n, q))
end

function GL(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.GL(_get_gap_function(T),n, q))
end

function SL(n::Int, q::Int)
  return MatrixGroup(GAP.Globals.SL(n, q))
end

function SL(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.SL(_get_gap_function(T), n, q))
end




function elements(G::T) where T <: Group
  els = GAP.gap_to_julia(GAP.Globals.Elements(G.X))
  elems = Vector{elem_type(G)}(undef, length(els))
  i = 1
  for x in els
    elems[i] = group_element(G, x)
    i += 1
  end
  return elems
end

# to be fixed later
Base.:hash(x::GroupElem) = 0

#shorter writing
#Base.:parent(x::GAPGroupElem) = x==one(symmetric_group(2)) ? symmetric_group(1) : symmetric_group(GAP.Globals.LargestMovedPointPerm(x.X))

function parent(x::GroupElem)
  return x.parent
end

function isfinite(G::PermGroup)
  return true
end

function isfinite(G::Group)
  return GAP.Globals.IsFinite(G)
end

function degree(x::PermGroup)
   return x.deg
end

function order(x::Group)
   return GAP.gap_to_julia(GAP.Globals.Size(x.X))
end

function order(x::GroupElem)
   return GAP.gap_to_julia(GAP.Globals.Order(x.X))
end

function order(::Type{T}, x::Union{GroupElem, Group}) where T<:Number
   return T(order(x))
end

#Base.:length(x::Group) = order(x)

function rand(x::Group)
   s = GAP.Globals.Random(x.X)
   return group_elem(x, s)
end

# one of the following should be non-parametric
rand_pseudo(G::Group) = rand(G)


function _maxgroup(x::T, y::T) where T <: Group
   if x == y
     return x
   elseif GAP.Globals.IsSubset(x, y)
     return x
   elseif GAP.Globals.IsSubset(y, x)
     return y
   else
     error("Not yet implemented")
   end
end

#We need a lattice of groups to implement this properly
function _prod(x::GroupElem, y::GroupElem)
  G = _maxgroup(parent(x), parent(y))
  return group_element(G, x.X*y.X)
end

Base.:*(x::GroupElem, y::GroupElem) = _prod(x, y)

function ==(x::Group, y::Group)
   return x.X == y.X
end

function ==(x::GroupElem, y::GroupElem)
   return x.X == y.X
end

Base.:one(x::Group) = group_element(x, GAP.Globals.Identity(x.X))
Base.:one(x::GroupElem) = one(parent(x))
one!(x::GroupElem) = one(parent(x))

Base.:show(io::IO, x::GroupElem) =  print(GAP.gap_to_julia(GAP.Globals.StringView(x.X)))
Base.:show(io::IO, x::Group) = print(GAP.gap_to_julia(GAP.Globals.StringView(x.X)))

Base.:isone(x::GroupElem) = x == one(parent(x))

Base.:inv(x::GroupElem) = group_element(parent(x), GAP.Globals.Inverse(x.X))

inv!(out::GroupElem, x::GroupElem) = inv(x)  #if needed later

Base.:^(x::GroupElem, y::Integer) = group_element(parent(x), x.X ^ y)

Base.:^(x::GroupElem, y::GroupElem) = group_element(_maxgroup(parent(x), parent(y), x.X ^ y.X))

#Is this useful?
#Base.:<(x::GroupElem, y::GroupElem) = x.X < y.X

Base.:/(x::GroupElem, y::GroupElem) = x*y^-1

mul(x::GroupElem, y::GroupElem) = x*y
mul!(out::GroupElem, x::GroupElem, y::GroupElem) = x*y

div_right(x::GroupElem, y::GroupElem) = x*inv(y)
div_left(x::GroupElem, y::GroupElem) = inv(y)*x
div_right!(out::GroupElem, x::GroupElem, y::GroupElem) = x*inv(y)
div_left!(out::GroupElem, x::GroupElem, y::GroupElem) = inv(y)*x

conj(x::GroupElem, y::GroupElem) = x^y
conj!(out::GroupElem, x::GroupElem, y::GroupElem) = x^y

comm(x::GroupElem, y::GroupElem) = x^-1*x^y
comm!(out::GroupElem, x::GroupElem, y::GroupElem) = x^-1*x^y

function iterate(G::Group)
   L=GAP.Globals.List(G.X)
   len=length(L)
   iszero(len) && return nothing
   return group_elem(G, L[1]), (1,len)
end

function iterate(G::Group, state)
   L = GAP.Globals.List(G.X)
   s,len = state
   s == len && return nothing
   return group_elem(G, L[s+1]), (s+1,len)
end

function collect(G::Group)
  L = GAP.gap_to_julia(GAP.Globals.List(G.X))
  return elem_type(G)[group_elem(G, x) for x in L]
end   

function perm(L::Array{Int64,1})
   return PermGroupElem(symmetric_group(length(L)), GAP.Globals.PermList(GAP.julia_to_gap(L)))
end

function perm(g::PermGroup, L::Array{Int64,1})
   x = GAP.Globals.PermList(GAP.julia_to_gap(L))
   if GAP.Globals.IN(x,g.X) 
     return PermGroupElem(g, GAP.Globals.PermList(GAP.julia_to_gap(L)))
   end
   throw(ArgumentError("the element does not embed in the group"))
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of arrays (not necessarly disjoint)
function cperm(L::Union{Array{Int64,1},UnitRange{Int64}}...)
   if length(L)==0
      return one(symmetric_group(1))
   else
      return prod([PermGroupElem(symmetric_group(maximum(y)), GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y)))) for y in L])
   end
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of arrays (not necessarly disjoint)
# WARNING: we allow e.g. PermList([2,3,1,4,5,6]) in Sym(3)
function cperm(g::PermGroup,L::Union{Array{Int64,1},UnitRange{Int64}}...)
   if length(L)==0
      return one(g)
   else
      x=GAP.Globals.Product(GAP.julia_to_gap([GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y))) for y in L]))
      if GAP.Globals.IN(x,g.X) return PermGroupElem(g, x)
      else throw(ArgumentError("the element does not embed in the group"))
      end
   end
end

function listperm(x::PermGroupElem)
   return [x(i) for i in 1:x.par.deg]
end

function gens(G::Group)
   L = GAP.Globals.GeneratorsOfGroup(G.X)
   res = Vector{elem_type(G)}(undef, length(L))
   for i = 1:length(res)
     res[i] = group_element(G, L[i])
   end
   return res
end

function gens(G::Group, i::Integer)
   L = GAP.Globals.GeneratorsOfGroup(G.X)
   return group_element(G, L[i])
end

ngens(G::Group) = length(gens(G))

Base.:sign(x::PermGroupElem) = GAP.Globals.SignPerm(x.X)

#Base.:isless(x::GAPGroupElem, y::GAPGroupElem) = x<y

#embedding of a permutation in permutation group
function (G::PermGroup)(x::PermGroupElem)
   if !GAP.Globals.IN(x.X,G.X)
      throw(ArgumentError("the element does not embed in the group"))
   end
   x1=deepcopy(x)
   x1.par=G
   return x1
end

#evaluation function
function (x::PermGroupElem)(n)
   return GAP.Globals.OnPoints(n,x.X)
end

include("./sub.jl")

function conjugate_subgroup(G::T, x::GroupElem) where T<:Group
  return T(GAP.Globals.ConjugateSubgroup(G.X,x.X))
end

################################################################################
#
# cosets and double cosets
#
################################################################################

# T=type of the group, S=type of the element
mutable struct GroupCoset{T,S} where {T<:Group, S<:GroupElem}
   X::T
   H::T
   repr::S
   right::String     # says if the coset is left or right
   coset::GapObj
end

function right_coset(H::Group,g::GroupElem)
   if !GAP.Globals.IsSubset(parent(g).X,H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return GroupCoset(parent(g),H,g,"right",GAP.Globals.RightCoset(H.X,g.X))
end

acting_domain(C::GroupCoset) = C.H

representative(C::GroupCoset) = C.repr

function elements(C::GroupCoset)
   L=GAP.gap_to_julia(GAP.Globals.AsList(C.coset))
   return elem_type(C.X)[group_element(C.X,x) for x in L]
end

order(C::GroupCoset) = GAP.Globals.Size(C.coset)

is_bicoset(C::GroupCoset) = GAP.Globals.IsBiCoset(C.coset)

function coset_decomposition(G::T, H::T) where T<:Group
   L=GAP.Globals.CosetDecomposition(G.X,H.X)
   return T[[group_element(G.X,x) for x in GAP.gap_to_julia(J)] for J in GAP.gap_to_julia(L)]
end

function right_transversal(G::T, H::T) where T<:Group
   L=GAP.Globals.RightTransversal(G.X,H.X)
   return T[group_elem(G.X,x) for x in GAP.gap_to_julia(L)]
end

# T=type of the group, S=type of the element
mutable struct GroupDoubleCoset{T,S} where {T<:Group, S<:GroupElem}
   X::T
   G::T
   H::T
   repr::S
   coset::GapObj
end

function double_coset(G::Group, g::GroupElem, H::Group)
   if !GAP.Globals.IsSubset(parent(g).X,G.X)
      throw(ArgumentError("G is not a subgroup of parent(g)"))
   end
   if !GAP.Globals.IsSubset(parent(g).X,H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return GroupDoubleCoset(parent(g),G,H,g,GAP.Globals.DoubleCoset(G,g,H))
end

function elements(C::GroupDoubleCoset)
   L=GAP.gap_to_julia(GAP.Globals.AsList(C.coset))
   return elem_type(C.X)[group_element(C.X,x) for x in L]
end

order(C::Union{GroupCoset,GroupDoubleCoset}) = GAP.Globals.Size(C.coset)

rand(C::Union{GroupCoset,GroupDoubleCoset}) = group_element(C.X, GAP.Globals.Random(C.coset))

################################################################################
#
#   Conjugacy Classes
#
################################################################################

struct GroupConjClass{T,S} where {T<:Group, S<:GroupElem}
   X::T
   repr::S
   CC::GapObj
end

conjugacy_class(G::Group, g::GroupElem) = GroupConjClass(G,g,GAP.Globals.ConjugacyClass(G.X,g.X))

order(C::GroupConjClass) = GAP.Globals.Size(C.CC)

rand(C::GroupConjClass) = group_elem(C.X,GAP.Globals.Random(C.CC))

function elements(C::GroupConjClass)
   L=GAP.gap_to_julia(GAP.Globals.AsList(C.CC))
   return elem_type(C.X)[group_element(C.X,x) for x in L]
end

representative(C::GroupConjClass) = C.repr

function conjugacy_classes(G::Group)
   L=GAP.gap_to_julia(GAP.Globals.ConjugacyClasses(G.X))
   return [GroupConjClass(G,group_element(G,GAP.Globals.Representative(cc)),cc) for cc in L]
end

rand(C::GroupConjClass) = group_element(C.X, GAP.Globals.Random(C.CC))

nr_conjugacy_classes(G::Group) = GAP.Globals.NrConjugacyClasses(G.X)

Base.:^(x::GroupElem, y::Integer) = group_element(x.parent, x.X ^ y)

Base.:^(x::GroupElem, y::GroupElem) = group_element(x.parent, x.X ^ y.X)

Base.:^(H::Group, y::GroupElem) = typeof(H)(H.X ^ y.X)

# for elements
function is_conjugate(G::Group, x::GroupElem, y::GroupElem)
   if GAP.Globals.IsConjugate(G.X, x.X, y.X)
      return true, group_element(G,GAP.Globals.RepresentativeAction(G.X, x.X, y.X))
   else
      return false, nothing
   end
end

# for subgroups
function is_conjugate(G::Group, H::Group, K::Group)
   if GAP.Globals.IsConjugate(G.X, H.X, K.X)
      return true, typeof(G)(GAP.Globals.RepresentativeAction(G.X, x.X, y.X))
   else
      return false, nothing
   end
end

