# further possible functions: similar, iterate, literal_pow, parent_type
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
   parent::GAPGroup
   X::GapObj       
end

struct PolycyclicGroup <: Group
  X::GapObj
end

mutable struct PolycyclicGroupElem <: GroupElem
   parent::GAPGroup
   X::GapObj
end

struct FPGroup <: Group
  X::GapObj
end

mutable struct FPGroupElem <: GroupElem
   parent::GAPGroup
   X::GapObj
end

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

function symmetric_group(n::Int64)
   if n<1
     throw(ArgumentError("it must be a positive integer"))
   end
   return PermGroup(GAP.Globals.SymmetricGroup(n), n)
end

function alternating_group(n::Int64)
  if n<1
    throw(ArgumentError("it must be a positive integer"))
  end
  return PermGroup(GAP.Globals.AlternatingGroup(n), n)
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

function parent(x::GAPGroupElem)
  return x.parent
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

function order(::Type{T}, x::Union{GAPGroupElem, GAPGroup}) where T<:Number
   return T(order(x))
end

Base.:length(x::GAPGroup) = order(x)

function rand(x::GAPGroup)
   s = GAP.Globals.Random(x.X)
   return group_elem(x, s)
end

# one of the following should be non-parametric
rand_pseudo(G::GAPGroup) = rand(G)

# maxgroup (Sym(n), Sym(m)) = Sym( max(m,n)); it serves for operations between different groups
function maxgroup(x::GAPGroup, y::GAPGroup)
   if x.deg < y.deg return y end
   return x
end

Base.:*(x::GAPGroupElem, y::GAPGroupElem) = GAPGroupElem(x.X * y.X, maxgroup(x.par,y.par))

function ==(x::GAPGroup, y::GAPGroup)
   return x.X == y.X
end

function ==(x::GAPGroupElem, y::GAPGroupElem)
   return x.X == y.X
end

Base.:one(x::GAPGroup) = GAPGroupElem(GAP.Globals.Identity(x.X),x)
Base.:one(x::GAPGroupElem) = one(parent(x))
one!(x::GAPGroupElem) = one(parent(x))

Base.:show(io::IO, x::GAPGroupElem) =  print(GAP.gap_to_julia(GAP.Globals.StringView(x.X)))
Base.:show(io::IO, x::GAPGroup) = print(GAP.gap_to_julia(GAP.Globals.StringView(x.X)))

Base.:isone(x::GAPGroupElem) = x == one(parent(x))

Base.:inv(x::GAPGroupElem) = GAPGroupElem(GAP.Globals.Inverse(x.X),x.par)

inv!(out::GAPGroupElem, x::GAPGroupElem) = inv(x)  #if needed later

Base.:^(x::GAPGroupElem, y::Integer) = GAPGroupElem(x.X ^ y,x.par)

Base.:^(x::GAPGroupElem, y::GAPGroupElem) = GAPGroupElem(x.X ^ y.X, maxgroup(x.par, y.par))

Base.:<(x::GAPGroupElem, y::GAPGroupElem) = x.X < y.X

Base.:/(x::GAPGroupElem, y::GAPGroupElem) = x*y^-1

mul(x::GAPGroupElem, y::GAPGroupElem) = x*y
mul!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x*y

div_right(x::GAPGroupElem, y::GAPGroupElem) = x*inv(y)
div_left(x::GAPGroupElem, y::GAPGroupElem) = inv(y)*x
div_right!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x*inv(y)
div_left!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = inv(y)*x

conj(x::GAPGroupElem, y::GAPGroupElem) = x^y
conj!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x^y

comm(x::GAPGroupElem, y::GAPGroupElem) = x^-1*x^y
comm!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x^-1*x^y

function iterate(G::GAPGroup)
   L=GAP.Globals.List(G.X)
   len=length(L)
   iszero(len) && return nothing
   return GAPGroupElem(L[1],G.deg), (1,len)
end
function iterate(G::GAPGroup, state)
   L=GAP.Globals.List(G.X)
   s,len = state
   s==len && return nothing
   return GAPGroupElem(L[s+1],G.deg), (s+1,len)
end

function collect(G::GAPGroup)
   if G.deg>10  throw(ArgumentError("the group is too big")) end  # to be discussed whether the group is too big or not
   L=GAP.gap_to_julia(GAP.Globals.List(G.X))
   return [GAPGroupElem(x,G) for x in L]
end   

#maybe in future add more checks
hasorder(x::GAPGroup) = true
hasorder(x::GAPGroupElem) = true
hasgens(x::GAPGroup) = true

function perm(L::Array{Int64,1})
   return GAPGroupElem(GAP.Globals.PermList(GAP.julia_to_gap(L)),symmetric_group(length(L)))
end

function perm(g::GAPGroup, L::Array{Int64,1})
   x=GAP.Globals.PermList(GAP.julia_to_gap(L))
   if GAP.Globals.IN(x,g.X) return GAPGroupElem(GAP.Globals.PermList(GAP.julia_to_gap(L)),g)
   else  throw(ArgumentError("the element does not embed in the group"))
   end
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of arrays (not necessarly disjoint)
function cperm(L::Union{Array{Int64,1},UnitRange{Int64}}...)
   if length(L)==0
      return one(symmetric_group(1))
   else
      x=GAP.Globals.Product(GAP.julia_to_gap([GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y))) for y in L]))
      return GAPGroupElem(x, symmetric_group(maximum([maximum(y) for y in L])))
#      return prod([GAPGroupElem(GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y))),symmetric_group(maximum(y))) for y in L])
   end
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of arrays (not necessarly disjoint)
# WARNING: we allow e.g. PermList([2,3,1,4,5,6]) in Sym(3)
function cperm(g::GAPGroup,L::Union{Array{Int64,1},UnitRange{Int64}}...)
   if length(L)==0
      return one(g)
   else
      x=GAP.Globals.Product(GAP.julia_to_gap([GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y))) for y in L]))
      if GAP.Globals.IN(x,g.X) return GAPGroupElem(x,g)
      else throw(ArgumentError("the element does not embed in the group"))
      end
   end
end

function listperm(x::GAPGroupElem)
   return [x(i) for i in 1:x.par.deg]
end

function gens(G::GAPGroup)
   L=GAP.Globals.GeneratorsOfGroup(G.X)
   l=length(L)
   return [GAPGroupElem(L[i],G) for i in 1:l]
end

function gens(G::GAPGroup, i::Integer)
   L=GAP.Globals.GeneratorsOfGroup(G.X)
   return GAPGroupElem(L[i],G)
end

ngens(G::GAPGroup) = length(gens(G))

Base.:sign(x::GAPGroupElem) = GAP.Globals.SignPerm(x.X)

Base.:isless(x::GAPGroupElem, y::GAPGroupElem) = x<y

#embedding of a permutation in permutation group
function (G::GAPGroup)(x::GAPGroupElem)
   if !GAP.Globals.IN(x.X,G.X)
      throw(ArgumentError("the element does not embed in the group"))
   end
   x1=deepcopy(x)
   x1.par=G
   return x1
end

#evaluation function
function (x::GAPGroupElem)(n)
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

