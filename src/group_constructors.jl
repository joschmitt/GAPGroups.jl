################################################################################
#
#  Some basic constructors
#  
################################################################################

export symmetric_group, alternating_group, small_group, small_groups_id, transitive_group,
       cyclic_group, abelian_group, mathieu_group, free_abelian_group, dihedral_group,
       quaternion_group, GL, SL, isquaternion_group


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

function issymmetric_group(G::Group)
  return GAP.Globals.IsSymmetricGroup(G.X)
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

function isalternating_group(G::Group)
  return GAP.Globals.IsAlternatingGroup(G.X)
end

function small_group(n::Int, m::Int)
  return PcGroup(GAP.Globals.SmallGroup(n, m))
end

function small_groups_id(G::Group)
  r = GAP.Globals.IdGroup(G.X)
  res = GAP.gap_to_julia(GAP.gap_to_julia(r))
  return (res[1], res[2])
end

function transitive_group(n::Int, m::Int)
  return PermGroup(GAP.Globals.TransitiveGroup(n, m))
end

function cyclic_group(n::Int)
  return PcGroup(GAP.Globals.CyclicGroup(n))
end

function cyclic_group(::Type{T}, n::Int) where T <: Group
  return T(GAP.Globals.CyclicGroup(_get_gap_function(T), n))
end

function iscyclic(G::Group)
  return GAP.Globals.IsCyclic(G.X)
end

function abelian_group(v::Vector{Int})
  for i = 1:length(v)
    iszero(v[i]) && error("Cannot represent an infinite group as a polycyclic group")
  end
  v1 = GAP.julia_to_gap(v)
  return PcGroup(GAP.Globals.AbelianGroup(v1))
end

function abelian_group(::Type{T}, v::Vector{Int}) where T <: Group
  v1 = GAP.julia_to_gap(v)
  return T(GAP.Globals.AbelianGroup(_get_gap_function(T), v1))
end

function isabelian(G::Group)
  return GAP.Globals.IsAbelian(G.X)
end

function mathieu_group(n::Int)
  @assert n in Int[9, 10, 11, 12, 21, 22, 23, 24]
  return PermGroup(GAP.Globals.MathieuGroup(n), n)
end


################################################################################
#
# begin FpGroups
#
################################################################################

function free_group(n::Int)
   return FPGroup(GAP.Globals.FreeGroup(n))
end

function free_group(L::String...)
   J=GAP.julia_to_gap([GAP.julia_to_gap(x) for x in L])
   return FPGroup(GAP.Globals.FreeGroup(J))
end

function free_abelian_group(n::Int)
  return FPGroup(GAP.Globals.FreeAbelianGroup(n))
end

# for the definition of group modulo relations, see the quo function in the sub.jl section

function free_group(G::FPGroup)
   return FPGroup(GAP.Globals.FreeGroupOfFpGroup(G.X))
end

################################################################################
#
# end FpGroups
#
################################################################################


function dihedral_group(n::Int)
  @assert iseven(n)
  return PcGroup(GAP.Globals.DihedralGroup(n))
end

function dihedral_group(::Type{T}, n::Int) where T <: Group
  @assert iseven(n)
  return T(GAP.Globals.DihedralGroup(_get_gap_function(T), n))
end

function isdihedral_group(G::Group)
  return GAP.Globals.IsDihedralGroup(G.X)
end

function quaternion_group(n::Int)
  @assert divisible(n, 4)
  return PcGroup(GAP.Globals.QuaternionGroup(n))
end

function quaternion_group(::Type{T}, n::Int) where T <: Group 
  @assert divisible(n, 4)
  return T(GAP.Globals.QuaternionGroup(_get_gap_function(T)), n)
end

function isquaternion_group(G::Group)
  return GAP.Globals.IsQuaternionGroup(G.X)
end

################################################################################
#
# start isometry groups
#
################################################################################

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

function symplectic_group(n::Int, q::Int)
  return MatrixGroup(GAP.Globals.Sp(n, q))
end

function symplectic_group(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.Sp(_get_gap_function(T), n, q))
end

function unitary_group(n::Int, q::Int)
  return MatrixGroup(GAP.Globals.GU(n, q))
end

function unitary_group(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.GU(_get_gap_function(T), n, q))
end

function special_unitary(n::Int, q::Int)
  return MatrixGroup(GAP.Globals.SU(n, q))
end

function special_unitary(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.SU(_get_gap_function(T), n, q))
end

function orthogonal_group(n::Int, q::Int)
  return MatrixGroup(GAP.Globals.GO(n, q))
end

function orthogonal_group(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.GO(_get_gap_function(T), n, q))
end

function orthogonal_group(e::Int, n::Int, q::Int)
  return MatrixGroup(GAP.Globals.GO(e, n, q))
end

function orthogonal_group(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.GO(_get_gap_function(T), e, n, q))
end

function special_orthogonal_group(n::Int, q::Int)
  return MatrixGroup(GAP.Globals.SO(n, q))
end

function special_orthogonal_group(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.SO(_get_gap_function(T), n, q))
end

function special_orthogonal_group(e::Int, n::Int, q::Int)
  return MatrixGroup(GAP.Globals.SO(e, n, q))
end

function special_orthogonal_group(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.SO(_get_gap_function(T), e, n, q))
end

function omega_group(n::Int, q::Int)
  return MatrixGroup(GAP.Globals.Omega(n, q))
end

function omega_group(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.Omega(_get_gap_function(T), n, q))
end

function omega_group(e::Int, n::Int, q::Int)
  return MatrixGroup(GAP.Globals.Omega(e, n, q))
end

function omega_group(::Type{T}, n::Int, q::Int) where T <: Group
  return T(GAP.Globals.Omega(_get_gap_function(T), e, n, q))
end

################################################################################
#
# end isometry groups
#
################################################################################
