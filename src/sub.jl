function Base.show(io::IO, x::GAPGroupHomomorphism)
  print(io, "Group homomorphism from
  \n $(domain(x)) \n to \n $(codomain(x))\n")
end

function hom(G::Group, H::Group, img::Function)
  #I create the gap function from the julia function
  #The julia function is supposed to be defined on GAPGroupElem
  #We need a function defined on the underlying GapObj
  function gap_fun(x::GapObj)
    el = group_element(G, x)
    img_el = img(el)
    return img_el.X  
  end
  return GAPGroupHomomorphism(G, H, GAP.julia_to_gap(gap_fun))
end

function hom(G::Group, H::Group, gensG::Vector, imgs::Vector)
  vgens = GAP.julia_to_gap([x.X for x in gensG])
  vimgs = GAP.julia_to_gap([x.X for x in imgs])
  mp = GAP.Globals.GroupHomomorphismByImages(G.X, H.X, vgens, vimgs)
  return GAPGroupHomomorphisms(G, H, mp)
end

function domain(f::GAPGroupHomomorphism)
  return f.domain
end

function codomain(f::GAPGroupHomomorphism)
  return f.codomain
end

(f::GAPGroupHomomorphism)(x::GroupElem) = image(f, x)

################################################################################
#
#  Image, Kernel, Cokernel
#
################################################################################

function image(f::GAPGroupHomomorphism, x::GroupElem)
  return group_element(codomain(f), f.image(x.X))
end

function kernel(f::GAPGroupHomomorphism)
  K = GAP.Globals.Kernel(f.image)
  return _as_subgroup(K, domain(f))
end

function image(f::GAPGroupHomomorphism)
  K = GAP.Globals.Image(f.image)
  return _as_subgroup(K, codomain(f))
end

function cokernel(f::GAPGroupHomomorphism)
  K, mK = image(f)
  return quo(codomain(f), K)
end

################################################################################
#
#  Subgroup function
#
################################################################################

function __as_subgroup(H::GapObj, G::T, ::Type{S}) where { T, S }
  function img(x::S)
    return group_element(G, x)
  end
  H1 = T(H)
  return H1, hom(H1, G, img)
end


function _as_subgroup(H::GapObj, G::T) where T <: Group
  S = elem_type(G)
  return __as_subgroup(H, G, S)
end

function sub(G::T, elements::Vector{S}) where T <: Group where S <: GroupElem
  @assert elem_type(G) == S
  elems_in_GAP = GAP.Globals.julia_to_gap(GapObj[x.X for x in elems])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(H, G)
end

###############################################################################
#
#  Index
#
###############################################################################

function index(G::Group, H::Group)
  i = index(G.X, H.X)
  return GAP.gap_to_julia(i)
end

###############################################################################
#
#  subgroups computation
#
###############################################################################

function normal_subgroups(G::Group)
  nsubs = GAP.gap_to_julia(GAP.Globals.NormalSubgroups(G.X))
  res = Vector{Tuple{GAPGroup, GAPGroupHomomorphism}}(undef, length(nsubs))
  for i = 1:length(res)
    N = nsubs[i]
    res[i] = _as_subgroup(N, G)
  end
  return res
end

function center(G::Group)
  Z = GAP.Globals.Center(G.X)
  return _as_subgroup(Z, G)
end

function centralizer(G::Group, H::Group)
  C = GAP.Globals.Centralizer(G.X, H.X)
  return _as_subgroup(C, G)
end

function centralizer(G::Group, x::GroupElem)
  C = GAP.Globals.Centralizer(G.X, x.X)
  return _as_subgroup(C, G)
end


################################################################################
#
#  IsNormal, IsCharacteristic, IsSolvable, IsNilpotent
#
################################################################################

function isnormal(G::Group, H::Group)
  return GAP.Globals.IsNormal(G.X, H.X)
end

function ischaracteristic(G::Group, H::Group)
  return GAP.Globals.IsCharacteristicSubgroup(G.X, H.X)
end

function issolvable(G::Group)
  return GAP.Globals.IsSolvable(G.X)
end

function isnilpotent(G::Group)
  return GAP.Globals.IsNilpotent(G.X)
end

################################################################################
#
#  Quotient function
#
################################################################################

function quo(G::T, elements::Vector{S}) where T <: Group where S <: GroupElem
  @assert elem_type(G) == S
  elems_in_gap = GAP.julia_to_gap(GapObj[x.X for x in elements])
  H = GAP.Globals.Group(elems_in_gap)
  @assert GAP.Globals.IsNormal(G.X, H)
  H1 = T(H)
  return quo(G, H)
end



function _get_type(G::GapObj)
  for i = 1:length(gap_group_types)
    if gap_group_types[i][1](G)
      return gap_group_types[i][2]
    end
  end
  error("Not a known type of group")
end

function quo(G::T, H::T) where T <: Group
  mp = GAP.Globals.NaturalHomomorphismByNormalSubgroup(G.X, H.X)
  cod = GAP.Globals.ImagesSource(mp)
  S = elem_type(G)
  S1 = _get_type(cod)
  codom = S1(cod)
  mp_julia = __create_fun(mp, codom, S)
  return codom, hom(G, codom, mp_julia)
end

function __create_fun(mp, codom, ::Type{S}) where S
  function mp_julia(x::S)
    el = GAP.Globals.Image(mp, x.X)
    return group_elem(codom, el)
  end
  return mp_julia
end
