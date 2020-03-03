################################################################################
#
#  Group Homomorphism
#
################################################################################

mutable struct GAPGroupHomomorphism
   domain::GAPGroup
   codomain::GAPGroup
   image::GapObj
   
   function GAPGroupHomomorphism(domain::GAPGroup, codomain::GAPGroup, image::GapObj)
     z = new()
     z.domain = domain
     z.codomain = codomain
     z.image = image
     return z
   end
end

function Base.show(io::IO, x::GAPGroupHomomorphism)
  print(io, "Group homomorphism from
  \n $(domain(x)) \n to \n $(codomain(x))\n")
end

function hom(G::GAPGroup, H::GAPGroup, img::Function)
  #I create the gap function from the julia function
  #The julia function is supposed to be defined on GAPGroupElem
  #We need a function defined on the underlying GapObj
  function gap_fun(x::GapObj)
    el = GAPGroupElem(x, G)
    img_el = img(el)
    return img_el.X  
  end
  return GAPGroupHomomorphism(G, H, GAP.julia_to_gap(gap_fun))
end

function hom(G::GAPGroup, H::GAPGroup, gensG::Vector{GAPGroupElem}, imgs::Vector{GAPGroupElem})
  vgens = GAP.julia_to_gap([x.X for x in gensG])
  vimgs = GAP.julia_to_gap([x.X for x in imgs])
  mp = GAP.Globals.GroupHomomorphismByImages(G.X, H.X, vgens, vimgs)
  return GAPGroupHomomorphisms(G, H, mp)
end

function haspreimage(f::GAPGroupHomomorphism, x::GAPGroupElem)
  y = 
end

function domain(f::GAPGroupHomomorphism)
  return f.domain
end

function codomain(f::GAPGroupHomomorphism)
  return f.codomain
end

(f::GAPGroupHomomorphism)(x::GAPGroupElem) = image(f, x)

################################################################################
#
#  Image, Kernel, Cokernel
#
################################################################################

function image(f::GAPGroupHomomorphism, x::GAPGroupElem)
  return GAPGroupElem(f.image(x.X), codomain(f))
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

function _as_subgroup(H::GapObj, G::GAPGroup)
  function img(x::GAPGroupElem)
    return GAPGroupElem(x, H)
  end
  H1 = GAPGroup(H)
  return H1, hom(H1, G, img)
end

function sub(G::GAPGroup, elements::Vector{GAPGroupElem})
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

function index(G::GAPGroup, H::GAPGroup)
  i = index(G.X, H.X)
  return GAP.gap_to_julia(i)
end

###############################################################################
#
#  subgroups computation
#
###############################################################################

function normal_subgroups(G::GAPGroup)
  nsubs = GAP.gap_to_julia(GAP.Globals.NormalSubgroups(G.X))
  res = Vector{Tuple{GAPGroup, GAPGroupHomomorphism}}(undef, length(nsubs))
  for i = 1:length(res)
    N = nsubs[i]
    res[i] = _as_subgroup(N, G)
  end
  return res
end

function center(G::GAPGroup)
  Z = GAP.Globals.Center(G.X)
  return _as_subgroup(Z, G)
end

function centralizer(G::GAPGroup, H::GAPGroup)
  C = GAP.Globals.Centralizer(G.X, H.X)
  return _as_subgroup(C, G)
end

function centralizer(G::GAPGroup, x::GAPGroupElem)
  C = GAP.Globals.Centralizer(G.X, x.X)
  return _as_subgroup(C, G)
end


################################################################################
#
#  IsNormal, IsCharacteristic, IsSolvable, IsNilpotent
#
################################################################################

function isnormal(G::GAPGroup, H::GAPGroup)
  return GAP.Globals.IsNormal(G.X, H.X)
end

function ischaracteristic(G::GAPGroup, H::GAPGroup)
  return GAP.Globals.IsCharacteristicSubgroup(G.X, H.X)
end

function issolvable(G::GAPGroup)
  return GAP.Globals.IsSolvable(G.X)
end

function isnilpotent(G::GAPGroup)
  return GAP.Globals.IsNilpotent(G.X)
end

################################################################################
#
#  Quotient function
#
################################################################################

function quo(G::GAPGroup, elements::Vector{GAPGroupElem})
  elems_in_gap = GAP.julia_to_gap(GapObj[x.X for x in elements])
  H = GAP.Globals.Group(elems_in_gap)
  @assert GAP.Globals.IsNormal(G.X, H)
  mp = GAP.Globals.NaturalHomomorphismByNormalSubgroup(G.X, H)
  cod = GAPGroup(GAP.Globals.ImagesSource(mp))
  mp_julia = function(x::GAPGroupElem)
    el = mp(x.X)
    return GAPGroupElem(el, cod)
  end
  return hom(G, cod, mp_julia)
end


function quo(G::GAPGroup, H::GAPGroup)
  @assert GAP.Globals.IsNormal(G.X, H.X)
  mp = GAP.Globals.NaturalHomomorphismByNormalSubgroup(G.X, H.X)
  cod = GAPGroup(GAP.Globals.ImagesSource(mp))
  mp_julia = function(x::GAPGroupElem)
    el = mp(x.X)
    return GAPGroupElem(el, cod)
  end
  return hom(G, cod, mp_julia)
end
