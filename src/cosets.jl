# T=type of the group, S=type of the element
mutable struct GroupCoset{T<: Group, S <: GroupElem} 
   X::T
   H::T
   repr::S
   side::Symbol     # says if the coset is left or right
   coset::GapObj
end

function _group_coset(X::Group, H::Group, repr::GroupElem, side::Symbol, coset::GapObj)
  return GroupCoset{typeof(X), typeof(repr)}(X, H, repr, side, coset)
end

function right_coset(H::Group, g::GroupElem)
   @assert elem_type(H) == typeof(g)
   if !GAP.Globals.IsSubset(parent(g).X, H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return _group_coset(parent(g), H, g, :right, GAP.Globals.RightCoset(H.X,g.X))
end

function left_coset(H::Group, g::GroupElem)
   @assert elem_type(H) == typeof(g)
   if !GAP.Globals.IsSubset(parent(g).X, H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return _group_coset(parent(g), H, g, :left, GAP.Globals.RightCoset(GAP.Globals.ConjugateSubgroup(H.X,GAP.Globals.Inverse(g.X)),g.X))
end

function show(io::IO, x::GroupCoset)
   if x.right == "right"
      print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.H.X))*" * "*GAP.gap_to_julia(GAP.Globals.StringView(x.repr.X)))
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.repr.X))*" * "*GAP.gap_to_julia(GAP.Globals.StringView(x.H.X)))
   end
   return nothing
end

acting_domain(C::GroupCoset) = C.H

representative(C::GroupCoset) = C.repr

function elements(C::GroupCoset)
  L = GAP.Globals.AsList(C.coset)
  l = Vector{elem_type(C.X)}(undef, length(L))
  for i = 1:length(l)
    l[i] = group_element(C.X, L[i])
  end
  return l
end

order(C::GroupCoset) = GAP.Globals.Size(C.coset)

is_bicoset(C::GroupCoset) = GAP.Globals.IsBiCoset(C.coset)

function right_cosets(G::Group, H::Group)
  L = GAP.Globals.RightCosets(G.X, H.X)
  l = Vector{GroupCoset{typeof(G), elem_type(G)}}(undef, length(L))
  for i = 1:length(l)
    l[i] = _group_coset(G, H, group_element(G, GAP.Globals.Representative(L[i])), :right, L[i])
  end
  return l
end


function right_transversal(G::T, H::T) where T<: Group
   L = GAP.Globals.RightTransversal(G.X,H.X)
   l = Vector{elem_type(G)}(undef, length(L))
   return elem_type(G)[group_element(G, x) for x in GAP.gap_to_julia(L)]
end


# T=type of the group, S=type of the element
mutable struct GroupDoubleCoset{T <: Group, S <: GroupElem}
   X::T
   G::T
   H::T
   repr::S
   coset::GapObj
end

Base.:show(io::IO, x::GroupDoubleCoset) =  print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.G.X))*" * "*GAP.gap_to_julia(GAP.Globals.StringView(x.repr.X))*" * "*GAP.gap_to_julia(GAP.Globals.StringView(x.H.X)))

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
