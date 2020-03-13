function empty_record()
  return GAP.EvalString("rec()")
end

_gap_family(x::FinFieldElem) = GAP.Globals.aageneric_ffe_family
_gap_family(x::fmpz) = GAP.Globals.aageneric_cyclotomic_family
_gap_family(x::fmpq) = GAP.Globals.aageneric_cyclotomic_family

function _gap_family(x::nf_elem)
  if Nemo.iscyclo_type(parent(x))
    return GAP.Globals.aageneric_cyclotomic_family
  else
    return GAP.Globals.aageneric_algebraic_elements_family
  end
end

_gap_family(x::Generic.NCRingElem) = GAP.Globals.aageneric_family

function _cyclotomic_filters(x::nf_elem)
  filters = String[]
  if !Nemo.iscyclo_type(parent(x))
    return filters
  end
  push!(filters, "IsCyc")
  push!(filters, "IsCyclotomic")

  for i = 1:Hecke.degree(parent(x)) - 1
    if !iszero(coeff(x, i))
      return filters
    end
  end
  push!(filters, "IsRat")
  if coeff(x, 0) > 0
    push!(filters, "IsPosRat")
  elseif coeff(x, 0) < 0
    push!(filters, "IsNegRat")
  end
  if iszero(x)
    push!(filters, "IsZeroCyc")
  end
  if Hecke.isintegral(coeff(x, 0))
    push!(filters, "IsInt")
  end
  return filters
end

function _cyclotomic_filters(x::fmpq)
  filters = [ "IsCyc", "IsCyclotomic", "IsRat" ]
  if x > 0
    push!(filters, "IsPosRat")
  elseif x < 0
    push!(filters, "IsNegRat")
  else iszero(x)
    push!(filters, "IsZeroCyc")
  end
  if Hecke.isintegral(x)
    push!(filters, "IsInt")
  end
  return filters
end

function _cyclotomic_filters(x::fmpz)
  filters = [ "IsCyc", "IsCyclotomic", "IsRat", "IsInt" ]
  if x > 0
    push!(filters, "IsPosRat")
  elseif x < 0
    push!(filters, "IsNegRat")
  else iszero(x)
    push!(filters, "IsZeroCyc")
  end
  return filters
end

function _gap_type(x::Generic.NCRingElem)
  filters = String[ "IsAAGenericRingElem", "IsAttributeStoringRep", "IsAdditivelyCommutativeElement", "IsAssociativeElement" ]
  if typeof(x) <: Generic.RingElem
    # Objects of type Ring are by definition commutative in AbstractAlgebra, as far as I understand.
    push!(filters, "IsCommutativeElement")
  end

  if typeof(x) <: FieldElem
    push!(filters, "IsZDFRE")
    push!(filters, "IsMultiplicativeElementWithOne")
    push!(filters, "IsMultiplicativeElementWithInverse")
  end

  if typeof(x) <: FinFieldElem
    push!(filters, "IsFFE")
  end

  # This looks weird, but integers and rationals are cyclotomics for GAP
  if ( typeof(x) <: nf_elem && Nemo.iscyclo_type(parent(x)) ) || typeof(x) <: fmpq || typeof(x) <: fmpz
    append!(filters, _cyclotomic_filters(x))
  end

  if ( typeof(x) <: nf_elem && !Nemo.iscyclo_type(parent(x)) )
    push!(filters, "IsAlgebraicElement")
  end

  if typeof(x) <: PolyElem
    push!(filters, "IsPolynomialFunction")
    push!(filters, "IsPolynomialFunctionsFamilyElement")
  end

  filters_string = filters[1]
  for i = 2:length(filters)
    filters_string *= " and "
    filters_string *= filters[i]
  end

  return GAP.Globals.NewType(_gap_family(x), GAP.EvalString(filters_string))
end

function _gap_type(R::Generic.Ring)
  return GAP.Globals.NewType(GAP.Globals.CollectionsFamily(GAP.Globals.aageneric_family), GAP.EvalString("IsAAGenericRing and IsAttributeStoringRep"))
end

function _gap_type(F::Generic.Field)
  field_elem_family = GAP.Globals.NewFamily(GAP.julia_to_gap("Field_Family"), GAP.Globals.IsObject, GAP.Globals.IsAAGenericRingElem)
  GAP.Globals.SetCharacteristic(field_elem_family, characteristic(F))
  return GAP.Globals.NewType(GAP.Globals.CollectionsFamily(field_elem_family), GAP.EvalString("IsAAGenericField and IsAttributeStoringRep"))
end

function _gap_type(M::MatElem)
  family = _gap_family(M[1, 1])
  return GAP.Globals.NewType(GAP.Globals.CollectionsFamily(GAP.Globals.CollectionsFamily(family)), GAP.EvalString("IsAAGenericMatrixObj and IsAttributeStoringRep and IsConstantTimeAccessList and IsMutable"))
end

function oscar_to_gap(x::Generic.RingElem)
  _set_comparison_number(x)
  xx = empty_record()
  GAP.Globals.ObjectifyWithAttributes(xx, _gap_type(x), GAP.Globals.ForeignPointer, x)
  return xx
end

function oscar_to_gap(R::Generic.Ring)
  _set_comparison_number(R)
  RR = empty_record()
  GAP.Globals.ObjectifyWithAttributes(RR, _gap_type(R), GAP.Globals.ForeignPointer, R)
  return RR
end

function oscar_to_gap(R::Generic.Field)
  _set_comparison_number(R)
  RR = empty_record()
  GAP.Globals.ObjectifyWithAttributes(RR, _gap_type(R), GAP.Globals.ForeignPointer, R)
  return RR
end

function oscar_to_gap(M::MatElem)
  _set_comparison_number(M)
  MM = empty_record()
  GAP.Globals.ObjectifyWithAttributes(MM, _gap_type(M), GAP.Globals.ForeignPointer, M, GAP.Globals.BaseDomain, oscar_to_gap(base_ring(M)))
  return MM
end

function oscar_matrix_row_to_gap(M::MatElem, i::Int)
  v = view(M, i:i, :)
  _set_comparison_number(v)
  vv = empty_record()
  matrix_row_type = GAP.Globals.NewType(GAP.Globals.CollectionsFamily(GAP.Globals.aageneric_family), GAP.EvalString("IsAAGenericMatrixRowObj and IsAttributeStoringRep"))
  GAP.Globals.ObjectifyWithAttributes(vv, matrix_row_type, GAP.Globals.ForeignPointer, v, GAP.Globals.BaseDomain, oscar_to_gap(base_ring(M)))
  return vv
end

function gap_to_oscar(x::Main.ForeignGAP.MPtr)
  if GAP.Globals.HasForeignPointer(x)
    return GAP.Globals.ForeignPointer(x)
  elseif GAP.Globals.IsCollection(x)
    if length(x) != 0
      if GAP.Globals.HasForeignPointer(x[1])
        # Check whether this is a list of MatrixRows, that is a matrix
        if GAP.Globals.IsAAGenericMatrixRowObj(x[1])
          return vcat( GAP.Globals.ForeignPointer(x[i]) for i = 1:length(x) )
        else
          return typeof(GAP.Globals.ForeignPointer(x[1]))[ GAP.Globals.ForeignPointer(x[i]) for i = 1:length(x) ]
        end
      end
    end
  end
  error("Can't convert this object (yet?)")
end

################################################################################
#
#  Some functions we need on the GAP side
#
################################################################################

# GAP needs '<' for all types. We can't pass this through to AbstractAlgebra.
# For the time being (that is, until somebody has a better idea), we give every
# AbstractAlgebra-object a running number saved in a dictionary. So, an object
# R is "smaller" than an object S, if their running numbers are smaller.
global comparison_counter = BigInt(0)
global GAP_comparison_dictionary = Dict{Any, BigInt}()

function _set_comparison_number(x::Any)
  if haskey(GAP_comparison_dictionary, x)
    return nothing
  end

  global comparison_counter
  comparison_counter += 1
  GAP_comparison_dictionary[x] = comparison_counter
  return nothing
end

_get_comparison_number(x::Any) = GAP_comparison_dictionary[x]

function _print_to_string(x)
  io = IOBuffer()
  print(io, x)
  s = String(take!(io))
  close(io)
  return s
end

_conductor(x::fmpz) = 1
_conductor(x::fmpq) = 1

function _conductor(x::nf_elem)
  @assert Nemo.iscyclo_type(parent(x))
  n = _multiplicative_order(x)
  if iseven(n) && isodd(mod(n, 2))
    return div(n, 2)
  end
  return n
end

function _multiplicative_order(x::nf_elem)
  @assert Nemo.iscyclo_type(parent(x))
  n = Nemo.get_special(parent(x), :cyclo)
  if isodd(n)
    n = 2*n
  end

  y = one(parent(x))
  for i = 1:n
    y *= x
    if isone(y)
      return i
    end
  end

  return 0
end

function _multiplicative_order(x::Union{ fmpz, fmpq })
  if isone(x)
    return 1
  elseif isone(-x)
    return 2
  end
  return 0
end

_galois_cyc(x::Union{ fmpz, fmpq }) = deepcopy(x)

function _galois_cyc(x::nf_elem, k::Int)
  K = parent(x)
  a = one(K)
  y = K()
  for i = 0:Hecke.degree(K) - 1
    y += coeff(x, i)*a^k
    a = a*gen(K)
  end
  return y
end

# Convert integer cyclotomics to GAP integers without asking stupid questions
_cyclotomic_to_gapint(x::fmpz) = GAP.julia_to_gap(BigInt(x))
_cyclotomic_to_gapint(x::fmpq) = GAP.julia_to_gap(BigInt(numerator(x)))
_cyclotomic_to_gapint(x::nf_elem) = GAP.julia_to_gap(BigInt(numerator(coeff(x, 0))))

function _cyclotomic_to_gapint(M::MatElem{T}) where { T <: Union{ fmpz, fmpq, nf_elem } }
  N = Array{BigInt, 2}(undef, nrows(M), ncols(M))
  for i = 1:nrows(M)
    for j = 1:ncols(M)
      N[i, j] = _cyclotomic_to_gapint(M[i, j])
    end
  end
  return GAP.julia_to_gap(N, Val(true))
end

function matrix_group(M::MatElem{T}...) where T
  return MatrixGroup(GAP.Globals.Group(( oscar_to_gap(MM) for MM in M )...))
end

function matrix_group(M::Vector{T}) where T <: MatElem
  return MatrixGroup(GAP.Globals.Group(GAP.julia_to_gap([ oscar_to_gap(MM) for MM in M ])))
end
