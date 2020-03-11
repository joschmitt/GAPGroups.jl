function empty_record()
  return GAP.EvalString("rec()")
end

function oscar_to_gap(R::Generic.Ring)
  _set_comparison_number(R)
  RR = empty_record()
  ring_type = GAP.Globals.NewType(GAP.Globals.CollectionsFamily(GAP.Globals.aageneric_family), GAP.EvalString("IsAAGenericRing and IsAttributeStoringRep"))
  GAP.Globals.ObjectifyWithAttributes(RR, ring_type, GAP.Globals.ForeignPointer, R)
  return RR
end

function oscar_to_gap(R::Generic.Field)
  _set_comparison_number(R)
  # One can only set the characteristic of a field via its family, as far as I can see...
  field_family = GAP.Globals.NewFamily(GAP.julia_to_gap("Field_Family"), GAP.Globals.IsObject, GAP.Globals.IsAAGenericField)
  field_type = GAP.Globals.NewType(field_family, GAP.Globals.IsAAGenericField)
  GAP.Globals.SetCharacteristic(field_family, characteristic(R))
  RR = empty_record()
  GAP.Globals.ObjectifyWithAttributes(RR, field_type, GAP.Globals.ForeignPointer, R)
  return RR
end

function oscar_to_gap(x::Generic.RingElem)
  _set_comparison_number(x)
  xx = empty_record()
  ring_elem_type = GAP.Globals.NewType(GAP.Globals.aageneric_family, GAP.EvalString("IsAAGenericRingElem and IsAttributeStoringRep"))
  GAP.Globals.ObjectifyWithAttributes(xx, ring_elem_type, GAP.Globals.ForeignPointer, x)
  return xx
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

function oscar_to_gap(M::MatElem)
  _set_comparison_number(M)
  MM = empty_record()
  matrix_type = GAP.Globals.NewType(GAP.Globals.CollectionsFamily(GAP.Globals.CollectionsFamily(GAP.Globals.aageneric_family)), GAP.EvalString("IsAAGenericMatrixObj and IsAttributeStoringRep"))
  GAP.Globals.ObjectifyWithAttributes(MM, matrix_type, GAP.Globals.ForeignPointer, M, GAP.Globals.BaseDomain, oscar_to_gap(base_ring(M)))
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

function matrix_group(M::MatElem{T}...) where T
  return MatrixGroup(GAP.Globals.Group(( oscar_to_gap(MM) for MM in M )...))
end

function matrix_group(M::Vector{T}) where T <: MatElem
  return MatrixGroup(GAP.Globals.Group(GAP.julia_to_gap([ oscar_to_gap(MM) for MM in M ])))
end
