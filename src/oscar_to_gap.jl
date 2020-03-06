function empty_record()
  return GAP.EvalString("rec()")
end

function oscar_to_gap(R::Generic.Ring)
  _set_comparison_number(R)
  RR = empty_record()
  GAP.Globals.ObjectifyWithAttributes(RR, GAP.Globals.aageneric_ring_type, GAP.Globals.JuliaPointer, R)
  return RR
end

function oscar_to_gap(R::Generic.Field)
  _set_comparison_number(R)
  # One can only set the characteristic of a field via its family, as far as I can see...
  field_family = GAP.Globals.NewFamily(GAP.julia_to_gap("Field_Family"), GAP.Globals.IsObject, GAP.Globals.IsAAGenericField)
  field_type = GAP.Globals.NewType(field_family, GAP.Globals.IsAAGenericField)
  GAP.Globals.SetCharacteristic(field_family, characteristic(R))
  RR = empty_record()
  GAP.Globals.ObjectifyWithAttributes(RR, field_type, GAP.Globals.JuliaPointer, R)
  return RR
end

function gap_to_oscar(x::Main.ForeignGAP.MPtr)
  if GAP.Globals.HasJuliaPointer(x)
    return GAP.Globals.JuliaPointer(x)
  elseif GAP.Globals.IsCollection(x)
    if length(x) != 0
      if GAP.Globals.HasJuliaPointer(x[1])
        # Check whether this is a list of MatrixRows, that is a matrix
        if GAP.Globals.IsAAGenericMatrixRowObj(x[1])
          return vcat( GAP.Globals.JuliaPointer(x[i]) for i = 1:length(x) )
        else
          return typeof(GAP.Globals.JuliaPointer(x[1]))[ GAP.Globals.JuliaPointer(x[i]) for i = 1:length(x) ]
        end
      end
    end
  end
  error("Can't convert this object (yet?)")
end

function oscar_to_gap(M::MatElem)
  _set_comparison_number(M)
  MM = empty_record()
  GAP.Globals.ObjectifyWithAttributes(MM, GAP.Globals.aageneric_matrix_type, GAP.Globals.JuliaPointer, M, GAP.Globals.BaseDomain, oscar_to_gap(base_ring(M)))
  return MM
end

function oscar_matrix_row_to_gap(M::MatElem, i::Int)
  v = view(M, i:i, :)
  _set_comparison_number(v)
  vv = empty_record()
  GAP.Globals.ObjectifyWithAttributes(vv, GAP.Globals.aageneric_matrix_row_type, GAP.Globals.JuliaPointer, v, GAP.Globals.BaseDomain, oscar_to_gap(base_ring(M)))
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
