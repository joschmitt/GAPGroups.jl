# AAGeneric as in AbstractAlgebra Generic
DeclareCategory("IsAAGenericObject", IsObject);
DeclareSynonym("IsAAGenericRing", IsAAGenericObject and IsRing);
DeclareSynonym("IsAAGenericField", IsAAGenericObject and IsField);
DeclareSynonym("IsAAGenericMatrixObj", IsAAGenericObject and IsMatrixObj and IsMatrix);
DeclareSynonym("IsAAGenericMatrixRowObj", IsAAGenericObject and IsVectorObj);

aageneric_family := NewFamily( "AAGeneric_Family", IsObject, IsAAGenericObject );
aageneric_ring_type := NewType(aageneric_family, IsAAGenericRing);

#aageneric_field_family := NewFamily( "AAGeneric_Field_Family", IsObject, IsAAGenericField );
#aageneric_field_type := NewType(aageneric_field_family, IsAAGenericField);

coll_aageneric_family := CollectionsFamily(aageneric_family);
aageneric_matrix_row_type := NewType(coll_aageneric_family, IsAAGenericMatrixRowObj and IsAttributeStoringRep);
aageneric_matrix_type := NewType(CollectionsFamily(coll_aageneric_family), IsAAGenericMatrixObj and IsAttributeStoringRep);

# The general comparing function using the global GAP_comparison_dictionary
InstallMethod( \<,
  ["IsAAGenericObject", "IsAAGenericObject"],
  function( x, y )
    return Julia.Base.\<(Julia._get_comparison_number(JuliaPointer(x)), Julia._get_comparison_number(JuliaPointer(y)));
  end
);

################################################################################
#
#  Basic functions
#
################################################################################

InstallMethod( String,
  ["IsAAGenericObject"],
  x -> String( JuliaPointer( x ) )
);

InstallMethod( PrintObj,
  ["IsAAGenericObject"],
  x -> Print( PrintString( x ) )
);

InstallMethod( ViewString,
  ["IsAAGenericObject"],
  x -> ViewString( JuliaPointer(x) )
);

################################################################################
#
#  Basic functions for matrices
#
################################################################################

InstallMethod( NumberRows,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Generic.nrows(JuliaPointer(x))
);

InstallMethod( Length,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Nemo.nrows(JuliaPointer(x))
);

InstallMethod( NumberColumns,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Generic.ncols(JuliaPointer(x))
);

InstallMethod( RankMat,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Generic.rank(JuliaPointer(x))
);

InstallMethod( \=,
  ["IsAAGenericMatrixObj", "IsAAGenericMatrixObj"],
  function( x, y )
    return Julia.Base.\=\=(JuliaPointer(x), JuliaPointer(y));
  end
);

InstallMethod( \+,
  ["IsAAGenericMatrixObj", "IsAAGenericMatrixObj"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\+(JuliaPointer(x), JuliaPointer(y)));
  end
);

InstallMethod( \-,
  ["IsAAGenericMatrixObj", "IsAAGenericMatrixObj"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\-(JuliaPointer(x), JuliaPointer(y)));
  end
);

InstallMethod( \*,
  ["IsAAGenericMatrixObj", "IsAAGenericMatrixObj"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\*(JuliaPointer(x), JuliaPointer(y)));
  end
);

InstallMethod( \^,
  ["IsAAGenericMatrixObj", "IsPosInt"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\^(JuliaPointer(x), y));
  end
);

InstallMethod( AdditiveInverse,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Base.\-(JuliaPointer(x)))
);

InstallMethod( Inverse,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.inv(JuliaPointer(x)))
);

InstallMethod( InverseMutable,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.inv(JuliaPointer(x)))
);

InstallMethod( IsZero,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Base.iszero(JuliaPointer(x))
);

InstallMethod( IsOne,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Base.isone(JuliaPointer(x))
);

InstallMethod( Zero,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.zero(Julia.Generic.parent(JuliaPointer(x))))
);

InstallMethod( One,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.one(Julia.Generic.parent(JuliaPointer(x))))
);

# Need ring elements in GAP for this...
#InstallMethod( TraceMat,
#  ["IsAAGenericMatrixObj"],
#  x -> Julia.Generic.tr(JuliaPointer(x))
#);
#
#InstallMethod( DeterminantMat,
#  ["IsAAGenericMatrixObj"],
#  x -> Julia.Generic.det(JuliaPointer(x))
#);

InstallMethod( TransposedMat,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.transpose(JuliaPointer(x)))
);

################################################################################
#
#  Matrix rows
#
################################################################################

# GAP likes it matrices to be lists of lists. We fake this up by returning views
# to the rows of the MatElems. But still, this is horrible...

InstallMethod( \[\],
  ["IsAAGenericMatrixObj", "IsPosInt"],
  function( x, y )
    return Julia.GAPGroups.oscar_matrix_row_to_gap(JuliaPointer(x), y);
  end
);

InstallMethod( \=,
  ["IsAAGenericMatrixRowObj", "IsAAGenericMatrixRowObj"],
  function( x, y )
    return Julia.Base.\=\=(JuliaPointer(x), JuliaPointer(y));
  end
);

InstallMethod( Length,
  ["IsAAGenericMatrixRowObj"],
  x -> Julia.Nemo.ncols(JuliaPointer(x))
);
