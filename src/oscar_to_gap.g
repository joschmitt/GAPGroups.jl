DeclareCategory("IsForeignObject", IsObject);
DeclareAttribute("ForeignPointer", IsForeignObject);

# AAGeneric as in AbstractAlgebra Generic
DeclareSynonym("IsAAGenericRing", IsForeignObject and IsRing);
DeclareSynonym("IsAAGenericRingElem", IsForeignObject and IsRingElement);
DeclareSynonym("IsAAGenericField", IsForeignObject and IsField);
DeclareSynonym("IsAAGenericMatrixObj", IsForeignObject and IsMatrixObj and IsMatrix and IsOrdinaryMatrix);
DeclareSynonym("IsAAGenericMatrixRowObj", IsForeignObject and IsVectorObj);

DeclareCategoryCollections("IsAAGenericRingElem");
DeclareCategoryCollections("IsAAGenericRingElemCollection");
DeclareCategoryCollections("IsAAGenericRingElemCollColl");
aageneric_family := NewFamily( "AAGeneric_Family", IsObject, IsAAGenericRingElem );

# The general comparing function using the global GAP_comparison_dictionary
InstallMethod( \<,
  ["IsForeignObject", "IsForeignObject"],
  function( x, y )
    return Julia.Base.\<(Julia.GAPGroups._get_comparison_number(ForeignPointer(x)), Julia.GAPGroups._get_comparison_number(ForeignPointer(y)));
  end
);

################################################################################
#
#  Basic functions
#
################################################################################

InstallMethod( String,
  ["IsForeignObject"],
  x -> String( ForeignPointer( x ) )
);

InstallMethod( PrintObj,
  ["IsForeignObject"],
  x -> Print( PrintString( x ) )
);

InstallMethod( ViewString,
  ["IsForeignObject"],
  x -> ViewString( ForeignPointer(x) )
);

################################################################################
#
#  Basic functions for matrices
#
################################################################################

InstallMethod( NumberRows,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Generic.nrows(ForeignPointer(x))
);

InstallOtherMethod( Length,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Nemo.nrows(ForeignPointer(x))
);

InstallMethod( NumberColumns,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Generic.ncols(ForeignPointer(x))
);

InstallMethod( \[\],
  ["IsAAGenericMatrixObj", "IsPosInt", "IsPosInt"],
  function( x, y, z )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.getindex(ForeignPointer(x), y, z));
  end
);

InstallOtherMethod( RankMat,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Generic.rank(ForeignPointer(x))
);

InstallMethod( \=,
  ["IsAAGenericMatrixObj", "IsAAGenericMatrixObj"],
  function( x, y )
    return Julia.Base.\=\=(ForeignPointer(x), ForeignPointer(y));
  end
);

InstallMethod( \+,
  ["IsAAGenericMatrixObj", "IsAAGenericMatrixObj"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\+(ForeignPointer(x), ForeignPointer(y)));
  end
);

InstallMethod( \-,
  ["IsAAGenericMatrixObj", "IsAAGenericMatrixObj"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\-(ForeignPointer(x), ForeignPointer(y)));
  end
);

InstallMethod( \*,
  ["IsAAGenericMatrixObj", "IsAAGenericMatrixObj"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\*(ForeignPointer(x), ForeignPointer(y)));
  end
);

InstallMethod( \^,
  ["IsAAGenericMatrixObj", "IsPosInt"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\^(ForeignPointer(x), y));
  end
);

InstallMethod( AdditiveInverse,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Base.\-(ForeignPointer(x)))
);

InstallMethod( Inverse,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.inv(ForeignPointer(x)))
);

InstallMethod( InverseMutable,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.inv(ForeignPointer(x)))
);

InstallMethod( IsZero,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Base.iszero(ForeignPointer(x))
);

InstallMethod( IsOne,
  ["IsAAGenericMatrixObj"],
  x -> Julia.Base.isone(ForeignPointer(x))
);

InstallMethod( Zero,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.zero(Julia.Generic.parent(ForeignPointer(x))))
);

InstallMethod( One,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.one(Julia.Generic.parent(ForeignPointer(x))))
);

InstallOtherMethod( TraceMat,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.tr(ForeignPointer(x)))
);

InstallMethod( DeterminantMat,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.det(ForeignPointer(x)))
);

InstallOtherMethod( TransposedMat,
  ["IsAAGenericMatrixObj"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.transpose(ForeignPointer(x)))
);

InstallMethod( NestingDepthA,
  ["IsAAGenericMatrixObj"],
  x -> 2
);

InstallMethod( NestingDepthM,
  ["IsAAGenericMatrixObj"],
  x -> 2
);

################################################################################
#
#  Matrix rows
#
################################################################################

# GAP likes it matrices to be lists of lists. We fake this up by returning views
# to the rows of the MatElems. But still, this is horrible...

InstallOtherMethod( \[\],
  ["IsAAGenericMatrixObj", "IsPosInt"],
  function( x, y )
    return Julia.GAPGroups.oscar_matrix_row_to_gap(ForeignPointer(x), y);
  end
);

InstallMethod( \=,
  ["IsAAGenericMatrixRowObj", "IsAAGenericMatrixRowObj"],
  function( x, y )
    return Julia.Base.\=\=(ForeignPointer(x), ForeignPointer(y));
  end
);

InstallMethod( Length,
  ["IsAAGenericMatrixRowObj"],
  x -> Julia.Nemo.ncols(ForeignPointer(x))
);

InstallMethod( NestingDepthA,
  ["IsAAGenericMatrixRowObj"],
  x -> 1
);

InstallMethod( NestingDepthM,
  ["IsAAGenericMatrixRowObj"],
  x -> 1
);
