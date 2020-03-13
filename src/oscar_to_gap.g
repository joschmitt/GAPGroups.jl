DeclareCategory("IsForeignObject", IsObject);
DeclareAttribute("ForeignPointer", IsForeignObject);

# AAGeneric as in AbstractAlgebra Generic
DeclareSynonym("IsAAGenericRing", IsForeignObject and IsRing);
DeclareSynonym("IsAAGenericRingElem", IsForeignObject and IsRingElement);
DeclareSynonym("IsAAGenericField", IsForeignObject and IsField);
DeclareSynonym("IsAAGenericMatrixObj", IsForeignObject and IsMatrixObj and IsMatrix and IsOrdinaryMatrix);
DeclareSynonym("IsAAGenericMatrixRowObj", IsForeignObject and IsVectorObj and IsList);

DeclareCategoryCollections("IsAAGenericRingElem");
DeclareCategoryCollections("IsAAGenericRingElemCollection");
DeclareCategoryCollections("IsAAGenericRingElemCollColl");
aageneric_family := NewFamily( "AAGeneric_Family", IsObject, IsAAGenericRingElem );
aageneric_cyclotomic_family := NewFamily( "AAGeneric_Cyclotomic_Family", IsObject, IsAAGenericRingElem and IsCyclotomic );
aageneric_algebraic_elements_family := NewFamily( "AAGeneric_Algebraic_Elements_Family", IsObject, IsAAGenericRingElem and IsAlgebraicElement );
aageneric_ffe_family := NewFamily( "AAGeneric_FFE_Family", IsObject, IsAAGenericRingElem and IsFFE );

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

InstallMethod( ViewObj,
  ["IsForeignObject"],
  function( x )
    Print( Julia.GAP.julia_to_gap(Julia.GAPGroups._print_to_string(ForeignPointer(x))) );
  end
);

################################################################################
#
#  Functions for RingElems
#
################################################################################

InstallMethod( \=,
  ["IsAAGenericRingElem", "IsAAGenericRingElem"],
  function( x, y )
    return Julia.Base.\=\=(ForeignPointer(x), ForeignPointer(y));
  end
);

InstallMethod( \=,
  ["IsAAGenericRingElem", "IsInt"],
  function( x, y )
    return Julia.Base.\=\=(ForeignPointer(x), Julia.GAP.gap_to_julia(y));
  end
);

InstallMethod( \=,
  ["IsInt", "IsAAGenericRingElem"],
  function( x, y )
    return Julia.Base.\=\=(Julia.GAP.gap_to_julia(x), ForeignPointer(y));
  end
);

InstallMethod( \+,
  ["IsAAGenericRingElem", "IsAAGenericRingElem"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\+(ForeignPointer(x), ForeignPointer(y)));
  end
);

InstallMethod( \-,
  ["IsAAGenericRingElem", "IsAAGenericRingElem"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\-(ForeignPointer(x), ForeignPointer(y)));
  end
);

InstallMethod( \*,
  ["IsAAGenericRingElem", "IsAAGenericRingElem"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\*(ForeignPointer(x), ForeignPointer(y)));
  end
);

InstallMethod( \*,
  ["IsAAGenericRingElem", "IsInt and IsZeroCyc"],
  SUM_FLAGS + 1, # Convince GAP that we really want this one...
  function( x, y )
    return Zero(x);
  end
);

InstallMethod( \*,
  ["IsInt and IsZeroCyc", "IsAAGenericRingElem"],
  SUM_FLAGS + 1, # Convince GAP that we really want this one...
  function( x, y )
    return Zero(x);
  end
);

InstallMethod( \^,
  ["IsAAGenericRingElem", "IsPosInt"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.\^(ForeignPointer(x), y));
  end
);

InstallMethod( AdditiveInverse,
  ["IsAAGenericRingElem"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Base.\-(ForeignPointer(x)))
);

InstallMethod( Inverse,
  ["IsAAGenericRingElem and IsMultiplicativeElementWithInverse"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.inv(ForeignPointer(x)))
);

InstallMethod( IsZero,
  ["IsAAGenericRingElem"],
  x -> Julia.Base.iszero(ForeignPointer(x))
);

InstallMethod( IsOne,
  ["IsAAGenericRingElem and IsMultiplicativeElementWithOne"],
  x -> Julia.Base.isone(ForeignPointer(x))
);

InstallMethod( Zero,
  ["IsAAGenericRingElem"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.zero(Julia.Generic.parent(ForeignPointer(x))))
);

InstallMethod( One,
  ["IsAAGenericRingElem and IsMultiplicativeElementWithOne"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Generic.one(Julia.Generic.parent(ForeignPointer(x))))
);

InstallOtherMethod( Conductor,
  ["IsForeignObject and IsCyclotomic and IsCyc"],
    x -> Julia.GAPGroups._conductor(ForeignPointer(x))
);

# This is just copied from fldabnum.gi.
# I had to replace the Conductor call by the Conductor method for lists in
# cyclotom.gi because GAP directly calls C for IsCyclotomicCollection.
InstallOtherMethod( FieldByGenerators,
  [ "IsAAGenericRingElemCollection and IsCyclotomicCollection" ],
  function( gens )
    local N, stab, entry;
    N := 1;
    for entry in gens do
      N := LcmInt( N, Conductor( entry ) );
    od;

    # Handle trivial cases.
    if N = 1 then
      return Rationals;
    elif N = 4 then
      return GaussianRationals;
    fi;

    # Compute the reduced stabilizer info.
    stab := Filtered( PrimeResidues( N ),
                   x -> ForAll( gens,
                                gen -> GaloisCyc( gen, x ) = gen ) );

    # Construct and return the field.
    return AbelianNumberFieldByReducedGaloisStabilizerInfo( Rationals, N,
             stab );
  end
);

InstallTrueMethod( IsIntegralCyclotomic, IsForeignObject and IsCyclotomic and IsInt);

InstallMethod( IsIntegralCyclotomic,
  ["IsForeignObject and IsCyclotomic"],
    x -> Julia.Hecke.isintegral(ForeignPointer(x))
);

InstallMethod( GaloisCyc,
  ["IsForeignObject and IsCyc", "IsInt"],
  function( x, k )
    return Julia.GAPGroups.oscar_to_gap(Julia.GAPGroups._galois_cyc(ForeignPointer(x), k));
  end
);

InstallMethod( Order,
  ["IsForeignObject and IsCyc and IsCyclotomic"],
  function( x )
    local n;
    n := Julia.GAPGroups._multiplicative_order(ForeignPointer(x));
    if IsZero(n) then
      return infinity;
    fi;
    return n;
  end
);

# Can't add methods to AbsInt because its a global function
#InstallGlobalFunction( AbsInt,
#  ["IsForeignObject and IsRat"],
#  x -> Julia.GAPGroups.oscar_to_gap(Julia.Hecke.abs(ForeignPointer(x)))
#);

InstallMethod( AbsoluteValue,
  ["IsForeignObject and IsRat"],
  x -> Julia.GAPGroups.oscar_to_gap(Julia.Hecke.abs(ForeignPointer(x)))
);

InstallMethod( \<,
  ["IsForeignObject and IsRat", "IsInt"],
  function( x, y )
    return Julia.Base.\<(ForeignPointer(x), Julia.GAP.gap_to_julia(y));
  end
);

InstallMethod( \<,
  ["IsInt", "IsForeignObject and IsRat"],
  function( x, y )
    return Julia.Base.\<(Julia.GAP.gap_to_julia(x), ForeignPointer(y));
  end
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
  x -> Julia.Hecke.nrows(ForeignPointer(x))
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

InstallMethod( \*,
  ["IsAAGenericMatrixObj and IsCyclotomicCollColl", "IsFFE"],
  function( M, z )
    local N;
    N := Julia.GAPGroups._cyclotomic_to_gapint(ForeignPointer(M));
    return N*z;
  end
);

InstallMethod( \*,
  ["IsFFE", "IsAAGenericMatrixObj and IsCyclotomicCollColl"],
  function( z, M )
    local N;
    N := Julia.GAPGroups._cyclotomic_to_gapint(ForeignPointer(M));
    return z*N;
  end
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

InstallOtherMethod( \[\],
  ["IsAAGenericMatrixRowObj", "IsPosInt"],
  function( x, y )
    return Julia.GAPGroups.oscar_to_gap(Julia.Base.getindex(ForeignPointer(x), 1, y));
  end
);

InstallMethod( \=,
  ["IsAAGenericMatrixRowObj", "IsAAGenericMatrixRowObj"],
  function( x, y )
    return Julia.Base.\=\=(ForeignPointer(x), ForeignPointer(y));
  end
);

InstallOtherMethod( Length,
  ["IsAAGenericMatrixRowObj"],
  x -> Julia.Hecke.ncols(ForeignPointer(x))
);

InstallMethod( NestingDepthA,
  ["IsAAGenericMatrixRowObj"],
  x -> 1
);

InstallMethod( NestingDepthM,
  ["IsAAGenericMatrixRowObj"],
  x -> 1
);

InstallMethod( \*,
  ["IsAAGenericMatrixRowObj", "IsAAGenericMatrixObj"],
  function( x, y )
    return Julia.GAPGroups.oscar_matrix_row_to_gap(Julia.Base.\*(ForeignPointer(x), ForeignPointer(y)), 1);
  end
);
