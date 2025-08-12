// This file contains two Magma functions:
// 1. `Equivalence`: Computes the equivalence of two quadratic forms using
//    a maximal isotropic subspace.
// 2. `IsoQuatAlg`: Computes an isomorphism between two quaternion algebras
//    based on the equivalence of quadratic forms.

// ============================================================================
// Function: Equivalence
//
// Computes a congruence matrix C such that S2 = C^T * S1 * C, using the
// computation of a maximal isotropic subspace.
//
// This function is based on the native Magma function `IsotropicSubspace(M)`
// for a symmetric matrix M. The method is inspired by the work of:
//
// - Denis Simon, "Solving Quadratic Forms in Dimension 4, 5 and more..."
// - John E. Watkins, "Some comments about indefinite LLL," Section 2.7.
//   (https://magma.maths.usyd.edu.au/~watkins/papers/illl.pdf)
//
// IMPORTANT ASSUMPTION:
// This implementation assumes that the function `IsotropicSubspace(M)`
// returns a **maximal** isotropic subspace. A future improvement would be
// to modify the source code of `IsotropicSubspace(M)` to allow for passing
// the factorization of the determinant, e.g., `IsotropicSubspace(M : FAC := ... )`.
// ============================================================================IsotropicSubspace(M) to allowed to give the factrozation of  
// det(M) like IsotropicSubspace(M : FAC := ... ).

function Equivalence(S1, S2 : FAC := [])
    if not (Nrows(S1) eq Nrows(S2)) then error "Matrices must have the same size."; end if;
    if not (IsSymmetric(S1) and IsSymmetric(S2)) then error "Matrices must be symmetric."; end if;
    F := BaseRing(S1);
    if not (F eq Rationals()) then error "Matrices must be over the Rational Field."; end if;

    n := Nrows(S1);

    // --- Paramters  FAC if there is an update for IsotropicSubspace(M) ---
    if FAC cmpeq [] then
        d1 := Determinant(S1);
        d2 := Determinant(S2);
        if d1 eq 0 or d2 eq 0 then
            error "One of the input matrices is singular.";
        end if;

        FAC:= Factorization(AbsoluteValue(d1))* Factorization(AbsoluteValue(d2));
    end if;


    M := DiagonalJoin(S1, -S2);
    W := IsotropicSubspace(M); // WILL BE BETTER to give FAC = FAC

    if Dimension(W) ne n then
        error "Failed to find a maximal isotropic subspace. Are the matrices congruent?";
    end if;
    L := BasisMatrix(W);

    A_int := Submatrix(L, 1, 1, n, n);
    B_int := Submatrix(L, 1, n + 1, n, n);

    A := ChangeRing(A_int, F);
    B := ChangeRing(B_int, F);

    if Determinant(A) eq 0 then
        error "The extracted block matrix A is singular. This can happen if the input forms are isotropic.";
    end if;
    
    C := A^-1 * B;

    return C;
end function;






// ============================================================================
// Function: IsoQuatAlg
//
// Computes an isomorphism between two quaternion algebras A and B.
// This function implements an improved method for finding an isomorphism
// when maximal orders and ramified primes are given.
//
// The method is based on the equivalence of quadratic forms, as described in
// the author's report:
// - Section 8: "Method 3: Isomorphism via Equivalence of Quadratic Forms"
// - Subsection 8.4: "Use of maximal orders."
// ============================================================================



function IsoQuatAlg(A,B : OA := false,OB := false,ram := false)
    a,b := -Norm(A.1),-Norm(A.2);
    c,d := -Norm(B.1),-Norm(B.2);

    DA := DiagonalMatrix(Rationals(),[-2*a,-2*b,2*a*b]);
    DB := DiagonalMatrix(Rationals(),[-2*c,-2*d,2*c*d]);


    if ram eq false then ram := RamifiedPrimes(A); end if;
    if OA eq false then OA := MaximalOrder(A); end if;
    if OB eq false then OB := MaximalOrder(B); end if;

    BTA := TraceZeroSubspace(OA); // actually a basis of trace zero subspace
    BTB := TraceZeroSubspace(OB); 

    rows := [Eltseq(e)[2..4] : e in BTA];
    NA := Matrix(rows);

    rows := [Eltseq(e)[2..4] : e in BTB];
    NB := Matrix(rows);

    SA := NA*DA*Transpose(NA);
    SB := NB*DB*Transpose(NB);


    fact := SeqFact([<r,4> : r in ram])*Factorization(4);


    print Determinant(SA)*Determinant(SB);
    print Facint(fact);


    M := Equivalence(SA, SB : FAC := fact);

    return Transpose(NA^-1 *M*NB);
end function;