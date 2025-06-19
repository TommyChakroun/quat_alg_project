load("utilities/algebra_type.sage")
load("utilities/utilities.sage")
load("src/maximal_orders/maximal_orders_utilities.sage")


#------------------------------------------------------------------------------------------
#
#            IDEMPOTENTS TO IDEAL
#      
#------------------------------------------------------------------------------------------

def center_perso(A):
    """
    INPUT :
        -- A -- an associative algebra
    OUTPUT :
        -- Z -- an finite dimensional algebra given by structure constant representing the center of A
        -- lift -- the inclusion Z to A
    """
    basis_Z = list(A.center_basis())
    F = A.base_ring()
    basis_A = A.basis()
    # Convert basis to a list of A-elements (if needed)

    table = []
    for e in basis_Z:
        rows = []
        for f in basis_Z:
            x = f*e
            rows.append(coordinate(x,A,basis_Z))
        M = Matrix(F,rows)
        table.append(M)

    cat = Algebras(F).Semisimple().WithBasis().FiniteDimensional().Commutative()
    Z = FiniteDimensionalAlgebra(F,table,category =cat)

    def lift(x):
        """
        x in Z -> view x in A
        """
        coords = get_coefficients(x,Z)
        n = len(basis_Z)
        res = sum(coords[i]*basis_Z[i] for i in range(n))
        return res
    
    return Z,lift


def idempotents_to_ideals(A,CentralIdempotents):
    """
    INPUT :
        -- A -- an algebra given by a basis and structure constants
        -- CentralIdempotents -- a list [e1,..,er] of element of A which are in the center of A
                                    1 = e1+..+er
                                    ei ^2 =ei 
                                    ei ej =0 if i !=j
    OUTPUT :
        -- MinIdealsList -- a list [BI1,..,BIr] of basis of Is := A es
    """
    basis_A = list(A.basis())
    F = A.base_ring()
    dim_A = dimension(A)
    
    MinIdealsList =  []
    for e in CentralIdempotents :
        SetGeneratorsIdeal = [a*e for a in basis_A]

        ## Write vector of SetGeneratorsIdeal in each row of a matrix
        rows = []
        for x in SetGeneratorsIdeal:
            rows.append(get_coefficients(x,A))
        M = Matrix(F,rows)
        ## Extract free family of SetGeneratorsIdeal

        basis_Ai = []
        independent_vectors = M.row_space().basis_matrix()
        for v in independent_vectors:
            basis_Ai.append(sum(v[i] * basis_A[i] for i in range(dim_A)))
    
        MinIdealsList.append(basis_Ai)
    return MinIdealsList