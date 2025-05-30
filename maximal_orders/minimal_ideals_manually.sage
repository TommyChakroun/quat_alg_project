load("maximal_orders/maximal_orders_utilities.sage")
load("utilities/algebra_type.sage")


#------------------------------------------------------------------------------------------
#
#            MINIMAL IDEALS OF A SEMI SIMPLE ALGEBRA
#
#------------------------------------------------------------------------------------------


def center(A):
    """
    INPUT :
        -- A -- a finite dimensional algebre given by structur constant
    OUTPUT :
        -- basis_Z -- a basis of the center of A
        -- Z -- the finite dimensional algebra Z= Z(A) given by structure constant associated to basis_Z
    """
    F = A.base_ring()
    N = dimension(A)
    C = A.table()
    e = A.basis()

    rows = []
    for j in range(N):
        for k in range(N):
            rows.append([C[j][i][k]-C[i][j][k] for i in range(N)])


    
    M = Matrix(F,rows)
    basis_Ker_M = M.right_kernel().basis()


    basis_Z = []
    for X in basis_Ker_M:
        basis_Z.append(sum(X[i]*e[i] for i in range(N)))

    table = []
    for e in basis_Z:
        rows = []
        for f in basis_Z:
            x = f*e
            rows.append(coordinate(x,B,basis_Z))
        M = Matrix(F,rows)
        table.append(M)

    Z = FiniteDimensionalAlgebra(F,table)
    return basis_Z,Z


def minimal_polynomial_over_field(a,A,basis_field):
    """
        -- A -- a finite dimensional algebra over a field F
        -- a -- an element of A
        -- basis_field -- a list of element of 
    """

    # Get basis of the algebra A
    basis = A.basis()
    n = len(basis)

    # Multiplication by a as a matrix
    M = Matrix(F, n)
    for j, bj in enumerate(basis):
        prod = a * bj
        coeffs = get_coefficients(prod,A)
        M.set_column(j, vector(F, coeffs))

    # Minimal polynomial of the multiplication matrix
    return M.minimal_polynomial()


def split_commutative_algebra(A):
    """
    INPUT :
        -- A -- a finite dimensional commutative algebra over a finite field given by structure constant  NOT NECESSARY UNITARY
    OUTPUT :
        -- basis_I,basis_J -- basis of two strict two sided ideal of A such that I oplus J = A, None if A is simple
    """
    pass


def minimal_ideals_commutative(A):
    """
    INPUT :
        -- A -- a finite dimensional commutative algebra over a finite field given by structure constant
    OUTPUT :
        -- MinIdealsList -- a list of list of element of A which are the basis of the minimal two sided ideal of A
    """
    if split_commutative_algebra(A) == "IsField":
        return [A.basis()]
    
    basis_I,basis_J = split_commutative_algebra(A)

    ## We compute recursively the minimal ideal of I and J.
    
    I = finite_dimensional_algebra_format(A, basis_I)
    J = finite_dimensional_algebra_format(A, basis_J)

    dim_I = I.dimension()
    dim_J = J.dimension()

    list_basis_ideal_I = minimal_ideals_commutative(I)
    list_basis_ideal_J = minimal_ideals_commutative(J)

    list_basis_ideal_A = []

    for b in list_basis_ideal_I:
        lift_b = []
        for e in b:
            cooefs = get_coefficients(e,I)
            lift_b.append(sum(coefs[i]*basis_I[i] for i in range(dim_I)))

    for b in list_basis_ideal_J:
        lift_b = []
        for e in b:
            cooefs = get_coefficients(e,J)
            lift_b.append(sum(coefs[i]*basis_J[i] for i in range(dim_J)))

    return list_basis_ideal_A


def minimal_ideals_manually(A):
    """
    INPUT : 
        -- A -- a finite dimensional algebra over a finite field given by structure constant
    OUTPUT :
        -- MinIdealsList -- a list of list of element of A which are the basis of the minimal two sided ideal of A
    """
    F = A.base_ring()

    basis_Z,Z = center(A)
    basis_A = A.basis()

    dim_Z = dimension(Z)
    dim_A = dimension(A)
    

    ListBasisMinIdealCenter = minimal_ideals_commutative(Z)
    MinIdealsList = []

    for basis_ZAi_in_Z in ListBasisMinIdealCenter:   # visit all basis of Z(Ai) where Z = Z(A1) ⊕ ...  ⊕ Z(Ar)
        ## View the basis of Z(Ai) in A 
        basis_ZAi=[]
        for e in basis_ZAi_in_Z :
            coords = get_coefficient(e,Z)
            basis_ZAi.append(sum(coords[i]*basis_Z[i] for i in range(dim_Z)))

        ## Construct the set of generators of ideal Ai = Z(Ai)A :  {bj*ak for bj a basis of Z(Ai) and ak a basis of A }
        dim_ZAi = len(basis_ZAi)
        SetGeneratorsIdeal = []

        for j in range(dim_ZAi):
            for k in range(dim_A):
                SetGeneratorsIdeal.append(basis_ZAi[j]*basis_A[k])

        ## Write vector of SetGeneratorsIdeal in each row of a matrix
        rows = []
        for x in SetGeneratorsIdeal:
            rows.append(get_coefficient(x,A))
        M = Matrix(F,rows)
        ## Extract free family of SetGeneratorsIdeal

        basis_Ai = []
        independent_vectors = M.row_space().basis_matrix()
        for v in independent_vectors:
            basis_Ai.append(sum(v[i] * basis_A[i] for i in range(dim_A)))
        
        ## Add the basis of ideal to the list

        MinIdealsList.append(basis_Ai)

    return MinIdealsList

