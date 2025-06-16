load("utilities/utilities.sage")

#------------------------------------------------------------------------------------------
#
#            UTILITIES ON ALGEBRAS TYPE
#
#------------------------------------------------------------------------------------------


##---List of functions---## 

# finite_dimensional_algebra_format(A, basis)
# tensor(A, B, basis_A, basis_B)
# opposite(B, basis_B)
# quat_alg_mixed_with_table(a,b,F,transformation = False)




def structure_constants(A, basis_A,mat_inv = None):
    """
    INPUT: 
        -- A -- an algebra
        -- basis -- a list of element of A which is a basis of A
    OUTPUT:
        -- table -- the table associate to basis_A
    """
    F = A.base_ring()
    n = dimension(A)

    if mat_inv == None :
        mat_new_basis = change_matrix(A,basis_A)  # matrix of basis_A in the natural basis A.basis() of A
        mat_inv = mat_new_basis.inverse()

    table = []
    for e in basis_A:
        rows = []
        for f in basis_A:
            product = f*e
            rhs = vector(F, get_coefficients(product,A))
            coords = rhs*mat_inv
            rows.append(coords)
        row_matrix = Matrix(F, rows)
        table.append(row_matrix)

    return table


def finite_dimensional_algebra_format(A, basis_A,semisimple = False):
    """
    INPUT: 
        -- A -- an algebra
        -- basis -- a list of element of A which is a basis of A
    OUTPUT:
        -- B -- an instance of FiniteDimensionalAlgebra with the mutlipication table 
                the table of the multiplication of A in basis_A
                (so B isomorphic to A)
    """
    F = A.base_ring()
    cat = AlgebrasWithBasis(F).Associative().FiniteDimensional()
    if semisimple:
        cat = cat = AlgebrasWithBasis(F).Associative().FiniteDimensional().Semisimple()
    return  FiniteDimensionalAlgebra(F, structure_constants(A, basis_A),category =cat )

## Utilities to fix type issue on the original Sage code

def structure_constants_subspace(A, basis_B):
    """
    INPUT: 
        -- A -- an algebra
        -- basis_B -- a list of element of A which is a basis of a subspace of B which is a subalgebra
    OUTPUT:
        -- table -- the table associate to basis_F
    """
    F = A.base_ring()
    n = dimension(A)

    mat_new_basis = change_matrix(A,basis_B)  # matrix of basis_A in the natural basis A.basis() of A
    table = []
    for e in basis_B:
        rows = []
        for f in basis_B:
            product = f*e
            rhs = vector(F, get_coefficients(product,A))
            coords = mat_new_basis.solve_left(rhs)
            rows.append(coords)
        row_matrix = Matrix(F, rows)
        table.append(row_matrix)

    return table


def subalgebra_from_subspace(A,basis_B):
    """
    INPUT :
        -- A -- an algebra given by structure constants
        -- Basis_B -- a basis of a subspace of A which is assume to be also a subalgebra
    OUTPUT :
        -- B -- the algebra F given by structure constant in the basis basis_B
        -- lift -- a function from F to A to lift element
    """
    F = A.base_ring()
    cat = Algebras(F).Semisimple().WithBasis().FiniteDimensional().Commutative()
    table = structure_constants_subspace(A,basis_B)
    B = FiniteDimensionalAlgebra(F,table, category = cat)
    dim_B = dimension(B)
    lift = lambda b: sum(c * v for c, v in zip(get_coefficients(b, B), basis_B))
    return B,lift


def tensor(A, B):
    """
    INPUT :
        -- A,B -- two algebra of one of the type MatrixSpace, QuaternionAlgebra,FiniteDimensionalAlgebra
                n = dim A
                m = dim B
    OUTPUT :
        -- C -- the tensor product of A and B given of the type FiniteDimesnionalAlgebra in the basis c0,..,c_{nm-1} given by
                e0⊗f0,e0⊗f1,...,e0⊗f_{m-1},e1⊗f0,..,e1⊗f_{m-1},..,e_{n-1}⊗f0,...,e_{n-1}⊗f_{m-1}
        -- tensor_index_map -- a dictionary mapping (i,j) to the index k in the tensor basis corresponding to basis_A[i]⊗basis_B[j]
    """
    F = A.base_ring()
    n = dimension(A)
    m = dimension(B)
    basis_A = A.basis()
    basis_B = B.basis()

    table = []

    for i in range(n):
        e = basis_A[i]
        for j in range(m):
            f = basis_B[j]
            rows = []
            for i1 in range(n):
                e1 = basis_A[i1]
                for j1 in range(m):
                    f1 = basis_B[j1]
                    # (e1⊗f1)(e⊗f)=e1 e ⊗ f1 f = sum_s,r c_s es ⊗ dr fr  (right multiplication by e⊗f )
                    x = e1 * e
                    y = f1 * f  
                    c = coordinate(x, A, basis_A)  
                    d = coordinate(y, B, basis_B)
                    coords = []
                    for s in range(n):
                        for r in range(m):  
                            coords.append(c[s] * d[r])
                    rows.append(coords)
            M = Matrix(F, rows)
            table.append(M)  # we don't transpose M because Sage is made weird

    C = FiniteDimensionalAlgebra(F, table)
    return C


def opposite(B):
    """
    INPUT : 
        -- B -- an algebra given in whatevertype
    OUTPUT :
        -- Bop -- the opposite algebra of B given in the type FiniteDimensionalAlgebra in the basis B.basis()
    """
    F = B.base_ring()
    n = dimension(B)
    basis_B = B.basis()

    table = []
    for e in basis_B:
        rows = []
        for f in basis_B:
            x = e * f  # right multiplication by e in the opposite algebra (multiplication reversed)
            coords = coordinate(x, B, basis_B)
            rows.append(coords)
        M = Matrix(F, rows)
        table.append(M)
    
    Bop = FiniteDimensionalAlgebra(F, table)
    return Bop


def quat_alg_mixed_with_table(F,a,b,transformation = False):
    """
    Return the quaternion algebra (a,b | F) view with a mixed basis
    """
    A = QuaternionAlgebra(F,a,b)
    P = random_invertible_matrix(F,3)

    one,i,j,k = A.basis()
    e0 = one
    e1 = P[0][0]*i + P[0][1]*j + P[0][2]*k
    e2 = P[1][0]*i + P[1][1]*j + P[1][2]*k
    e3 = P[2][0]*i + P[2][1]*j + P[2][2]*k

    new_basis = [e0,e1,e2,e3]
    B = finite_dimensional_algebra_format(A, new_basis)
    
    if transformation:
        return B,P
    return B
