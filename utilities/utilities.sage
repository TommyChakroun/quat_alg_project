#---------------------------------
#
#           UTILITIES
#
#---------------------------------

##---List of functions---##

# get_coefficients(x,B)
# is_in_field(x,B)
# dimension(A)
# coordinate(x,B,basis_subspace)
# evaluate(map_dict,A,B,x)
# is_linear_iso(map_dict,A,B)
# is_algebra_homomorphism(map_dict, A, B)
# random_invertible_matrix(F, n)
# matrix_in_basis(V, WW, bv, bw, g)
# read_line(string)

def get_coefficients(x, B):
    """
    Return the coefficients of x in the natural basis of B for SageMath
    """
    # Handle QuaternionAlgebra elements
    if hasattr(x, 'coefficient_tuple'):
        return list(x.coefficient_tuple())

    # Handle FiniteDimensionalAlgebra and other modules
    if hasattr(x, 'vector'):
        return list(x.vector())

    # Handle MatrixSpace elements
    if hasattr(x, 'list'):
        return x.list()

    raise TypeError(f"Don't know how to extract coefficients from object of type " + str(x)+str(B))


def is_in_field(x,B):
    """ 
    INPUT: 
        -- x -- an element of algebra B
        -- B -- an algebra over a field F, with first basis element e0 neutral.
    OUTPUT: 
        -- boolean : True if x is in F, False otherwise
    """
    n=B.dimension()
    for i in range(1,n):
        if x.vector()[i] !=0 :
            return False
    return True


def dimension(A):
    """
    Return the dimension of the algebra A. Work at least for A of the type :
    FiniteDimensionalAlgebra, MatrixSpace, QuaternionAlgebra
    INPUT :
        -- A -- an algebra
    OUTPUT :
        -- n -- the dimension of A.
    """
    try:
        return A.dimension()
    except:
        return 4


def change_matrix(B,basis_subspace):
    """
    Return the coordinate vector of x in the basis basis_subspace of a subspace of B, assume x lies in this subspace.
    INPUT:
        -- B -- an algebra with natural basis e1,..,en
        -- basis_subspace -- a list of lineary independant element of B, f1,..,fr
    OUTPUT :
        -- mat -- with r rows and n columne such that rows[i] is the coordinates of fi in e1,..,en
    """
    F = B.base_ring()
    mat = matrix(F, [vector(F, get_coefficients(b, B)) for b in basis_subspace])
    return mat


def coordinate(x, B, basis_subspace):
    """
    Return the coordinate vector of x in the basis basis_subspace of a subspace of B, assume x lies in this subspace.
    INPUT:
        -- x -- an element of B
        -- B -- an algebra
        -- basis_subspace -- a list of lineary independant element of B.
    OUTPUT :
        -- coords -- a list of element of F such that x = sum coords[i]*basis_subspace[i]
    """
    F = B.base_ring()
    mat = change_matrix(B,basis_subspace)
    rhs = vector(F, get_coefficients(x, B))
    return mat.solve_left(rhs)


def evaluate(map_dict,A,B,x):
    """
    Evaluate a linear map f : A -> B in x in A

    INPUT : 
        -- map_dict -- dictionary with key 0,..,n-1 where n=dim(A) and value element of B, map_dict[i] represent f(A.basis()[i]).
        -- x -- element in A
    OUTPUT
        -- f(x) -- in B
    """
    basis_A = A.basis()
    N = len(basis_A)
    coord = get_coefficients(x, A)
    result = B.zero()
    for i in range(N):
        result = result + coord[i] * map_dict[i]
    return result


def is_linear_iso(map_dict,A,B):
    """
    Check if the linear map f : A ->B is a linear isomorphism.
    INPUT : 
        -- map_dict -- dictionary with key 0,..,n-1 where n=dim(A) and value element of B, map_dict[i] represent f(A.basis()[i]).
        -- A -- algebra
        -- B -- algebra
    OUTPUT :
        -- boolean -- True if f is a linear isomorphism, False otherwise
    """
    F = A.base_ring()
    n = dimension(A)
 
    image = [map_dict[i] for i in range(n)]
    rows = []
    for b in image:
        rows.append(get_coefficients(b, B))
    M = Matrix(F, rows)
    return M.is_invertible()


def is_algebra_homomorphism(map_dict, A, B):
    """
    Check if the linear map f : A ->B is an algebra homomorphism.
    INPUT : 
        -- map_dict -- dictionary with key 0,..,n-1 where n=dim(A) and value element of B, map_dict[i] represent f(A.basis()[i]).
        -- A -- algebra
        -- B -- algebra
    OUTPUT :
        -- boolean -- True if f is an algebra homomrophism, False otherwise
    """

    basis_A = A.basis()
    N = dimension(A)

    f = lambda x : evaluate(map_dict, A, B, x)
    
    # Check if the map preserves the unit element
    if f(A.one()) != B.one():
        return False
    
    # Check if the map preserves multiplication
    for i in range(N):
        for j in range(N):
            ei = basis_A[i]
            ej = basis_A[j]
            if f(ei)*f(ej) != f(ei * ej):
                return False
    return True


def random_invertible_matrix(F, n):
    """
    Return a random invertible matrix of size n with coefficient in F.
    """
    i=0
    while True and i<1000:
        M = random_matrix(F, n)
        if M.is_invertible():
            return M
        i=i+1
    raise TypeError("invertible matrix not found")


def matrix_in_basis(V, W, bv, bw, g):
    """
    INPUT : 
        -- V -- vector space
        -- W -- vector space
        -- bv -- list 
        -- bw -- bases of E,F
        -- g -- real function define on bv taking value in W
    OUTPUT : 
        -- M -- the matrix of g in the bases bv,bw
    """
    F = V.base_ring()
    rows = []
    for e in bv:
        coords = coordinate(g(e), W, bw)
        rows.append(coords)
    M = Matrix(F, rows).transpose()
    return M


def minimal_polynomial(a,A):
    """
        -- A -- a finite dimensional algebra over a field F
        -- a -- an element of A
        -- basis_field -- a list of element of 
    """

    # Get basis of the algebra A
    basis = A.basis()
    n = len(basis)
    F = A.base_ring()

    # Multiplication by a as a matrix
    M = Matrix(F, n)
    for j, bj in enumerate(basis):
        prod = a * bj
        coeffs = get_coefficients(prod,A)
        M.set_column(j, vector(F, coeffs))

    # Minimal polynomial of the multiplication matrix
    return M.minimal_polynomial()



def extract_basis_from_generators(A, generators):
    """
    Given an algebra A over QQ and a list of elements (generators),
    return a list of linearly independent elements of A that span the same space.

    INPUT:
    - A: an algebra over QQ with finite dimension and structure constants.
    - generators: list of elements in A.

    OUTPUT:
    - List of elements in A forming a basis of the subspace generated by the input.
    """
    from sage.matrix.constructor import Matrix

    F = A.base_ring()
    basis_A = list(A.basis())
    dim_A = len(basis_A)

    # Convert each generator to its coordinate vector in the basis of A
    rows = [vector(F, get_coefficients(x,A)) for x in generators]

    # Create a matrix with each generator vector as a row
    M = Matrix(F, rows)

    # Extract a basis of the row space
    row_basis = M.row_space().basis_matrix()

    # Convert each basis vector back to an element of A
    result_basis = []
    for v in row_basis.rows():
        a = sum(v[i] * basis_A[i] for i in range(dim_A))
        result_basis.append(a)

    return result_basis
