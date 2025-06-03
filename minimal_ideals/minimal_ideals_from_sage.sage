load("maximal_orders/maximal_orders_utilities.sage")
load("utilities/algebra_type.sage")
load("utilities/utilities.sage")
load("minimal_ideals/idempotent_to_ideals.sage")


#------------------------------------------------------------------------------------------
#
#            MINIMAL IDEALS FROM INDEMPOTENTS
#         COPY OF THE CODE OF https://github.com/sagemath/sage/blob/develop/src/sage/categories/finite_dimensional_semisimple_algebras_with_basis.py#L66
#         with fix
#-----------------------------------------------------------------------------------------








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





## Correcttion of the copy of the original Sage Code

def orthogonal_decomposition_sage(A, generators=None):
    if A.dimension() == 1:
        return A.basis().list()

    if generators is None:
        generators = A.basis().list()

    BA = list(A.basis())
    dimA = dimension(A)
    # Searching for a good generator ...
    for gen in generators:
        # Computing the eigenspaces of the linear map x ↦ gen*x

        phi = A.module_morphism(on_basis=lambda i: gen * A.term(i), codomain=A)
        eigenspaces = phi.matrix().eigenspaces_right(format = 'galois')

        if len(eigenspaces) >= 2:
            # Split the algebra according to the eigenspaces
            subalgebras = [
                subalgebra_from_subspace(A,[sum(c[i]*BA[i] for i in range(dimA)) for c in eigenspace.basis()])
                for eigenvalue,eigenspace in eigenspaces
            ]
            # Recursively decompose each eigenspace
            return tuple([
                lift(idempotent)
                for subalgebra,lift in subalgebras
                for idempotent in orthogonal_decomposition_sage(subalgebra)
            ])

    raise Exception("Unable to fully decompose %s!" % A)


def idempotents_commutative_sage(A):
    return tuple([
        (e.leading_coefficient() / (e * e).leading_coefficient()) * e
        for e in orthogonal_decomposition_sage(A)
    ])


def central_idempotents_sage(A):
    Z,lift = center_perso(A)
    return tuple([lift(x) for x in idempotents_commutative_sage(Z)])





## Minimal ideals From Sage (don't work for all semisimple algebra)


def minimal_ideals_sage(A):
    return idempotents_to_ideals(A,central_idempotents_sage(A))
