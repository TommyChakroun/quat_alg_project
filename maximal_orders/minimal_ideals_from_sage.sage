load("maximal_orders/maximal_orders_utilities.sage")
load("utilities/algebra_type.sage")
load("utilities/utilities.sage")


#------------------------------------------------------------------------------------------
#
#            MINIMAL IDEALS FROM INDEMPOTENTS
#         COPY OF THE CODE OF https://github.com/sagemath/sage/blob/develop/src/sage/categories/finite_dimensional_semisimple_algebras_with_basis.py#L66
#         with fix
#-----------------------------------------------------------------------------------------

from sage.categories.finite_dimensional_semisimple_algebras_with_basis import FiniteDimensionalSemisimpleAlgebrasWithBasis


def orthogonal_decomposition_perso(A, generators=None):
    if A.dimension() == 1:
        return A.basis().list()

    category = Algebras(A.base_ring()).Semisimple().WithBasis().FiniteDimensional().Commutative().Subobjects()

    if generators is None:
        generators = A.basis().list()

    # Searching for a good generator ...
    for gen in generators:
        # Computing the eigenspaces of the linear map x ↦ gen*x
        print("HERE !")
        phi = A.module_morphism(on_basis=lambda i: gen * A.term(i), codomain=A)
        eigenspaces = phi.matrix().eigenspaces_right(format='galois')
        print("EIgen sapce compute :")
        print(eigenspaces)

        if len(eigenspaces) >= 2:
            # Split the algebra according to the eigenspaces
            subalgebras = [
                A.submodule(map(A.from_vector, eigenspace.basis()), category=category)
                for eigenvalue, eigenspace in eigenspaces
            ]
            # Recursively decompose each eigenspace
            return tuple([
                idempotent.lift()
                for subalgebra in subalgebras
                for idempotent in orthogonal_decomposition_perso(subalgebra)
            ])

    raise Exception("Unable to fully decompose %s!" % A)


def central_orthogonal_idempotents_commutative_perso(A):
    return tuple([
        (e.leading_coefficient() / (e * e).leading_coefficient()) * e
        for e in orthogonal_decomposition_perso(A)
    ])


def center_perso(A):
    """
    Return the center of a finite dimensional semisimple algebra A
    as a submodule of A, avoiding hashability issues.
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


def central_orthogonal_idempotents_perso(A):
    Z,lift = center_perso(A)
    return tuple([lift(x) for x in central_orthogonal_idempotents_commutative_perso(Z)])


def minimal_ideals_sage(A):
    basis_A = list(A.basis())
    F = A.base_ring()
    dim_A = dimension(A)

    CentralIdempotents = central_orthogonal_idempotents_perso(A)
    
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
