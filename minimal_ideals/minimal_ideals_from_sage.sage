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
