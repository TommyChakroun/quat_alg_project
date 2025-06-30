load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/minimal_ideals/minimal_ideals_manually.sage")
load("src/minimal_ideals/minimal_ideals_from_sage.sage")
load("src/minimal_ideals/minimal_ideals_from_magma.sage")

#------------------------------------------------------------------------------------------
#  Test of my personal function for central orthogonal idempotent in semisimpel algebra
#------------------------------------------------------------------------------------------




def cartesian_product_algebra(F,n,basis):
    V = VectorSpace(F,n)
    mat_F_can = change_matrix(V,basis)
    mat_can_F = mat_F_can.inverse()
    table = []
    for e in basis:
        rows = []
        for f in basis:
            x = f.pairwise_product(e)
            coords = x*mat_can_F
            rows.append(coords)
        M = Matrix(F,rows)
        table.append(M)
    cat = Algebras(F).WithBasis().FiniteDimensional()
    return FiniteDimensionalAlgebra(F,table,category = cat)




def test_16(F, n, split = True,nb_ite = 100):
    """
    Test function to:
    - Generate a random invertible basis change matrix P for F^n.
    - Define the cartesian product algebra A over a new basis.
    - Compute idempotents in the algebra.
    - Express those idempotents back in the original vector space F^n.
    """
    def print_title(title):
        """Helper function to print formatted titles"""
        separator = "=" * len(title)
        print(separator)
        print(f" {title}")
        print(separator)
    
    V = VectorSpace(F, n)
    P = random_invertible_matrix(F, n)
    new_basis = [P * e for e in V.basis()]
    A = cartesian_product_algebra(F, n, new_basis)
    
    print_title("Test 16: Cartesian Product Algebra Analysis")
    
    print_title("New Basis")
    for i, basis_vector in enumerate(new_basis):
        print(f"b_{i+1} = {basis_vector}")
    
    print_title("Algebra A")
    print(A)
    
    if split:
        idempotents = idempotents_commutative_perso_split(A)
    else :
        idempotents = idempotents_commutative_perso(A,nb_ite)
    
    print_title("Idempotents in A")
    for i, e in enumerate(idempotents):
        print(f"e_{i+1} = {e}")
    
    print_title("Idempotents viewed in F^n")
    idempotents_in_Fn = [
        sum(c * b for c, b in zip(get_coefficients(e, A), new_basis))
        for e in idempotents
    ]
    for i, v in enumerate(idempotents_in_Fn):
        print(f"e_{i+1} in F^{n} = {v}")
    
    print("=" * 50)  # Final separator