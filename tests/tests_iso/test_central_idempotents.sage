load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/minimal_ideals/minimal_ideals_manually.sage")
load("src/minimal_ideals/minimal_ideals_from_sage.sage")
load("src/minimal_ideals/minimal_ideals_from_magma.sage")

#------------------------------------------------------------------------------------------
# DEBUGGING SAGE'S CENTRAL ORTHOGONAL IDEMPOTENTS FUNCTION
# Issue: A.central_orthogonal_idempotents() fails due to center computation problems
#------------------------------------------------------------------------------------------

# Set the finite field - change this line to use different fields
#F = GF(5)

# Define algebra categories over finite field F
#finite_dim_category = Algebras(F).FiniteDimensional().WithBasis()
#semisimple_category = Algebras(F).Semisimple().FiniteDimensional().WithBasis()

# Create 4x4 matrix algebra over finite field F
#matrix_space = MatrixSpace(F, 4)
#structure_table = structure_constants(matrix_space, matrix_space.basis())

# Test different algebra constructions
#algebra_basic = FiniteDimensionalAlgebra(F, structure_table)
#basis_basic = algebra_basic.basis()

#algebra_finite_dim = FiniteDimensionalAlgebra(F, structure_table, category=finite_dim_category)
#basis_finite_dim = algebra_finite_dim.basis()

#algebra_semisimple = FiniteDimensionalAlgebra(F, structure_table, category=semisimple_category)
#basis_semisimple = algebra_semisimple.basis()

# Note: center_basis() works but returns tuple (hence the comma in output)
# Note: center() fails when converting center_basis() result to submodule
# This causes central_orthogonal_idempotents() to fail

#------------------------------------------------------------------------------------------
# WORKING CASE: COMMUTATIVE SEMISIMPLE ALGEBRA
# central_orthogonal_idempotents() works for commutative algebras
#------------------------------------------------------------------------------------------

# Define commutative semisimple category
#commutative_category = Algebras(F).Semisimple().FiniteDimensional().WithBasis().Commutative()

# Create simple commutative algebra (1x1 matrix)
#simple_structure = [Matrix(F, [[1]])]
#commutative_algebra = FiniteDimensionalAlgebra(F, simple_structure, category=commutative_category)


#B,random_basis = mixed_matrix_space(QQ,4)

#Zbasis_I = B.basis()

#Zbasis_O = left_order(B,Zbasis_I)

#A,pi = finite_algebra_from_order(B,Zbasis_O,5)



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




def test_16(F, n, split = False,nb_ite = 100):
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