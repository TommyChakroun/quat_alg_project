load("utilities/algebra_type.sage")
load("maximal_orders/maximal_orders_utilities.sage")
load("maximal_orders/minimal_ideals_from_sage.sage")
load("tests_max_order/test_diverse_max_orders.sage")

#------------------------------------------------------------------------------------------
# DEBUGGING SAGE'S CENTRAL ORTHOGONAL IDEMPOTENTS FUNCTION
# Issue: A.central_orthogonal_idempotents() fails due to center computation problems
#------------------------------------------------------------------------------------------

# Set the finite field - change this line to use different fields
F = GF(5)

# Define algebra categories over finite field F
finite_dim_category = Algebras(F).FiniteDimensional().WithBasis()
semisimple_category = Algebras(F).Semisimple().FiniteDimensional().WithBasis()

# Create 4x4 matrix algebra over finite field F
matrix_space = MatrixSpace(F, 4)
structure_table = structure_constants(matrix_space, matrix_space.basis())

# Test different algebra constructions
algebra_basic = FiniteDimensionalAlgebra(F, structure_table)
basis_basic = algebra_basic.basis()

algebra_finite_dim = FiniteDimensionalAlgebra(F, structure_table, category=finite_dim_category)
basis_finite_dim = algebra_finite_dim.basis()

algebra_semisimple = FiniteDimensionalAlgebra(F, structure_table, category=semisimple_category)
basis_semisimple = algebra_semisimple.basis()

# Note: center_basis() works but returns tuple (hence the comma in output)
# Note: center() fails when converting center_basis() result to submodule
# This causes central_orthogonal_idempotents() to fail

#------------------------------------------------------------------------------------------
# WORKING CASE: COMMUTATIVE SEMISIMPLE ALGEBRA
# central_orthogonal_idempotents() works for commutative algebras
#------------------------------------------------------------------------------------------

# Define commutative semisimple category
commutative_category = Algebras(F).Semisimple().FiniteDimensional().WithBasis().Commutative()

# Create simple commutative algebra (1x1 matrix)
simple_structure = [Matrix(F, [[1]])]
commutative_algebra = FiniteDimensionalAlgebra(F, simple_structure, category=commutative_category)


B,random_basis = mixed_matrix_space(QQ,4)

Zbasis_I = B.basis()

Zbasis_O = left_order(B,Zbasis_I)

A,pi = finite_algebra_from_order(B,Zbasis_O,5)

