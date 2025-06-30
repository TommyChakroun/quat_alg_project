# SageMath Script
#
# TITLE: Formal Isomorphism of Quaternion Algebras and Representation in a Matrix Ring

# =============================================================================
#  0. LOAD EXTERNAL UTILITIES
# =============================================================================
# The user must ensure the following utility files are present and contain
# the necessary function definitions for:
# - matrix_ring_iso_from_algebra_iso
# - is_algebra_homomorphism (optional, for verification)
# - is_linear_iso (optional, for verification)
load("utilities/algebra_type.sage")
load("utilities/utilities.sage")
load("src/isomorphism/explicit_iso_matrix_ring.sage")

# =============================================================================
#  1. SETUP FORMAL RING AND RELATIONS
# =============================================================================

R = PolynomialRing(QQ, 'a,b,c,d,c11,c12,c13,c21,c22,c23')
a,b,c,d,c11,c12,c13,c21,c22,c23 = R.gens()

# Define the relations that must hold for the homomorphism to be valid.
relation_1 = c*c11^2 + d*c12^2 - c*d*c13^2 - a
relation_2 = c*c21^2 + d*c22^2 - c*d*c23^2 - b
relation_3 = c*c11*c21 + d*c12*c22 - c*d*c13*c23

# Create a quotient ring where these relations are zero.
I = R.ideal([relation_1, relation_2, relation_3])
R_q = R.quotient(I,'a,b,c,d,c11,c12,c13,c21,c22,c23')

# The field of fractions of the quotient ring is our base field.
F = R_q.fraction_field()

# =============================================================================
#  2. DEFINE ALGEBRAS AND THE FORMAL HOMOMORPHISM
# =============================================================================
a,b,c,d,c11,c12,c13,c21,c22,c23 = F.gens()

# Define the quaternion algebras over the final field F.
A = QuaternionAlgebra(F, a, b)
B = QuaternionAlgebra(F, c, d)

oneA, iA, jA, kA = A.basis()
oneB, iB, jB, kB = B.basis()

# Define the images of the basis vectors of A in B.
f_iA = c11*iB + c12*jB + c13*kB
f_jA = c21*iB + c22*jB + c23*kB
f_kA = f_iA * f_jA # f(k) is determined by f(i)*f(j)

# The isomorphism dictionary maps basis elements of A to their images in B.
isom_dict = {
    0: oneB,
    1:   f_iA,
    2:   f_jA,
    3:   f_kA
}

# =============================================================================
#  3. COMPUTE AND DISPLAY THE MATRIX REPRESENTATION
# =============================================================================
# Optional: Verify that the map is a valid homomorphism.
try:
    print("Is the map a valid algebra homomorphism?", is_algebra_homomorphism(isom_dict, A, B))
    print("Is the map a linear isomorphism?", is_linear_iso(isom_dict, A, B))
except NameError:
    print("Verification functions (is_algebra_homomorphism, is_linear_iso) not found. Skipping checks.")

# Compute the matrix representation of phi: A tensor B^op -> M4(F).
phi_dict = matrix_ring_iso_from_algebra_iso(A, B, isom_dict)

# Print the final dictionary of matrices.
if phi_dict:
    print("\n--- Resulting Matrix Dictionary ---")
    for (basis_A, basis_B), matrix_rep in sorted(phi_dict.items(), key=lambda x: (str(x[0][0]), str(x[0][1]))):
        print(f"\nMatrix for phi({basis_A} tensor {basis_B}):")
        pretty_print(matrix_rep)
else:
    print("\nResulting dictionary is empty or could not be computed.")

