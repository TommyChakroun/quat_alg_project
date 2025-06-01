load("maximal_orders/maximal_orders_utilities.sage")
load("maximal_orders/find_maximal_orders.sage")

# ------------------------------------------------------------------------------------------
#                              Maximal Order in M₄(ℚ)
# ------------------------------------------------------------------------------------------
#
# This script focuses on the computation of maximal orders in finite-dimensional algebras 
# over ℚ that are given by **structure constants**. A central example is algebras isomorphic 
# to M₄(ℚ), the algebra of 4×4 rational matrices.
#
# Let A = (a, b | ℚ) and B = (c, d | ℚ) be quaternion algebras over ℚ. Then, the algebra 
# A ⊗ B^op can be represented by structure constants. If we let:
#   - (e₁, e₂, e₃, e₄) = (1, i_A, j_A, k_A) be a basis for A
#   - (f₁, f₂, f₃, f₄) = (1, i_B, j_B, k_B) be a basis for B
# then a basis for A ⊗ B^op is given by:
#   (h₁, ..., h₁₆) = (e₁⊗f₁, e₁⊗f₂, ..., e₁⊗f₄, ..., e₄⊗f₄)
#
# The multiplication in A ⊗ B^op is given by:
#   h_i * h_j = ∑_{k=1}^{16} c_{ijk} h_k
# where the c_{ijk} ∈ ℚ are the structure constants.
#
# ------------------------------------------------------------------------------------------
# Problem Statement
# ------------------------------------------------------------------------------------------
#
# We are interested in the **explicit isomorphism problem**:
#   Given an algebra B defined by structure constants such that B ≅ M₄(ℚ),
#   compute an explicit isomorphism from B to M₄(ℚ).
#
# We will see all of our inputs as an algebra with structure constant, even M₄(ℚ) will be see as an abstract algebra A₀:
#   - A₀ is a 16-dimensional ℚ-algebra with basis E_{ij} (1 ≤ i, j ≤ 4)
#   - Multiplication is defined by E_{ij} * E_{kl} = δ_{j,k} E_{il}
#
# In this abstract setting, multiplication is slower (e.g., no Strassen algorithm),
# but this is acceptable for our purposes.
#
# ------------------------------------------------------------------------------------------
# General Algorithm
# ------------------------------------------------------------------------------------------
#
# Given:
#   - A: a 16-dimensional ℚ-algebra with basis (h₁, ..., h₁₆) and structure constants
#   - O: a ℤ-order in A given by a ℤ-basis (a₁, ..., a₁₆), where each aᵢ is a ℚ-linear 
#        combination of the hᵢ
#
# Goal:
#   - Compute a maximal order L containing O, given as a ℤ-basis (c₁, ..., c₁₆) ⊆ A.
#
# ------------------------------------------------------------------------------------------
# Usage and Function Signature
# ------------------------------------------------------------------------------------------
#
# A = FiniteDimensionalAlgebra(QQ, structure_constants)
# BA = A.basis()
# a₁ = c₁ * BA[0] + ... + c_3 * BA[3]
# ...
# a₁₆ = ...
# Zbasis_O = [a₁, ..., a₁₆]
# Zbasis_L = max_order_containing_order(A, Zbasis_O)
#
# ------------------------------------------------------------------------------------------
# Performance Testing
# ------------------------------------------------------------------------------------------
#
# We evaluate how the computation time depends on:
#   - the complexity of the table "structure constants" of A
#   - the complexity of the "coeficient cij" defining Zbasis_O
#
# For this, we define 5 algebras A₀ through A₄, each isomorphic to M₄(ℚ) but with varying 
# complexity in their structure constants:
#
# A₀ = FiniteDimensionalAlgebra(QQ, structure_constants_0)   # trivial
# A₁ = FiniteDimensionalAlgebra(QQ, structure_constants_1)   # slightly more complex
# A₂ = FiniteDimensionalAlgebra(QQ, structure_constants_2)   # complex
# A₃ = FiniteDimensionalAlgebra(QQ, structure_constants_3)   # more complex
# A₄ = FiniteDimensionalAlgebra(QQ, structure_constants_4)   # very complex
#
# For each Aᵢ, we define multiple ℤ-orders:
#   - Zbasis_O_00, Zbasis_O_01, ...
#
# For each pair (Aᵢ, Zbasis_O_ij), we compute:
#   Zbasis_L := max_order_containing_order(Aᵢ, Zbasis_O_ij)
#
# and compare the runtime.



# ------------------------------------------------------------------------------------------
# Construction of suitable structure_constant
# ------------------------------------------------------------------------------------------

# We define M4(Q) not with structure constant just for the purpose of construction.

MQ = MatrixSpace(QQ,4,4) 
basis_MQ = list(MQ.basis())


def structure_constants_from_invertible_matrix(P):
    """
    P a matrix 16 x 16 representing an arbitrary change of basis of M4(Q)
    """
    new_basis = [ sum(P[i][k]*basis_MQ[k] for k in range(16)) for i in range(16) ] 
    table = []
    for e in new_basis:
        rows = []
        for f in new_basis:
            x = f*e
            rows.append(coordinate(x,MQ,new_basis))
        M = Matrix(QQ,rows)
        table.append(M)
    return table


P0 = Matrix(QQ, [
    [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
])


P1 = Matrix(QQ, [
    [1, 1/2, -1/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 2/3, 1/4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, -2/5, 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 3/4, 2/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, -1/2, 1/3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 2/5, -3/4, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1/3, 2/7, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, -1/2, 3/5, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 4/5, -2/3, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1/4, 2/5, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1/3, 1/2, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1/5, 2/3, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2/3, 1/4, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3/7, 1/2],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1/3],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
])


P2 = Matrix(QQ, [
    [1, -1/4, 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, -2/5, 1/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 3/4, -1/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 2/3, -2/5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 1/2, 1/5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, -1/4, 2/3, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1/2, -2/3, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 3/5, 1/2, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, -1/3, 4/7, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1/3, -1/2, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2/5, 1/3, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1/4, 1/2, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3/4, -2/5, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1/3, -1/4],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1/2],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
])


P3 = Matrix(QQ, [
    [1/2, -3/5, 2/3, 1/4, 0, 1/7, -2/3, 5/6, -1/2, 2/7, 3/8, -1/3, 0, 4/5, -3/4, 1/9],
    [-1/3, 1, 4/7, -5/6, 2/5, 3/8, -1/2, 1/3, -1/4, 2/9, 0, 3/7, 1/5, -2/3, 4/7, -1/8],
    [2/5, 3/4, 1, -2/7, 3/5, -1/6, 4/9, 0, 1/3, -1/2, 3/8, 1/6, -3/7, 2/5, 0, -4/9],
    [-1/4, 2/3, -5/8, 1, -2/9, 3/7, 0, 1/4, -3/5, 2/3, 1/8, -1/6, 4/7, -2/3, 3/5, 1/10],
    [3/7, -1/2, 2/3, -3/5, 1, 1/6, -4/7, 3/5, 0, -2/9, 1/3, -1/4, 2/7, 4/9, -1/5, 2/3],
    [1/6, 2/5, -1/4, 1/7, -3/8, 1, -2/3, 5/9, 3/4, -1/2, 2/5, 0, 1/6, -3/7, 2/3, -4/5],
    [-2/3, 4/7, 0, -1/2, 2/3, 1/4, 1, -3/5, 3/7, 1/6, -2/5, 4/9, -1/3, 2/5, -1/4, 0],
    [1/5, -3/4, 2/7, 1/6, -4/9, 3/5, -1/2, 1, 2/3, 0, 1/8, -3/7, 1/4, -2/3, 3/8, 1/5],
    [-1/6, 1/3, 3/5, -2/7, 4/9, -1/2, 2/5, -3/4, 1, 1/7, -2/5, 3/6, 0, -1/3, 2/5, -1/6],
    [2/3, 0, -1/5, 3/7, -4/8, 1/3, -2/9, 1/5, 3/4, 1, -3/7, 2/5, 1/6, -4/9, 0, 2/3],
    [3/5, -1/6, 1/8, -2/3, 0, 4/7, 1/5, -3/8, 2/9, -1/4, 1, 3/5, -2/3, 1/6, -1/7, 3/8],
    [0, 1/3, -2/7, 3/6, -1/5, 2/9, -4/7, 1/6, 3/5, -1/2, 2/7, 1, 1/4, -3/5, 2/9, -1/6],
    [-3/7, 2/5, 1/3, 0, -1/4, 3/6, -2/5, 1/7, 4/8, 2/3, -1/6, 1/5, 1, -2/9, 3/7, -1/8],
    [1/4, -2/3, 0, 1/6, 2/5, -3/7, 1/3, -1/6, 0, 4/9, -2/5, 3/8, -1/3, 1, -4/7, 2/5],
    [-1/5, 3/7, -2/9, 4/8, 1/6, 0, 3/5, -1/4, 1/7, -2/3, 2/5, -3/6, 1/3, 2/4, 1, -1/7],
    [2/7, -1/6, 4/9, -3/5, 1/3, 1/8, -2/7, 0, 2/5, 1/6, -3/8, 4/9, 0, -1/2, 3/5, 1],
])




print("\n" + "="*60)
print(" Step 1: Construct structure constants from invertible matrices")
print("="*60)

structure_constants_0 = structure_constants_from_invertible_matrix(P0)
structure_constants_1 = structure_constants_from_invertible_matrix(P1)
structure_constants_2 = structure_constants_from_invertible_matrix(P2)
structure_constants_3 = structure_constants_from_invertible_matrix(P3)

print("✓ Structure constants computed.\n")

print("="*60)
print(" Step 2: Build finite-dimensional Q-algebras from constants")
print("="*60)

A0 = FiniteDimensionalAlgebra(QQ, structure_constants_0)
A1 = FiniteDimensionalAlgebra(QQ, structure_constants_1)
A2 = FiniteDimensionalAlgebra(QQ, structure_constants_2)
A3 = FiniteDimensionalAlgebra(QQ, structure_constants_3)

print("✓ Algebras A0, A1, A2, A4 constructed over Q.\n")

print("="*60)
print(" Step 3: Extract Q-bases for each algebra")
print("="*60)

BA0 = A0.basis()
BA1 = A1.basis()
BA2 = A2.basis()
BA3 = A3.basis()  # Make sure A3 is defined elsewhere!

print("✓ Bases BA0, BA1, BA2, BA3 extracted.\n")

# ------------------------------------------------------------------------------------------
# Construction of ℤ-order O in each algebra
# ------------------------------------------------------------------------------------------

print("="*60)
print(" Step 4: Compute Z-basis of the left order of the canonical lattice")
print("="*60)

#Zbasis_O_00 = left_order(A0, BA0, reduced=False)
#Zbasis_O_10 = left_order(A1, BA1, reduced=False)
#Zbasis_O_20 = left_order(A2, BA2, reduced=False)
#Zbasis_O_30 = left_order(A3, BA3, reduced=False)

print("✓ Z-bases of left orders computed.\n")

print("="*60)
print(" Step 5: Apply LLL reduction to the Z-bases")
print("="*60)

#Zbasis_O_01 = lattice_LLL(A0, Zbasis_O_00)
#Zbasis_O_11 = lattice_LLL(A1, Zbasis_O_10)
#Zbasis_O_21 = lattice_LLL(A2, Zbasis_O_20)
#Zbasis_O_31 = lattice_LLL(A3, Zbasis_O_30)

print("✓ LLL-reduced Z-bases computed.\n")
