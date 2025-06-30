#------------------------------------------------------------------------------------------
#
#            Some computation of number of matrices in M_n(F_p)
#             satisfying some property
#
# Note: This code is intended for a SageMath environment.
#------------------------------------------------------------------------------------------


#==========================================================================================
# Part 1: Number of irreducible polynomials
#==========================================================================================

def count_irr_pol(d, p):
    """
    Brute-force counting of monic irreducible polynomials in F_p[X] of degree d.
    This method is very slow and should only be used for verification.
    """
    F = GF(p)
    R = PolynomialRing(F, 'x')
    x = R.gen()
    cpt = 0
    # Iterate through all possible coefficient vectors for a monic polynomial of degree d
    for coeffs in F**d:
        # Construct polynomial P = x^d + c_{d-1}x^{d-1} + ... + c_0
        P = x**d + R.sum(coeffs[i] * x**i for i in range(d))
        if P.is_irreducible():
            cpt += 1
    return cpt


def formula_irr_pol(d, p):
    """
    Computes the number of monic irreducible polynomials of degree d in F_p[X]
    using the formula: (1/d) * sum_{k|d} mu(d/k) * p^k.
    """
    # The sum is guaranteed to be divisible by d. Use integer division for correctness.
    return sum(moebius(d // k) * p**k for k in divisors(d)) // d


#==========================================================================================
# Part 2: Number of matrices with an irreducible minimal polynomial
#==========================================================================================

def count_irr_mat(n, p):
    """
    Brute-force counting of matrices in M_n(F_p) with an irreducible minimal polynomial.
    This is computationally expensive and only feasible for very small n and p.
    """
    cpt = 0
    # Iterate through all p^(n^2) matrices in the space
    for M in MatrixSpace(GF(p), n):
        # M.minpoly() is also an alias for M.minimal_polynomial()
        pi = M.minimal_polynomial()
        if pi.is_irreducible():
            cpt += 1
    return cpt


def formula_irr_mat(n, p):
    """
    Computes the number of matrices in M_n(F_p) with an irreducible
    minimal polynomial using the derived formula.
    """
    # The total number of irreducible matrices is the sum over divisors d of n.
    total_count = 0
    
    # Pre-calculate the order of GL(n, F_p), as it's a common factor.
    # GL(n, p).order() is the canonical way to do this in Sage.
    order_gln_p = GL(n, p).order()
    
    for d in divisors(n):
        # Term 1: Number of monic irreducible polynomials of degree d
        num_irr_pol_d = formula_irr_pol(d, p)
        
        # Term 2: Number of matrices for each such polynomial.
        # This is |GL_n(F_p)| / |C(D_pi)|, where C(D_pi) is isomorphic to GL_{n/d}(F_{p^d}).
        
        # The field for the centralizer group is F_{p^d}, which has size q = p^d
        q = p**d
        # The dimension for the centralizer group is k = n/d
        k = n // d
        
        # The order of GL(k, F_q)
        order_glk_qd = GL(k, q).order()
        
        # The number of matrices similar to the canonical form for a poly of degree d
        num_matrices_per_poly = order_gln_p // order_glk_qd
        
        # Add the contribution for this degree d to the total count
        total_count += num_irr_pol_d * num_matrices_per_poly
        
    return total_count


#==========================================================================================
# Part 3: Probability Functions
#==========================================================================================


def prob_reducible_mat(n, p):
    """
    Calculates the probability that a random matrix in M_n(F_p) has a
    reducible minimal polynomial.
    """
    # Total number of matrices in the space is p^(n^2)
    total_matrices = p**(n**2)
    # Number of matrices with an IRREDUCIBLE minimal polynomial
    num_irr_mat = formula_irr_mat(n, p)
    # The probability of being irreducible is the ratio
    prob_irr = num_irr_mat / total_matrices
    # The probability of being reducible is 1 minus that
    return 1 - prob_irr


#==========================================================================================
# Example Usage
#==========================================================================================

# Script to generate the data for your LaTeX tables

print("="*40)
print("Data for the Detailed First Table")
print("="*40)
# List of (n,p) pairs you want to compute
table1_pairs = [(2,2), (2,3), (3,2), (3,3)]

for (n, p) in table1_pairs:
    num_irr = formula_irr_mat(n, p)
    total_mat = p**(n**2)
    prob_irr = num_irr / total_mat
    prob_red = 1 - prob_irr
    print(f"(n,p)=({n},{p}) | I(n,p)={num_irr} | Total={total_mat} | P(irr)={prob_irr.n(digits=4)} | P(red)={prob_red.n(digits=4)}")

print("\n" + "="*40)
print("Data for the Probability Table (P(reducible))")
print("="*40)
n_values = [2, 3, 4]
p_values = [2, 3, 5]

# Print header
header = f"{'p \\ n':<8}" + "".join(f"{n:<10}" for n in n_values)
print(header)
print("-" * len(header))

for p in p_values:
    row_str = f"{p:<8}"
    for n in n_values:
        # Using the function to get the probability of a reducible matrix
        prob = prob_reducible_mat(n, p)
        row_str += f"{prob.n(digits=4):<10}"
    print(row_str)


print("--------------")

for n in range(2,10):
    for p in primes_first_n(100):
        r = prob_reducible_mat(n,p).n()
        if r<0.6:
            print(f"(n,p) = {(n,p)} : {r}")

          

