load("utilities/utilities.sage")
load("utilities/algebra_type.sage")


#------------------------------------------------------------------------------------------
#
#        FIND RANK-ONE ELEMENT / IDEMPOTENT IN AN ALGEBRA ISOMORPHIC TO Mn(Fp)
#                             (using randomization)
#
#------------------------------------------------------------------------------------------


#==========================================================================================
#  THEORETICAL IDEA
#==========================================================================================
#
# The probability that a matrix M in Mn(Fp) has a reducible minimal polynomial
# can be computed by the formula:
#
#     1/p^(n^2) * ( p^(n^2) - sum_{k|n}{N_k(p) * |GL_n(Fp)| / |GL_{n/k}(Fp^k)|} )
#
# where N_k(p) is the number of irreducible polynomials of degree k in Fp[X].
#
# One can show that the probability that a random matrix M in Mn(Fp) has a
# reducible minimal polynomial is non-zero and reasonably high.
#
# The overall strategy is as follows:
#
#   1. By iterating random choices, we will eventually find a matrix `c` with a
#      reducible minimal polynomial, which allows us to construct a zero divisor `x`.
#
#   2. By repeating this process, we have a good chance of finding a zero divisor `x`
#      with matrix rank n-1.
#
#   3. For such an `x`, its minimal polynomial is of the form t*g(t). It follows that
#      y := g(x) has rank 1.
#
#   4. After finding a rank-one element `y`, we can solve the equation `aya = y` for `a`
#      in the algebra A.
#
#   5. The element e := ay will be a rank-one idempotent (e^2 = e).
#
#==========================================================================================


#==========================================================================================
#  IMPLEMENTATION
#==========================================================================================

#--------------------------------------------------------------------
#  Step 1: Find a Zero Divisor
#--------------------------------------------------------------------

def zero_divisor_MnFp(A):
    """
    Finds a zero divisor in a matrix algebra over a finite field.

    INPUT:
        - A: An algebra over a finite field Fp, isomorphic to Mn(Fp).

    OUTPUT:
        - An element x in A that is a zero divisor.

    ALGORITHM:
        - We find an element with a reducible minimal polynomial by random sampling.
          If a polynomial pi(t) is reducible, e.g., pi = f*g, then f(c) and g(c)
          are zero divisors.
        - For a deterministic method, see:
          Lajos Ronyai, "Simple Algebras Are Difficult" (1987), Section 5.
    """
    BA = list(A.basis())
    F = A.base_ring()
    ite_max = 100

    for _ in range(ite_max):
        # Pick a random element in the algebra
        c = sum(F.random_element() * b for b in BA)
        pi = minimal_polynomial(c, A)

        if not pi.is_irreducible():
            # If the minimal polynomial is reducible, we can construct a zero divisor.
            factors = pi.factor()
            g = factors[0][0]  # Take the first factor
            return g(c)

    raise ValueError("Zero divisor not found. The probability is low, try increasing ite_max.")


#--------------------------------------------------------------------
#  Step 2: Find a Rank-One Element
#--------------------------------------------------------------------

def rank_one_MnFp(A):
    """
    Finds a rank-one element in a matrix algebra over a finite field.

    INPUT:
        - A: An algebra over a finite field Fp, isomorphic to Mn(Fp).

    OUTPUT:
        - An element x in A with rank 1.

    ALGORITHM:
        - We repeatedly find a zero divisor `x` until it has rank n-1.
          (Note: rank is defined via dim(xA), where right_rank(x,A) = n*rank(x)).
        - For such an `x`, its minimal polynomial is divisible by t. We can write it
          as pi(t) = t * Q(t). The element y = Q(x) will have rank 1.
    """
    # Note: dimension() and right_rank() are assumed to be defined in utility files.
    dim_A = dimension(A)
    n = int(sqrt(dim_A))
    F = A.base_ring()
    ite_max = 100

    for _ in range(ite_max):
        x = zero_divisor_MnFp(A)
        x_rank_dim = right_rank(x, A)

        # Case 1: The zero divisor has rank n-1.
        if x_rank_dim == n * (n - 1):
            pi = minimal_polynomial(x, A)
            R = PolynomialRing(F, 'x')
            t = R.gen()
            # The minimal polynomial of a singular matrix has a factor of t.
            Q = pi // t
            return Q(x)

        # Case 2: The zero divisor already has rank 1.
        if x_rank_dim == n:
            return x

    raise ValueError("Rank-one element not found. The probability is low, try increasing ite_max.")


#--------------------------------------------------------------------
#  Step 3: Find a Rank-One Idempotent
#--------------------------------------------------------------------

def rank_one_idempotent_MnFp(A):
    """
    Finds a rank-one idempotent in a matrix algebra over a finite field.

    INPUT:
        - A: An algebra over a finite field Fp, isomorphic to Mn(Fp).

    OUTPUT:
        - An element e in A such that e is a rank-one idempotent (e^2 = e).
    """
    # Note: solve_xax_eq_x() is assumed to be defined in a utility file.
    x = rank_one_MnFp(A)
    # For a rank-one element x, there exists 'a' such that xax = x.
    # The element e = ax is then a rank-one idempotent.
    a = solve_xax_eq_x(x, A)
    return a * x
