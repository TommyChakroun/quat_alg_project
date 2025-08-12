load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/iso_splitting_algebra/maximal_orders/maximal_orders_utilities.sage")
load("src/iso_splitting_algebra/maximal_orders/find_maximal_orders.sage")
load("src/iso_splitting_algebra/rank_one/rank_one_MnFp.sage")


# =========================================================================================
#
#            FINDING A RANK-ONE ELEMENT IN AN ALGEBRA ISOMORPHIC TO Mn(Q)
#
# =========================================================================================


# =========================================================================================
#
# OVERVIEW:
# This script aims to find a rank-one element in a central simple algebra A over Q,
# where A is known to be isomorphic to a full matrix algebra Mn(Q). The process
# involves several steps:
#
# 1. Approximate Idempotent: First, find an approximate rank-one idempotent `e` by
#    working modulo a prime `p` (and later, modulo powers of `p` and multiple primes).
#    This element `e` behaves like a rank-one idempotent in the finite algebra A/pA.
#
# 2. Exact Zero Divisor: Use the approximate idempotent to construct an exact
#    zero divisor in a maximal order O of A (isomorphic to Mn(Z)). This provides a
#    non-invertible element within the integer-like structure of the algebra.
#
# 3. Rank-One Element: Transform the zero divisor into a true rank-one element. This
#    is a recursive process. Given a zero-divisor `x`, we can find an idempotent `e`
#    and consider the algebra `eAe`, which is isomorphic to `Mr(Q)` for `r < n`.
#    We then recursively find a rank-one element in this smaller algebra and lift
#    it back to the original algebra `A`.
#
# =========================================================================================


# ------------------------------------------------------------------------------------------
#
#   STEP 1: APPROXIMATE IDEMPOTENT IN Mn(Z)
#
# ------------------------------------------------------------------------------------------

def rank_one_idempotent_mod_p(A, p=5,Op = None):
    """
    Finds a rank-one idempotent in the algebra modulo a prime p.

    INPUT:
        -- A -- A Q-algebra given by structure constants and a basis (a1, ..., aN)
                such that O = Z*a1 + ... + Z*aN is a maximal order, and thus
                isomorphic to Mn(Z).
        -- p -- A prime number.
    OUTPUT:
        -- e -- An element `e` in O such that `e` modulo `p` is a rank-one idempotent
                when viewed in the finite algebra Op, which is isomorphic to Mn(Fp).
    """
    dim_A = dimension(A)
    BA = list(A.basis())

    if Op == None :
        Op, _ = finite_algebra_from_order(A, BA, p)

    e_mod_p = rank_one_idempotent_MnFp(Op)

    coords_bar = get_coefficients(e_mod_p, Op) # Get coefficients of e in Op's basis

    # Lift the element from Op back to the original algebra A
    e0 = A.sum(Integer(c) * b for c, b in zip(coords_bar, BA))
    return e0


def improve_idempotent_mod_p(A, p=5, n=4,Op = None):
    """
    Improves an approximate idempotent using Hensel's Lemma.

    INPUT:
        -- A -- An algebra over Q isomorphic to Mn(Q), assumed to be given by a natural
                basis of a maximal order, so A = Q*a1 + ... + Q*aN with
                O = Z*a1 + ... + Z*aN being a maximal order isomorphic to Mn(Z).
        -- p -- A prime number.
        -- n -- A natural number for the number of iterations.
    OUTPUT:
        -- en -- An element in A such that e_n mod p has rank one and e_n^2 - e_n = 0 mod (p^(2^n)).

    ALGORITHM:
        Starts from e0, a rank-one idempotent modulo p, and improves the solution
        using the iterative formula: e_{k+1} = 3*e_k^2 - 2*e_k^3.
    """
    e0 = rank_one_idempotent_mod_p(A, p,Op)
    ek = e0
    for _ in range(n):
        ek_squared = ek**2
        ek_cubed = ek_squared * ek
        ek = 3 * ek_squared - 2 * ek_cubed

    return ek


def idempotent_mod_several_primes(A, primes, n,Op):
    """
    Constructs an element that is an approximate idempotent modulo several primes.

    INPUT:
        -- A -- An algebra over Q isomorphic to Mn(Q), assumed to be given by a natural
                basis of a maximal order, so A = Q*a1 + ... + Q*aN with
                O = Z*a1 + ... + Z*aN being a maximal order isomorphic to Mn(Z).
        -- primes -- a list of primes
        -- n -- A natural number (the number of Hensel iterations).
    OUTPUT:
        -- e -- An element in A such that for each of the first `r` primes `pi`:
                - `e` mod `pi` is a rank-one idempotent in Opi (isomorphic to Mn(Fpi)).
                - e^2 - e = 0 mod (pi^(2^n)).

    ALGORITHM:
        Uses the Chinese Remainder Theorem to combine solutions.
    """
    dim_A = dimension(A)
    BA = list(A.basis())

    moduli = [p**(2**n) for p in primes]
    list_coords = []

    for p in primes:
        e = improve_idempotent_mod_p(A, p, n,Op)
        list_coords.append(get_coefficients(e, A))

    crt_coords = []
    for k in range(dim_A):
        remainders = [c[k] for c in list_coords]
        crt_coords.append(crt(remainders, moduli))

    x = A.sum(crt_coords[k] * BA[k] for k in range(dim_A))
    return x


# ------------------------------------------------------------------------------------------
#
#   STEP 2: EXACT ZERO DIVISOR IN Mn(Z)
#
# ------------------------------------------------------------------------------------------

def heuristic_zero_divisor_MnZ(A,primes,n_lifting_steps = 4,Op=None):
    """
    Constructs an exact zero divisor in an order isomorphic to Mn(Z).
    This function implements the Hensel lifting and LLL method.
    """
    BA = A.basis()
    N = dimension(A)

    # --- Step 1: Find an approximate idempotent ---

    # The result 'e' is a vector of N integer coefficients (x'_1, ..., x'_N).
    if len(primes)>1:
        e = idempotent_mod_several_primes(A, primes = primes, n=n_lifting_steps,Op=Op)
    else:
        e = improve_idempotent_mod_p(A,primes[0],n=n_lifting_steps,Op=Op)
    
    # This is the modulus K = p^(2^n) from our discussion.
    modulus_K = prod(p**(2**n_lifting_steps) for p in primes)

    # --- Step 2: Construct the basis for the lattice L in ZZ^(N+1) ---
    # The lattice is designed to find small integers c_1,...,c_N, d
    # satisfying c_k = d * x'_k (mod K).
    rows = []

    # Add the "approximation" vector: (x'_1, ..., x'_N, 1)
    row_x = list(e) + [1]
    rows.append(row_x)

    # Add the N "congruence" vectors: (K, 0,...), (0, K, 0,...), etc.
    for i in range(N):
        row_k = [0] * (N + 1)
        row_k[i] = modulus_K
        rows.append(row_k)
        
    # --- Step 3: Run LLL to find a short vector ---
    # This vector corresponds to the small integer solution (c_1, ..., c_N, d)
    L_basis = Matrix(ZZ, rows)
    L_reduced_basis = L_basis.LLL()


    for vect in L_reduced_basis:
        denominator = vect[N]
        if gcd(denominator,modulus_K) == 1:
            coeffs = [vect[i] / denominator for i in range(N)]
            e_exact = A.sum(c * b for c, b in zip(coeffs, BA))
            return e_exact

    raise ValueError("LLL returned a zero denominator. Try increasing precision (n_lifting_steps).")


def zero_divisor_MnZ(A,primes = [5],n=4):
    dim_A = dimension(A)
    BA = A.basis()

    p = primes[0]
    Op,_ = finite_algebra_from_order(A, BA, p)

    for ite in range(100):
        e = heuristic_zero_divisor_MnZ(A,primes=primes,n_lifting_steps=n,Op=Op)
        if right_rank(e,A)<dim_A:
            return e
    raise ValueError("Not found")



# ------------------------------------------------------------------------------------------
#
#   STEP 3: ZERO DIVISOR IN Mn(Q) VIA MAXIMAL ORDER
#
# ------------------------------------------------------------------------------------------

def zero_divisor_MnQ(A, BO=None):
    """
    Finds a zero divisor in an algebra isomorphic to Mn(Q).

    INPUT:
        -- A -- An algebra isomorphic to Mn(Q).
        -- BO -- A basis of a maximal order in A. If None, it will be computed.
    OUTPUT:
        -- x -- A zero divisor in A.
    """
    dim_A = dimension(A)

    if BO is None:
        BO = max_order(A)

    # Create a new algebra A_tilde with respect to the maximal order basis BO.
    A_tilde = FiniteDimensionalAlgebra(QQ, structure_constants(A, BO))

    x_in_A_tilde = zero_divisor_MnZ(A_tilde)

    coords_x_in_A_tilde = get_coefficients(x_in_A_tilde, A_tilde)

    # Convert the element back to the original algebra A
    return A.sum(coords_x_in_A_tilde[i] * BO[i] for i in range(dim_A))


# ------------------------------------------------------------------------------------------
#
#   STEP 4: RANK-ONE ELEMENT IN Mn(Q) FROM A ZERO DIVISOR
#
# ------------------------------------------------------------------------------------------

def rank_one_MnQ(A,zero_divisor = None):
    """
    Computes a rank-one element in an algebra isomorphic to Mn(Q).

    INPUT:
        -- A -- An algebra over Q isomorphic to Mn(Q).
    OUTPUT:
        -- x -- An element in A of "rank one", meaning dim(xA) = n.

    ALGORITHM:
        This is a recursive algorithm. Choose a zero divisor `x` in A, compute `y`
        such that `xy=0`. Replace `x` by `y` if `rank(y) < rank(x)`. Solve `xax=x`
        and let `e := ax` which has the same rank as `x`. Compute the algebra
        B := eAe. By theory, we know B is isomorphic to Mr(Q). By recursion,
        find a rank-one element `z` in B. Lift `z` from B to A.
    """
    dim_A = dimension(A)
    n = int(sqrt(dim_A))
    BA = list(A.basis())
    F = A.base_ring()

    if zero_divisor == None:
        x = zero_divisor_MnQ(A)
    else:
        x = zero_divisor

    # Let pi(t) be the minimal polynomial of (right multiplication by) x.
    # Since x is a zero divisor, pi(t) = t * Q(t). Let y = Q(x). Then xy=0.
    # If rank(x) > n^2/2, then rank(y) < n^2/2. We can replace x with y.
    pi = minimal_polynomial(x, A)
    R = PolynomialRing(F, 'x')
    t = R.gen()
    Q = pi // t
    if right_rank(x, A) > n * n / 2: # right_rank(A,x) is dim(xA)
        x = Q(x)

    if right_rank(x, A) == n:
        return x

    # Find an idempotent e=ax such that rank(e)=rank(x)
    a = solve_xax_eq_x(x, A)
    e = a * x

    # Construct the Peirce subalgebra B = eAe
    basis_eAe = extract_basis_from_generators(A, [e * b * e for b in BA])
    structure_constants_eAe = structure_constants_subspace(A, basis_eAe)
    B = FiniteDimensionalAlgebra(F, structure_constants_eAe)

    # Recursively find a rank-one element z in B
    z = rank_one_MnQ(B)

    # Lift z from B back into A.
    # CORRECTED: Changed get_coefficients(z, A) to get_coefficients(z, B)
    # as z is an element of B.
    z_coeffs_in_B = get_coefficients(z, B)
    z_lift_in_A = sum(c * b for c, b in zip(z_coeffs_in_B, basis_eAe))

    return z_lift_in_A
