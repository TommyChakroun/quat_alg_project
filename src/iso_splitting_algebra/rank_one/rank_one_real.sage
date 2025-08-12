#
# This SageMath script implements algorithms to find an explicit isomorphism
# from a given central simple algebra A over Q to a matrix ring Mn(Q).
# The core strategy, based on the work of Ivanyos and Rónyai, involves:
#   1. Extending the algebra to a real number field E to find a rank-one element.
#   2. Using this rank-one element to define an embedding of A into Mn(E).
#   3. Approximating this embedding and using lattice reduction (LLL) to find a
#      "reduced" basis for a maximal order in A.
#   4. Searching for a rank-one element within this reduced basis, which corresponds
#      to a zero-divisor in the original algebra A.
#
# The script assumes the existence of helper functions from loaded files:
# - "utilities/utilities.sage": General helper functions.
# - "utilities/algebra_type.sage": Functions for algebra properties.
# - "src/isomorphism/explicit_iso_matrix_ring.sage": Functions for constructing isomorphisms.
#

load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/iso_splitting_algebra/explicit_iso_matrix_ring.sage")


# =========================================================================================
#
#   SECTION 1: FINDING A RANK-ONE ELEMENT IN AN EXTENDED ALGEBRA A ⊗ E
#
# =========================================================================================

def rank_one_MnR(A):
    """
    Finds a rank-one element in an algebra extended by a real number field.
    This is based on Lemma 2.4 from:
    https://www.ams.org/journals/mcom/1990-55-192/S0025-5718-1990-1035925-1/S0025-5718-1990-1035925-1.pdf

    INPUT:
        - A: A finite-dimensional algebra over Q, assumed to be isomorphic to Mn(Q).

    OUTPUT:
        - x: An element in the extended algebra B of rank one when viewed as a matrix in Mn(E).
        - B: The algebra A ⊗_Q E, where E = Q(gamma) is a real number field.
        Returns "Not found" if no suitable element is found within the iteration limit.

    ALGORITHM:
        Let the basis of A be {a_1, ..., a_N}.
        1. By the Schwartz-Zippel lemma, a random linear combination
           a = h_1*a_1 + ... + h_N*a_N (with small integer coefficients h_i)
           will, with high probability, have a minimal polynomial π_a that is
           separable, of degree n, and has at least one real root.
        2. We iterate to find such an element 'a'.
        3. We take an irreducible factor 'g' of its minimal polynomial π_a that has a real root.
        4. Let gamma be a real root of g, and define the number field E := Q(gamma).
        5. We extend the algebra to B := A ⊗_Q E. In the polynomial ring E[t], the
           minimal polynomial π_a now has (t-gamma) as a factor. We can write
           π_a(t) = (t-gamma) * Q(t).
        6. The element x := Q(a) in B is then a rank-one element.
    """
    N = A.dimension()
    BA = list(A.basis())
    n = int(sqrt(N))
    H = 2 * n * (n - 1)

    nb_ite = 1000
    for ite in range(nb_ite):
        h = [randint(1, H) for l in range(N)]
        a = A.sum(h[i] * BA[i] for i in range(N))

        pi = minimal_polynomial(a, A)

        if pi.degree() == n and gcd(pi, pi.derivative()) == 1:
            fact = [g for g, _ in pi.factor()]
            for g in fact:
                if len(g.real_roots()) > 0:
                    E = NumberField(g, 'gamma')
                    gamma = E.gen()
                    R = PolynomialRing(E, 't')
                    t = R.gen()
                    B = FiniteDimensionalAlgebra(E, A.table())
                    pi_E = R(pi)
                    Q = pi_E // (t - gamma)
                    a_in_B = B.sum(c * b for c, b in zip(get_coefficients(a, A), B.basis()))
                    x = Q(a_in_B)
                    return x, B
    return "Not found"


# =========================================================================================
#
#   SECTION 2: CONSTRUCTING AN EMBEDDING A -> Mn(E)
#
# =========================================================================================

def embed_R_N(A):
    """
    Constructs an embedding A -> Mn(E) from a rank-one element.

    INPUT:
        - A: A Q-algebra isomorphic to Mn(Q).

    OUTPUT:
        - BL: An N x N matrix (where N=n^2) with coefficients in a number field E.
              The i-th row of BL represents the image of the i-th basis element of A
              under an embedding φ: A -> Mn(E), flattened into a vector of length N.

    ALGORITHM:
        1. Call rank_one_MnR(A) to get a rank-one element 'x' in an extended
           algebra B = A ⊗ E.
        2. Use the helper function `matrix_ring_iso_from_rank_one` (assumed to be
           defined elsewhere) with B and x to get an explicit isomorphism φ from B
           to Mn(E).
        3. Apply this isomorphism to each basis element of A to compute its matrix
           image in Mn(E).
        4. Flatten each resulting n x n matrix into a row vector of length N=n^2.
        5. Collect these vectors into an N x N matrix BL.
    """
    dim_A = A.dimension()
    x, B = rank_one_MnR(A)
    E = B.base_ring()
    dico = matrix_ring_iso_from_rank_one(B, x)
    BL = Matrix(E, [vector(dico[i].list()) for i in range(dim_A)])
    return BL


# =========================================================================================
#
#   SECTION 3: RATIONAL APPROXIMATION AND LATTICE REDUCTION
#
# =========================================================================================

def approx_number_field_element(x, precision=1e-10):
    """
    Approximates an element of a real number field with a rational number.

    INPUT:
        - x: An element in a number field E = Q(gamma).
        - precision: The desired precision for the approximation.

    OUTPUT:
        - r: A rational number (element of QQ) that approximates the real
             value of x under a fixed real embedding of E.

    ALGORITHM:
        1. Get the defining polynomial 'g' of the number field E.
        2. Find the real roots of 'g'. An error is raised if none exist.
        3. Choose the first real root as the specific real value for the generator gamma.
           This fixes a particular real embedding of E into R.
        4. Evaluate the polynomial representation of 'x' at this real value of gamma.
        5. Convert the resulting floating-point number into a Sage rational number (QQ)
           with a number of digits determined by the desired precision.
    """
    E = x.parent()
    g = E.polynomial()
    real_roots = g.real_roots()
    if not real_roots:
        raise ValueError("The defining polynomial of E has no real roots.")

    gamma_real = real_roots[0]

    poly_rep = x.polynomial()
    real_val = poly_rep(gamma_real)

    num_digits = int(-log(precision, 10))
    r = QQ(real_val.n(digits=num_digits))
    return r


def frob_norm(a, A, BL):
    """
    Approximates the Frobenius (Euclidean) norm of the image of an algebra element.

    INPUT:
        - a: An element of the algebra A.
        - A: A Q-algebra isomorphic to Mn(Q).
        - BL: An N x N matrix representing the embedding A -> Mn(E), where the
              i-th row is the image of the i-th basis element.

    OUTPUT:
        - norm: A numerical approximation of the Frobenius norm of the matrix
                image of 'a' in Mn(E).

    ALGORITHM:
        1. Get the coordinates of 'a' with respect to the basis of A.
        2. Compute the image of 'a' under the embedding by taking a linear
           combination of the rows of BL, weighted by the coordinates of 'a'.
        3. Approximate each component of the resulting vector in E with a rational number.
        4. Calculate the square root of the sum of the squares of these components.
    """
    dim_A = A.dimension()
    coords = get_coefficients(a, A)
    a_embed = sum(coords[i] * BL[i] for i in range(dim_A))
    
    sum_sq = sum(approx_number_field_element(c)^2 for c in a_embed)
    return sqrt(sum_sq).n()


def approx_det(L):
    """
    Approximates the determinant of a matrix with number field entries.

    INPUT:
        - L: A square matrix with entries in a number field E.

    OUTPUT:
        - The determinant of the matrix of real-valued approximations of L.

    ALGORITHM:
        1. Create a new matrix by replacing each entry of L with its real
           approximation using `approx_number_field_element`.
        2. Compute the determinant of this new matrix of real numbers.
    """
    L_approx = [[approx_number_field_element(x).n() for x in row] for row in L]
    M = Matrix(RR, L_approx)
    return M.determinant()


def approx_LLL_transform(BL):
    """
    Computes the LLL transformation matrix for a rational approximation of BL.

    INPUT:
        - BL: An N x N matrix with entries in a number field E.

    OUTPUT:
        - T: An N x N integer matrix with determinant +-1 such that T*M is
             the LLL-reduced form of M, where M is the rational approximation of BL.

    ALGORITHM:
        1. Get the dimensions of the matrix BL.
        2. Create a rational approximation M of the matrix BL using the
           `approx_number_field_element` function for each entry.
        3. Apply Sage's LLL algorithm to the rational matrix M, requesting that
           the transformation matrix T be returned.
        4. Return T.
    """
    N = BL.nrows()
    E = BL[0, 0].parent()

    BL_approx_list = [[approx_number_field_element(BL[i, j]) for j in range(N)] for i in range(N)]
    M = Matrix(QQ, BL_approx_list)

    _, T = M.LLL(transformation=True)
    return T


# =========================================================================================
#
#   SECTION 4: FINDING A RATIONAL RANK-ONE ELEMENT
#
# =========================================================================================

# This is an absolute constant from the literature for the case n=4 (N=16).
c16 = 2^116 * 3^16

def reduce(A, BL):
    """
    Applies an LLL-based transformation to the basis of an algebra.

    INPUT:
        - A: A Q-algebra isomorphic to Mn(Q) with a given basis.
        - BL: A matrix representing an embedding A -> Mn(E), where the rows
              correspond to the images of the basis elements of A.

    OUTPUT:
        - BA: A new basis for the algebra A. This basis is obtained by applying an
              integer transformation matrix T (from LLL) to the original basis.
              It is expected to be "more reduced".
    """
    N = A.dimension()
    T = approx_LLL_transform(BL)
    BA0 = A.basis()
    BA = [A.sum(T[i, j] * BA0[j] for j in range(N)) for i in range(N)]
    return BA



def is_reduced(BA, BL):
    """
    Checks if a basis satisfies the "reduced" condition from Ivanyos & Rónyai.

    INPUT:
        - BA: A basis of the algebra A.
        - BL: A matrix representing an embedding A -> Mn(E).

    OUTPUT:
        - boolean: True if the product of the norms of the basis images is less
                   than or equal to cN * |det(BL)|. False otherwise. This check
                   is specifically for n=4 (N=16) using the constant c16.

    ALGORITHM:
        1. Calculate the product of the Frobenius norms of the images of each
           element in the basis BA. This is the Left-Hand Side (LHS).
        2. Calculate the Right-Hand Side (RHS) as c16 * |det(BL_approx)|.
        3. Print both LHS and RHS for inspection.
        4. Return true if LHS <= RHS.
    """
    A = BA[0].parent()
    LHS = product(frob_norm(a, A, BL) for a in BA)
    RHS = c16 * abs(approx_det(BL))
    print(f"Reduction Check: LHS (Product of Norms) = {LHS.n(digits=10)}")
    print(f"Reduction Check: RHS (c16 * |det(BL)|) = {RHS.n(digits=10)}")
    return LHS <= RHS




import heapq

def rank_one_in_max_order(A, max_tests=10000):
    """
    Searches for a rank-one element in a maximal order of A.

    INPUT:
        - A: A Q-algebra isomorphic to Mn(Q). It is assumed that the Z-span of
             the input basis of A forms a maximal order.
        - max_tests: The maximum number of candidate elements to test.

    OUTPUT:
        - a: An element in A which has rank one when viewed in Mn(Q).
        - Returns None if no such element is found within the iteration limit.

    ALGORITHM (Circular Search):
        1. Compute an embedding A -> Mn(E) and find a reduced basis BA.
        2. The theory suggests a rank-one element can be found as a short integer
           linear combination of the reduced basis elements.
        3. We perform a deterministic search for x = Σ x_i * b_i by testing
           coefficient vectors (x_0, ..., x_{N-1}) in increasing order of their
           squared L2 norm (x_0^2 + ... + x_{N-1}^2).
        4. A priority queue (min-heap) is used to generate coefficient vectors
           in this "circular" order, starting from (0,...,0) and expanding outwards.
        5. For each generated element `x`, check if it is rank-one. If so, return it.
        6. The search stops after `max_tests` candidates have been checked.
    """
    N = A.dimension()
    n = int(sqrt(N))

    print("Step 1: Computing embedding A -> Mn(E)...")
    BL = embed_R_N(A)

    print("Step 2: Computing reduced basis using LLL...")
    BA = reduce(A, BL)

    print("Step 3: Checking if the new basis is reduced...")
    if not is_reduced(BA, BL):
        print("Warning: The basis does not satisfy the reduction condition. The search may be less efficient.")
    else:
        print("Basis is reduced. Proceeding with search.")

    print(f"Step 4: Starting circular search for rank-one element (up to {max_tests} tests).")
    
    # Priority queue stores tuples of (squared_norm, coeffs_tuple)
    pq = [(0, tuple([0] * N))]
    # Set to keep track of coefficient vectors we have already visited
    visited = {tuple([0] * N)}
    
    test_count = 0
    while pq and test_count < max_tests:
        squared_norm, coeffs_tuple = heapq.heappop(pq)
        
        # Test the element (skip the zero vector, which is always first)
        if squared_norm > 0:
            test_count += 1
            coeffs = list(coeffs_tuple)
            x = A.sum(coeffs[i] * BA[i] for i in range(N))
            
            print(f"Test {test_count}: coeffs={coeffs} -> norm={frob_norm(x,A,BL)}")

            # The rank check `right_rank(x, A) == n` is based on the user's provided code.
            if right_rank(x, A) == n:
                print(f"\nSuccess: Found rank-one element after {test_count} tests.")
                print(f"Element: {x}")
                print(f"Coefficients: {coeffs}")
                return x

        # Add all neighbors of the current vector to the queue
        for i in range(N):
            base_coeffs = list(coeffs_tuple)
            
            # Neighbor by incrementing coordinate i
            base_coeffs[i] += 1
            neighbor_p_tuple = tuple(base_coeffs)
            if neighbor_p_tuple not in visited:
                visited.add(neighbor_p_tuple)
                new_norm_sq = sum(c*c for c in neighbor_p_tuple)
                heapq.heappush(pq, (new_norm_sq, neighbor_p_tuple))
            
            # Neighbor by decrementing coordinate i
            base_coeffs[i] -= 2 # From +1 to -1
            neighbor_m_tuple = tuple(base_coeffs)
            if neighbor_m_tuple not in visited:
                visited.add(neighbor_m_tuple)
                new_norm_sq = sum(c*c for c in neighbor_m_tuple)
                heapq.heappush(pq, (new_norm_sq, neighbor_m_tuple))

    print(f"\nSearch failed: No rank-one element found within {max_tests} tests.")
    return None
