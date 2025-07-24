
def my_qfsolve(S, fact_det=None):
    """
    Solves x^T S x = 0 using PARI's qfsolve.
    This is optimized by allowing for a pre-computed factorization of the determinant of S.

    INPUT:
        - S: A symmetric matrix with integer coefficients.
        - fact_det: (Optional) The pre-computed factorization of abs(det(S)).

    OUTPUT:
        - A vector x that is a non-trivial solution to x^T S x = 0, or None if no such solution exists.
    """
    n = S.nrows()
    if S.determinant() == 0:
        return ker = S.kernel().basis()[0]
        
    if fact_det is None:
        fact_det = factor(abs(S.determinant()))
    
    # Format the factors and matrix into PARI strings.
    pari_F_str = f"[{'; '.join(f'{p},{e}' for p, e in sorted(fact_det))}]"
    pari_mat_str = ";".join(str(row) for row in S.rows())

    # Call PARI. qfsolve returns 0 if there is no solution.
    # The algorithm relies on local information at primes dividing the determinant[cite: 12].
    try:
        sol = pari(f"qfsolve([{pari_mat_str}, {pari_F_str}])")
        if sol == -1:
            return None
        # Convert the PARI column vector to a SageMath vector of rationals.
        return Vector(QQ, sol.list())
    except Exception as e:
        print(f"Error calling PARI: {e}")
        return None





def isotropic_subspace(S, fact_det=None):
    """
    Finds a basis for a maximal totally isotropic subspace for the quadratic form defined by S.
    This implements the recursive strategy from D. Simon's paper.

    INPUT:
        - S: A symmetric matrix with integer coefficients.
        - fact_det: (Optional) Pre-computed factorization of abs(det(S)).

    OUTPUT:
        - A list [u1,...,ur] of vectors forming a basis for a maximal isotropic subspace.
          Returns an empty list [] if the form is anisotropic.
    """

    if fact_det is None:
        if S.determinant() == 0:
            raise ValueError("The input matrix must be non-degenerate.")
        fact_det = factor(abs(S.determinant()))

    # --- Step 1: Find a single isotropic vector (Base Case) ---
    # Use the algorithm from Simon's paper (implemented in PARI's qfsolve)[cite: 56].
    u = my_qfsolve(S, fact_det)
    
    # If no vector is found, the form is anisotropic. This is the base case for the recursion.
    if u is None:
        return []

    # --- Step 2: Construct a Hyperbolic Plane H = span{u, v} ---
    # We need to find a vector v such that u^T*S*v = 1 and v^T*S*v = 0.
    
    # First, find any vector w such that u^T*S*w is not zero.
    # Since u is not in the kernel of S, u^T*S is not the zero vector.
    w = None
    for i in range(S.nrows()):
        # Check if the i-th component of the row vector u^T*S is non-zero
        if (u.transpose() * S)[0, i] != 0:
            w = Vector(QQ, S.nrows())
            w[i] = 1
            break
            
    # Scale w to get a temporary v with u^T*S*v_temp = 1
    scaling_factor = (u.transpose() * S * w)[0, 0]
    v_temp = w / scaling_factor
    
    # Adjust v_temp to make it isotropic, creating the final v.
    v_quadratic_form = (v_temp.transpose() * S * v_temp)[0, 0]
    v = v_temp - (v_quadratic_form / 2) * u

    # --- Step 3: Restrict the form to the orthogonal complement of H ---
    # The complement W = H^perp is the space of all vectors x where u^T*S*x = 0 and v^T*S*x = 0.
    # We find a basis for this space by computing the kernel of the map defined by these two conditions.
    ortho_map = Matrix([u.transpose() * S, v.transpose() * S])
    W_basis = ortho_map.right_kernel().basis()

    # If the complement is empty (dim n=2), our basis is just u.
    if not W_basis:
        return [u]
        
    # Construct the matrix S' for the quadratic form restricted to the basis of W.
    S_prime_rows = []
    for w1 in W_basis:
        new_row = [(w1.transpose() * S * w2)[0, 0] for w2 in W_basis]
        S_prime_rows.append(new_row)
    S_prime = Matrix(QQ, S_prime_rows)

    # --- Step 4: Recursive call on the smaller space ---
    # The prime factors of det(S') are a subset of those of det(S), so we reuse fact_det.
    RecursiveBasis = isotropic_subspace(S_prime, fact_det)

    # --- Step 5: Lift the result back to the original space ---
    # The vectors in RecursiveBasis are coordinates wrt W_basis. Convert them back.
    ChangeOfBasisMatrix = Matrix(W_basis).transpose()
    LiftedBasis = [ChangeOfBasisMatrix * b for b in RecursiveBasis]

    # The new maximal basis is the vector u from this step plus the lifted basis from the recursion.
    return [u] + LiftedBasis