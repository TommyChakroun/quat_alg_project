# SageMath requires these external files for utility functions and type definitions.
load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/isomorphism/explicit_iso_matrix_ring.sage")

#------------------------------------------------------------------------------------------
#
#           EXPLICIT ISOMORPHISM PROBLEM FOR QUATERNION ALGEBRAS
#               BY SOLVING ALGEBRAIC EQUATIONS
#
#------------------------------------------------------------------------------------------

# This file provides algorithms to solve the explicit isomorphism problem for
# quaternion algebras over number fields. The core of the problem relies on the
# ability to find rational points on ternary conics over a number field K or Q.
# Specifically, we need to solve equations of the form:
#
#   a*X^2 + b*Y^2 + c*Z^2 = 0, for X,Y,Z in K
#
# The problem of finding a rational point on a conic is a major topic in
# computational number theory. The general strategy involves:
# 1. Checking for trivial solutions (e.g., if -a/b is a square in K).
# 2. Reducing the conic to a norm equation of the form N_{L/K}(z) = d, where L/K
#    is a quadratic extension.
#
# In this script, we assume that such a "conic solver" is available and we build
# the isomorphism algorithms upon it.


# In general if A =(a,b|K), B =(c,d |K) are two quaternion algebra over K, to sove the explicit isomorphism problem 
# we need to solve 2 ternary conics over K:
# 1. The first conic is used to find an element mu in B such that
#    mu^2 = a, where a is the invariant of A. This equation is solved over quadratic extension of K.
# 2. The second conic is used to find an element nu in B that anticommutes with mu. This equation is solved over K.

# However if A = (a,b|K) isom matrix ring (1,1 |K) isom M_2(K), then we need only 1 ternary conic over K:
# 1. The conic is used to find a non-zero element e in A such e^2=0. This is solved over K.

# This is why we sepparate the split case and the general case.












#------------------------------------------------------------------------------------------
#           EXPLICIT ISOMORPHISM PROBLEM OF SPLIT QUATERNION ALGEBRAS
#------------------------------------------------------------------------------------------



def matrix_ring(A, e=None):
    """
    Finds an explicit isomorphism from a quaternion algebra A to the 2x2 matrix ring M_2(K).

    The algorithm relies on finding a non-zero element e in A such that e^2 = 0.
    The algebra A then acts on the 2-dimensional left ideal I = Ae, which gives an
    isomorphism A -> End_K(I) ~= M_2(K).

    INPUT:
        - ``A`` -- A quaternion algebra over a number field K.
        - ``e`` -- (Optional) A non-zero element of A such that e^2 = 0. If not
                   provided, the function will attempt to find one.

    OUTPUT:
        - A tuple ``(is_isomorphic, matrices)`` where:
            - ``is_isomorphic`` is `True` if A is a matrix ring, `False` otherwise.
            - ``matrices`` is a list of four 2x2 matrices representing the images
              of the basis elements 1, i, j, k under the isomorphism.
    """
    if not A.is_matrix_ring():
        return False, []

    K = A.base_ring()
    i, j, k = A.gens()
    a, b = A.invariants()

    if e is None:
        X, Y, Z = K['X, Y, Z'].gens()
        C = Conic(a*X^2 + b*Y^2 - a*b*Z^2)
        has_sol, sol = C.has_rational_point(point=True)
        e = sol[0]*i + sol[1]*j + sol[2]*k

    BAe = [e]
    for s in [i,j,k]:
        x = s*e
        if x.reduced_trace() != 0:
            BAe.append(x)
            break
    
    matrices = []
    BA = A.basis()
    for a_elem in BA:
        rows = []
        for z in BAe:
            y = a_elem * z 
            coords = coordinate(y, A, BAe)
            rows.append(coords)
        M = Matrix(K, rows).transpose()
        matrices.append(M)

    return True, matrices



#------------------------------------------------------------------------------------------
#           UTILITIE : QUATERNION COMPLEMENT
#------------------------------------------------------------------------------------------



def quadratic_complement(B, mu):
    """
    Finds an element nu in B that anticommutes with a given pure quaternion mu.

    Given mu in B with mu^2 in K, the set of elements nu satisfying
    mu*nu + nu*mu = 0 forms a 2-dimensional K-subspace of the pure quaternions in B.
    This function finds one non-zero element in that space.

    INPUT:
        - ``B`` -- A quaternion algebra (a,b | K).
        - ``mu`` -- A pure quaternion in B (i.e., trace is 0).

    OUTPUT:
        - ``nu`` -- A non-zero pure quaternion in B such that nu*mu = -mu*nu.
    """
    if mu.reduced_trace() != 0:
        raise ValueError("Input element mu must be a pure quaternion.")

    K, (i, j, k), (c, d) = B.base_ring(), B.gens(), B.invariants()
    mu_coords = mu.coefficient_tuple()
    m1, m2, m3 = mu_coords[1], mu_coords[2], mu_coords[3]

    if c*m1 != 0:
        n1, n2, n3 = -d*m2 / (c*m1), 1, 0
    elif d*m2 != 0:
        n1, n2, n3 = 1, -c*m1 / (d*m2), 0
    elif c*d*m3 != 0:
        n1, n2, n3 = 1, 0, c*m1 / (c*d*m3)
    else:
        return j

    return n1*i + n2*j + n3*k









#------------------------------------------------------------------------------------------
#           GENERAL EXPLICIT ISOMORPHISM FOR QUATERNION ALGEBRAS
#------------------------------------------------------------------------------------------


## Solving a quadratic equation of dimension 4 over Q : aX^2+bY^2+cZ^2+dW^2 = 0
## is more efficient than solving a quadratic equation of dimension 3 over a quadratic numver field : a U^2 + bV^2+ c W^2=0


def iso_quat_alg(A, B):
    """
    Finds an explicit isomorphism between two quaternion algebras A and B.

    ALGORITHM:
    1. Handle the case where A and B are matrix rings.
    2. To embed A = (a,b | K) into B = (c,d | K), we first find an element
       `mu` in B such that `mu^2 = a`. This is the image of `i_A`. This step
       requires finding a rational point on the conic `c*Y^2 + d*Z^2 - c*d*W^2 + a = 0`.
    3. We then find an element `nu_0` in B that anticommutes with `mu`.
    4. We find the image of `j_A`, which will be of the form `j' = (x + y*mu)*nu_0`.
       We solve for `x, y` such that `(j')^2 = b`. This reduces to solving the
       norm equation `x^2 - a*y^2 = b / nu_0^2`, which is done by finding a point
       on the conic `X^2 - a*Y^2 - (b/nu_0^2)*Z^2 = 0`.
    5. The isomorphism is then defined by `i_A -> mu` and `j_A -> j'`.

    INPUT:
        - ``A``, ``B`` -- Two quaternion algebras over the same number field K.

    OUTPUT:
        - A tuple ``(is_isomorphic, isom_map)`` where:
            - ``is_isomorphic`` is `True` if an isomorphism exists, `False` otherwise.
            - ``isom_map`` is a list of four elements in B which are the images
              of the basis elements 1, i, j, k of A.
    """
    if A.base_ring() != B.base_ring():
        raise ValueError("Algebras must be over the same base field.")
    if not A.is_isomorphic(B):
        return False, []

    K = A.base_ring()
    a, b = A.invariants()
    c, d = B.invariants()
    i_A, j_A, k_A = A.gens()
    i_B, j_B, k_B = B.gens()

    if A.is_matrix_ring():
        # This case requires a more complex implementation of the inverse map
        # from M_2(K) to B. We focus on the non-split case here.
        print("Did'n implement for matrix ring case yet.")
        return False, [] # Placeholder for matrix ring implementation


    # Step 2: Find mu in B such that mu^2 = a
    # We have to solve c x^2+dy^2-cdz^2 = a
    # we solve the homogenous  c x^2+dy^2-cdz^2 -a w^2
    from sage.quadratic_forms.qfsolve import qfsolve
    D = diagonal_matrix(QQ,[c,d,-c*d,-a])
    sol = qfsolve(D)
    mu = sol[0]/sol[3]*i_B+sol[1]/sol[3]*j_B+sol[2]/sol[3]*k_B


    # Step 3: Find an element nu that anticommutes with mu
    nu = quadratic_complement(B, mu)
    beta = -nu.reduced_norm()

    # Step 3: Find j' = (x + y*mu)*nu_0 such that (j')^2 = b
    D = diagonal_matrix(QQ,[1,-a,-b/beta])
    sol2 = qfsolve(D)
    x = sol2[0] / sol2[2]
    y = sol2[1] / sol2[2]

    # Step 5: Construct the isomorphism map
    isom_i = mu
    isom_j = (x + y*mu) * nu
    isom_map = [B.one(), isom_i, isom_j, isom_i * isom_j]

    return True, isom_map




#------------------------------------------------------------------------------------------
#           GENERAL EXPLICIT ISOMORPHISM FOR QUATERNION ALGEBRAS
#------------------------------------------------------------------------------------------




def iso_quat_alg_old(A, B):
    """
    Finds an explicit isomorphism between two quaternion algebras A and B.

    ALGORITHM:
    1. Handle the case where A and B are matrix rings.
    2. To embed A = (a,b | K) into B = (c,d | K), we first find an element
       `mu` in B such that `mu^2 = a`. This is the image of `i_A`. This step
       requires finding a rational point on the conic `c*Y^2 + d*Z^2 - c*d*W^2 + a = 0`.
    3. We then find an element `nu_0` in B that anticommutes with `mu`.
    4. We find the image of `j_A`, which will be of the form `j' = (x + y*mu)*nu_0`.
       We solve for `x, y` such that `(j')^2 = b`. This reduces to solving the
       norm equation `x^2 - a*y^2 = b / nu_0^2`, which is done by finding a point
       on the conic `X^2 - a*Y^2 - (b/nu_0^2)*Z^2 = 0`.
    5. The isomorphism is then defined by `i_A -> mu` and `j_A -> j'`.

    INPUT:
        - ``A``, ``B`` -- Two quaternion algebras over the same number field K.

    OUTPUT:
        - A tuple ``(is_isomorphic, isom_map)`` where:
            - ``is_isomorphic`` is `True` if an isomorphism exists, `False` otherwise.
            - ``isom_map`` is a list of four elements in B which are the images
              of the basis elements 1, i, j, k of A.
    """
    if A.base_ring() != B.base_ring():
        raise ValueError("Algebras must be over the same base field.")
    if not A.is_isomorphic(B):
        return False, []

    K = A.base_ring()
    a, b = A.invariants()
    c, d = B.invariants()
    i_A, j_A, k_A = A.gens()
    i_B, j_B, k_B = B.gens()

    if A.is_matrix_ring():
        # This case requires a more complex implementation of the inverse map
        # from M_2(K) to B. We focus on the non-split case here.
        print("Did'n implement for matrix ring case yet.")
        return False, [] # Placeholder for matrix ring implementation


    # Step 2: Find mu in B such that mu^2 = a
    x = K['x'].gen()
    L = K.extension(x**2 - a, 't')

    is_c_square, r = L(c).is_square(root=True)

    if is_c_square:
        if K.characteristic() == 2:
            raise NotImplementedError("Special case for characteristic 2 not implemented.")
        # Since c=r^2, we solve (U-rV)(U+rV) = d by setting up a
        # simple linear system. This forces a solution where W=1.
        U_prime = (L(d) + 1) / 2
        V_prime = (1 - L(d)) / (2 * r)

    else :
        # General case: we solve the conic c*Y^2 + d*Z^2 - c*d*W^2 + a = 0
        P_L, (U, V, W) = L['U, V, W'].objgens()
        C = Conic(U**2 - c*V**2 - d*W**2)
        has_sol_L, sol_L = C.has_rational_point(point=True)

        U_prime = sol_L[0] / sol_L[2]
        V_prime = sol_L[1] / sol_L[2]
    
    u1,u2 = U_prime.vector()
    v1,v2 = V_prime.vector()

    numerator = u1 + v1*i_B + j_B
    denominator = u2 + v2*i_B 

    mu = numerator/denominator 

    # Step 3: Find an element nu that anticommutes with mu
    nu = quadratic_complement(B, mu)
    beta = -nu.reduced_norm()

    # Step 3: Find j' = (x + y*mu)*nu_0 such that (j')^2 = b
    P2, (X, Y, Z) = K['X, Y, Z'].objgens()
    C = Conic(X**2 - a*Y**2 - b/beta*Z**2)
    has_sol, sol2 = C.has_rational_point(point=True)
    x = sol2[0] / sol2[2]
    y = sol2[1] / sol2[2]

    # Step 5: Construct the isomorphism map
    isom_i = mu
    isom_j = (x + y*mu) * nu
    isom_map = [B.one(), isom_i, isom_j, isom_i * isom_j]

    return True, isom_map
