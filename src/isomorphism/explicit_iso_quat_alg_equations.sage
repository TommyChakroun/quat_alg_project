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



def quaternionic_complement(B, mu):
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

    (i, j, k), (a, b) = B.gens(), B.invariants()
    x0,x1,x2,x3 = mu.coefficient_tuple()


    if x1 != 0:
        v1 = (b*x2,-a*x1,0)
        v2 = (b*x3,0,x1)
    elif x2 != 0:
        v1 = (b*x2,-a*x1,0)
        v2 = (0,a*x3,x2)
    elif x3!=0 :
        v1 = (1,0,0)
        v2 = (0,1,0)

    M = Matrix(QQ,[v1,v2])
    N = M.LLL()

    return N[0][0]*i + N[0][1]*j + N[0][2]*k



def complement_with_prime_norm(B,mu):
    (i, j, k), (a, b) = B.gens(), B.invariants()
    N = quaternionic_complement_basis(B,mu)

    nu1 = N[0][0]*i + N[0][1]*j + N[0][2]*k
    nu2 = N[1][0]*i + N[1][1]*j + N[1][2]*k

    nb_ite = int(10000*log(max(abs(nu1.reduced_norm()),abs(nu2.reduced_norm()))))
 
    for ite in range(nb_ite):
        nu = randint(-50,50)*nu1+randint(-50,50)*nu2
        gamma = abs(nu.reduced_norm())
        if is_pseudoprime(gamma):
            return nu

    return nu1





#------------------------------------------------------------------------------------------
#           GENERAL EXPLICIT ISOMORPHISM FOR QUATERNION ALGEBRAS
#------------------------------------------------------------------------------------------


## Solving a quadratic equation of dimension 4 over Q : aX^2+bY^2+cZ^2+dW^2 = 0
## is more efficient than solving a quadratic equation of dimension 3 over a quadratic numver field : a U^2 + bV^2+ c W^2=0


def diagonal_qfsolve(a,factors = None):
    """
    Solves x^T G x = 0 for a diagonal integer matrix G = diag(a),
    optimized by precomputing the factorization of the determinant of G.
    INPUT :
        -- a = [a1,..,an] -- a list of integers
        -- factors = [F1,..,Fn] -- the list of their factorization if given
    OUTPUT :
        -- x1,..,xn -- a solution to x^T G x =0
    """
    if factors == None :
        factors = [factor(abs(e)) for e in a]
    
    fact_det = product(f for f in factors)
    combined_factors =dict(fact_det)

    # Format the factors and matrix into PARI strings.
    pari_F_str = f"[{'; '.join(f'{p},{e}' for p, e in sorted(combined_factors.items()))}]"

    # Call PARI and return the raw result.
    sol = pari(f"qfsolve([matdiagonal({a}), {pari_F_str}])")
    try :
        return [QQ(e) for e in sol]
    except:
        print("Pari was unable to solve the digonal quadratic :")
        print(a)
        print(f"sol : {sol}")



def iso_quat_alg(A, B):
    """
    Finds an explicit isomorphism between two quaternion algebras A and B.

    ALGORITHM:
    1. Handle the case where A and B are matrix rings.
    2. To embed A = (alpha,beta | K) into B = (a,b | K), we first find an element
       `mu` in B such that `mu^2 = alpha`. This is the image of `i_A`. This step
       requires finding a rational point on the conic `a*Y^2 + b*Z^2 - a*b*W^2 + alpha = 0`.
    3. We then find an element `nu_0` in B that anticommutes with `mu`.
    4. We find the image of `j_A`, which will be of the form `j' = (x + y*mu)*nu_0`.
       We solve for `x, y` such that `(j')^2 = beta`. This reduces to solving the
       norm equation `x^2 - alpha*y^2 = beta / nu_0^2`, which is done by finding a point
       on the conic `X^2 - alpha*Y^2 - (beta/nu_0^2)*Z^2 = 0`.
    5. The isomorphism is then defined by `i_A -> mu` and `j_A -> j'`.

    INPUT:
        - ``A``, ``B`` -- Two quaternion algebras over the same number field K.

    OUTPUT:
        - A tuple ``(is_isomorphic, isom_map)`` where:
            - ``is_isomorphic`` is `True` if an isomorphism exists, `False` otherwise.
            - ``isom_map`` is a list of four elements in B which are the images
              of the basis elements 1, i, j, k of A.
    """
    if A.base_ring()!=QQ or B.base_ring()!=QQ:
        raise ValueError("Algebras must be over rational field")
    ram_A = A.ramified_primes()
    ram_B = B.ramified_primes()
    if ram_A != ram_B:
        return False,[]
    if ram_A ==[]:
        raise ValueError("A,B are split, use the other function")

    alpha, beta = A.invariants()
    a, b = B.invariants()
    i_A, j_A, k_A = A.gens()
    i_B, j_B, k_B = B.gens()

    
    # We will need all this factorization
    Falpha = factor(abs(alpha))
    Fbeta = factor(abs(beta))
    Fa = factor(abs(a))
    Fb = factor(abs(b))


    # Step 2: Find mu in B such that mu^2 = alpha
    # We have to solve a x^2+by^2-abz^2 = alpha
    # we solve the homogenous  a x^2+by^2-abz^2 -alpha w^2
    sol = diagonal_qfsolve([a,b,-a*b,-alpha],factors =[Fa,Fb,Fa*Fb,Falpha])
    mu = sol[0]/sol[3]*i_B+sol[1]/sol[3]*j_B+sol[2]/sol[3]*k_B

    # Step 3: Find an element nu that anticommutes with mu
    nu0 = quaternionic_complement(B, mu)
    gamma0 = -nu0.reduced_norm() 
    nu = gamma0.denominator()*nu0
    gamma = -nu.reduced_norm() 

    # We will need this other factorization

    Fgamma = factor(gamma)

    # Step 4: Find j' = (x + y*mu)*nu_0 such that (j')^2 = beta
    sol2 = diagonal_qfsolve([gamma,-alpha*gamma,-beta],factors = [Fgamma,Falpha*Fgamma,Fbeta])
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
    nu = quaternionic_complement(B, mu)
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
