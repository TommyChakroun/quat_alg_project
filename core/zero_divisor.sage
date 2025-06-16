load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("core/rank_one.sage")
load("maximal_orders/maximal_orders_utilities.sage")
load("maximal_orders/find_maximal_orders.sage")
load("explicit_embeding_real_number/explicit_embedding.sage")


#------------------------------------------------------------------------------------------
#
#            ZERO DIVISOR IN ALGEBRA ISOMORPHIC TO Mn(Q)
#
#------------------------------------------------------------------------------------------

## Utilities for approximation of real number or p adic number 

def approx_number_field_element(E,x, precision=1e-10):
    """
    Approximates a list of matrices with entries in a real number field by matrices with rational entries.

    INPUT:
        - E -- a real number field, defined as a quotient of Q[x] by an irreducible polynomial g.
        - Zbasis_J -- a list of matrices with entries in E.
        - precision -- a positive real number specifying the desired approximation accuracy.

    OUTPUT:
        - Zbasis_J_approx -- a list of matrices with rational entries that approximate the input matrices.
    """

    # Get the defining polynomial of the number field
    g = E.defining_polynomial()

    # Find the real roots of the polynomial g
    real_roots = g.real_roots()
    
    if not real_roots:
        raise ValueError("The defining polynomial of E has no real roots.")
        
    # Choose one of the real roots for the embedding
    gamma_real = real_roots[0]

    poly_rep = x.polynomial()
            
    # Substitute the abstract root gamma with its real numerical value
    real_val = poly_rep(gamma_real)
                
    # Approximate the real number by a rational number
    return (QQ(real_val.n(digits=int(-log(precision, 10)))))
        


## Zero divisor in a finite field or a number field ( the zero divisor has probably rank 1)

def zero_divisor_finite_matrix_ring(A):
    """
    INPUT : 
        -- A -- an algebra over a finite field Fp isomorphic to Mn(Fp)
    OUTPUT :
        -- x -- in A a zero divisor
         
    ALGORITHM :
        - Lajos Ronyai, Simple Algebras Are Difficult (1987) Section 5

    """
    BA = list(A.basis())
    dim_A = dimension(A)
    n = int(sqrt(dim_A))
    F = A.base_ring()
    p = F.order()

    n_ite = 100

    for ite in range(n_ite):
        c = sum(F.random_element()*a for a in BA)
        pi = minimal_polynomial(c,A)
        if not pi.is_irreducible():
            fact = [g for g,_ in pi.factor()]
            g = fact[0]
            return g(c)
    return "Not found"


def zero_divisor_real_matrix_ring(A):
    """
    INPUT : 
        -- A -- a Q-algebra given by structure constant assumed to be isomorphic to Mn(Q)
    OUTPUT :
        -- embedding_dict -- A dictionnary with key 0,.., N-1 where N = dim_A and value matrices in Mn(E)
                            for E some number field E = Q(gamma) representing an emmeding i : A -> Mn(E)
                            which map the i th element of the basis of A onto the matrix embedding_dict[i]
    """

    N = dimension(A)
    BA = list(A.basis())
    n = int(sqrt(N))
    H = 2*n*(n-1)

    n_ite = 100

    for ite in range(100):
        t = [randint(1,H) for l in range(N)]
        a = sum( t[i]*BA[i] for i in range(N))
        pi = minimal_polynomial(a,A)
        if pi.degree()==n and gcd(pi,pi.derivative())==1:
            fact = [g for g,_ in pi.factor()]
            for g in fact:
                if len(g.real_roots()) >0:
                    E = NumberField(g, 'gamma')
                    gamma = E.gen()

                    R = PolynomialRing(E, 't')
                    t = R.gen()

                    B = FiniteDimensionalAlgebra(E,A.table()) # B = A⊗E over Q

                    pi_E = R(pi)
                    Q = pi_E // (t-gamma)

                    a_in_B = sum(c*b for c,b in zip(get_coefficients(a,A),B.basis()))
                    
                    b = Q(a_in_B)

                    return B,b
    return "Not found"


## Element e of ranke 1 and such that e^2=e form zero divisor.

def solve_xax_eq_x(A, x):
    """
    Optimized version to solve x * a * x = x in algebra A given structure constants.
    Assumes A is defined over QQ and has basis e_1,...,e_n.
    """

    # Basis and structure constants
    BA = list(A.basis())
    n = len(BA)
    F = A.base_ring()
    C = structure_constants(A,BA)  # table: C[j][i][k] = c_{i j k}

    # lambda_i from x = sum lambda_i e_i
    lambda_vec = vector(F, get_coefficients(x,A))

    # Precompute v^(k)[s] = sum_i lambda_i * c_{i k s}
    V = [vector(F, [
        sum(lambda_vec[i] * C[k][i][s] for i in range(n))
        for s in range(n)
    ]) for k in range(n)]

    # Precompute w^(t)[s] = sum_j lambda_j * c_{s j t}
    W = [vector(F, [
        sum(lambda_vec[j] * C[j][s][t] for j in range(n))
        for s in range(n)
    ]) for t in range(n)]

    # Build M[t,k] = dot_product(V[k], W[t])
    M = Matrix(F, n, n, lambda t, k: V[k].dot_product(W[t]))


    # Solve M * mu = lambda
    mu = M.solve_right(lambda_vec)

    # Return a = sum mu_i * e_i
    return sum(mu[i] * BA[i] for i in range(n))


def rank_one_from_zero_divisor(A):
    """
    INPUT : 
        -- A -- an algebra over Q isomorphic to Mn(Q)
    OUTPUT :
        -- x -- in A an element of "rank one", that is dim(Ax)=n.
    """
    dim_A = dimension(A)
    n = int(sqrt(dim_A))
    BA = list(A.basis())
    F = A.base_ring()

    x = zero_divisor_finite_matrix_ring(A)

    # write pi_x = t Q(t) and let y = Q(x)
    # we have xy=0 and y !=0
    # if rank(x)>n/2, dim ker x <= n/2 and since Im y in Ker x, rank(y) < n/2.
    
    pi = minimal_polynomial(x,A)
    R = PolynomialRing(F,'x')
    t = R.gen()
    Q = pi//t
    if right_rank(A,x)>n*n/2:
        x = Q(x)
    if right_rank(A,x)==n:
        return x


    a = solve_xax_eq_x(A,x)
    e = a*x
    basis_eAe = extract_basis_from_generators(A,[e*BA[i]*e for i in range(dim_A)])
    structure_constants_eAe = structure_constants_subspace(A, basis_eAe)

    B = FiniteDimensionalAlgebra(F,structure_constants_eAe)

    z = rank_one_from_zero_divisor(B)

    z_lift_in_A = sum(c*b for c,b in zip(get_coefficients(z,A),basis_eAe))
    
    return z_lift_in_A



## Extension to the ratitional matrix ring


def zero_divisor_matrix_ring_mod_pn(A,BO = None,p=5,n=4):
    """
    INPUT : 
        -- A -- an algebra over Q isomorphic to Mn(Q)
    OUTPUT :
        -- x -- in A a zero divisor
    
    ALGORITHM :
    1. Case n=2 : 
        - J. E. Cremona and D. Rusin, Efficient solution of rational conics, Math. Comp. 72 (2003)
        -Gabor Ivanyos and Agnes Szanto, Lattice basis reduction for indefinite forms and an application, Proceedings of the 5th Conference on Formal Power Series and Algebraic Combinatorics (Florence, 1993), Discrete Math
        - Denis Simon, Solving quadratic equations using reduced unimodular quadratic forms, Math.Comp. 74 (2005)
    2. Case n=4:
        - Jana Pílniková, Trivializing a central simple algebra of degree 4 over the rational numbers (2007)
    
    3. Potential idea in general : 
        start computing O = Ze1⊕... ⊕ZeN in A. We know O isom Mn(Z) as a ring. So for all prime p 
        Op :=  F_p e1⊕... ⊕F_p eN isom Mn(F_p)
        try to compute an element of rank 1 in Op on the form s = a1 e1 + ..+ aN eN 
        lift each ai in ni Z and hope the lift of s, that is S = n1 e1 + ..+ nN eN in O
        is a zero divisor in A.
        So far we know that det(S)=det(s)=0 mod p so p|det(S). [determinant is invariander automorphism of Mn(F)]
        But we want det(S) = 0

    """
    dim_A = dimension(A) # = 16 in for our purpose
    BA = A.basis()

    if BO == None:
        BO = max_order(A)
    

    Op,_ = finite_algebra_from_order(A,BO,p)
    BOp = Op.basis()

    x = rank_one_from_zero_divisor(Op)
    temp = solve_xax_eq_x(Op, x)
    e = temp*x # e satisfy rank e =1 view in M4(F_p) and e^2=e

    coords_bar = get_coefficients(e,Op) # natural coefficient of e in Op


    a0 = A.sum(Integer(c) * b for c, b in zip(coords_bar, BO))

    ak = a0
    k_iterations = n
    
    for k in range(k_iterations):
        ak_squared = ak**2
        ak_cubed = ak_squared * ak
        ak = 3 * ak_squared - 2 * ak_cubed
        ak_coeffs_int = coordinate(ak**2-ak,A,BO)

    M = p**(2**(k_iterations+4))
    

    return ak
    

def zero_divisor_matrix_ring_aprox_from_real(A):
    dim_A = dimension(A)
    n = int(sqrt(dim_A))
    BA = A.basis()

    B,b = zero_divisor_real_matrix_ring(A)
    e = solve_xax_eq_x(B, b)*b
    print("e : ")
    print(e)
    print("CHECK e^2 =e in A extended to Q(gamma) : ")
    print(e**2==e)
    print("CHECK e has rank 1 in  A extended to Q(gamma) : ")
    print(right_rank(B,e) == n)
    real_coords  = get_coefficients(e,B)
    E = B.base_ring()
    approx_coords = [approx_number_field_element(E,x) for x in real_coords]

    return sum(approx_coords[i]*BA[i] for i in range(dim_A))



def zero_divisor_matrix_ring(A, BO, p=5, n=4):
    """
    Trouve un zéro-diviseur dans A en utilisant une approximation p-adique 
    et la réduction de réseau (LLL).
    """
    # 1. Obtenir un idempotent approximatif e tel que e^2 - e = 0 (mod p^k)
    #    (La fonction auxiliaire est renommée pour plus de clarté)
    k = p**(2**n)
    e = zero_divisor_matrix_ring_mod_pn(A, BO, p, n)

    print("e : ")
    print(e)
    print("e^2-e ")
    print(e**2-e)

    # 2. Construire les générateurs du réseau J = e*BO + k*BO
    #    Les lignes de la matrice M sont les coordonnées de ces générateurs.
    rows = []
    for basis_vec in BO:
        print("e*b for b in BO")
        rows.append(coordinate(e * basis_vec, A, BO))
    for basis_vec in BO:
        print("p^(2^n)*b for b in BO")
        rows.append(coordinate(k * basis_vec, A, BO))
    
    M = Matrix(ZZ, rows)

    print("M : ")
    print(M)


    print("COmpute LLL ... ")
    L_base_reduite = M.LLL()
    print("LLL done.")
    print( L_base_reduite)


    
    one_A = A.one() 

    for b_coords in L_base_reduite:
        f_candidat = sum(b_coords[i]*BO[i] for i in range(len(BO)))
        if f_candidat.is_zero():
            continue
        print("f**2-f")
        if f_candidat**2 == f_candidat and f_candidat != one_A:
            print("Find non trivial idempotent via LLL.")
            return f_candidat

    print("Fail to find non trivial idempotent")
    return None






























## TESTS

def test_17(n, p, P, target_x=None):
    if target_x==None:
        target_x = Matrix(GF(p),[[1 for _ in range(n)]for _ in range(n)])

    A = MatrixSpace(QQ, n, n)
    BA = list(A.basis())
    new_BA = [sum(P[i][k] * BA[k] for k in range(n**2)) for i in range(n**2)]
    D = MatrixSpace(GF(p), n, n)
    new_BA_mod_p = [Matrix(GF(p), M) for M in new_BA]
    x = target_x
    coords = coordinate(x, D, new_BA_mod_p)
    lift_x = sum(Integer(coords[i]) * new_BA[i] for i in range(n**2))
    det_A = lift_x.det()
    rank_A = lift_x.rank()
    return lift_x


def rank_one_matrix(p,v1,v2):
    V1 = Matrix(GF(p),4,1,v1)
    V2 = Matrix(GF(p),1,4,v2)
    return V1*V2