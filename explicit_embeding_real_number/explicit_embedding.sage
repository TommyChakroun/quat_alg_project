load("utilities/utilities.sage")
load("minimal_ideals/minimal_ideals_manually.sage")
load("core/explicit_iso_matrix_ring.sage")
load("utilities/algebra_type.sage")

#------------------------------------------------------------------------------------------
#
#            EXPLICIT EMBEDING  A ->  A⊗R -> Mn(R)
#
#------------------------------------------------------------------------------------------

"""
Let A be a central simple algebra over Q of dimension N = n^2
Suppose that A is isomrophic to Mn(Q).
Then there exist alpha in R such that if K = Q(alpha) we can compute an explicit 
isomorphism
"""


def real_splitting_element(A,nb_ite = 10):
    """
    see https://www.ams.org/journals/mcom/1990-55-192/S0025-5718-1990-1035925-1/S0025-5718-1990-1035925-1.pdf
    Lemma 2.4.

    INPUT : 
        -- A -- central simple algebra over Q of dimension N=n^2
    OUTPUT :
        -- a -- spliiting element of A, that is the minimal polynomial pi of a over Q 
                is separable (gcd(pi,pi')=1 i.e. pi has simple root in bar(Q)) 
                and deg(pi)=n
        -- pi -- the minimal polynomial of a
        -- g -- an irreducible factor of pi which has a real root

    REMARK : In particular a is a cycliq matrix view in Mn(Q)
            not sure if that is useful.
    """

    N = dimension(A)
    BA = list(A.basis())
    n = int(sqrt(N))
    H = 2*n*(n-1)
    for ite in range(nb_ite):
        t = [randint(1,H) for l in range(N)]
        a = sum( t[i]*BA[i] for i in range(N))
        pi = minimal_polynomial(a,A)
        if pi.degree()==n and gcd(pi,pi.derivative())==1:
            fact = [g for g,_ in pi.factor()]
            for g in fact:
                if len(g.real_roots()) >0:
                    return a,pi,g
    return "Not found"


def explicit_real_embedding(A):
    """
    INPUT : 
        -- A -- a Q-algebra given by structure constant assumed to be isomorphic to Mn(Q)
    OUTPUT :
        -- embedding_dict -- A dictionnary with key 0,.., N-1 where N = dim_A and value matrices in Mn(E)
                            for E some number field E = Q(gamma) representing an emmeding i : A -> Mn(E)
                            which map the i th element of the basis of A onto the matrix embedding_dict[i]
    """

    a,pi,g = real_splitting_element(A)
    E = NumberField(g, 'gamma')
    gamma = E.gen()

    R = PolynomialRing(E, 't')
    t = R.gen()

    B = FiniteDimensionalAlgebra(E,A.table()) # B = A⊗E over Q

    pi_E = R(pi)
    Q = pi_E // (t-gamma)

    a_in_B = sum(c*b for c,b in zip(get_coefficients(a,A),B.basis()))
    
    b = Q(a_in_B)

    return E,matrix_ring_iso_from_rank_one(B, b)


def push_lattice_embedding(A,Zbasis_I,E,embedding_dict):
    """
    INPUT :
        -- A -- algebra over Q isomorphic to Mn(Q)
        -- E -- real number field 
        -- embedding_dict -- representing a embedding A -> Mn(E)
    OUTPUT :
        -- Zbasis_J -- a list of n^2 vector in R^{n^2} representing the matrices
    """
    dim_A = dimension(A)
    Zbasis_J = []
    for i in range(dim_A):
         coords = get_coefficients(Zbasis_I[i],A)
         Zbasis_J.append(sum(coords[i]*embedding_dict[i] for i in range(dim_A)).list())
    return Zbasis_J


def approx_lattice(E, Zbasis_J, precision=1e-10):
    """
    Approximates a list of matrices with entries in a real number field by matrices with rational entries.

    INPUT:
        - E -- a real number field, defined as a quotient of Q[x] by an irreducible polynomial g.
        - Zbasis_J -- a list of matrices with entries in E.
        - precision -- a positive real number specifying the desired approximation accuracy.

    OUTPUT:
        - Zbasis_J_approx -- a list of matrices with rational entries that approximate the input matrices.
    """

    N = len(Zbasis_J[0])

    # Get the defining polynomial of the number field
    g = E.defining_polynomial()

    # Find the real roots of the polynomial g
    real_roots = g.real_roots()
    
    if not real_roots:
        raise ValueError("The defining polynomial of E has no real roots.")
        
    # Choose one of the real roots for the embedding
    gamma_real = real_roots[0]

    Zbasis_J_approx = []
    for vect in Zbasis_J:
        # Create a new matrix with the same dimensions, but with rational entries
        approx_vect = []
        for k in range(N):
            # Get the element from the number field
            element_in_E = vect[k]
                
            # Create a polynomial representation of the element
            poly_rep = element_in_E.polynomial()
                
            # Substitute the abstract root gamma with its real numerical value
            real_val = poly_rep(gamma_real)
                
            # Approximate the real number by a rational number
            approx_vect.append((QQ(real_val.n(digits=int(-log(precision, 10))))))
        
        Zbasis_J_approx.append(approx_vect)

    return Zbasis_J_approx