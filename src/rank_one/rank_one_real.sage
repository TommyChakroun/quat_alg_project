load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/isomorphism/explicit_iso_matrix_ring.sage")


#------------------------------------------------------------------------------------------
#
#            RANK-ONE ELEMENT IN A ⊗ R
#
#------------------------------------------------------------------------------------------



#--------------------------------------------------------------------
#  Step 1: Find a rank-one element in the algebra extended to a number field
#--------------------------------------------------------------------


def rank_one_MnR(A):
    """
    Finds a rank-one element in an algebra extended by a real number field.
    This is based on Lemma 2.4 from:
    https://www.ams.org/journals/mcom/1990-55-192/S0025-5718-1990-1035925-1/S0025-5718-1990-1035925-1.pdf

    INPUT: 
        - A: An algebra over Q, assumed to be isomorphic to Mn(Q).

    OUTPUT:
        - x: An element in B of rank one when viewed as a matrix in Mn(E).
        - B: The algebra A ⊗_Q E, where E is a number field E = Q(gamma).
         

    ALGORITHM:
        Let the basis of A be {a_1, ..., a_N}.
        By the Schwartz-Zippel lemma, a random linear combination 
        a = h_1*a_1 + ... + h_N*a_N (with small integer coefficients h_i)
        will, with high probability, have a minimal polynomial π_a that is 
        separable, of degree n, and has at least one real root.

        We find such an element 'a' and take an irreducible factor 'g' of its 
        minimal polynomial π_a that has a real root. Let gamma be a real root 
        of g, and define the number field E := Q(gamma).

        We then extend the algebra to B := A ⊗_Q E. In the polynomial ring 
        E[t], the minimal polynomial π_a now has (t-gamma) as a factor. We can 
        write π_a(t) = (t-gamma) * Q(t).
        
        The element x := Q(a) in B is then a rank-one element.
    """

    N = dimension(A)
    BA = list(A.basis())
    n = int(sqrt(N))
    H = 2*n*(n-1)

    nb_ite = 100
    for ite in range(nb_ite):
        h = [randint(1,H) for l in range(N)]
        a = A.sum( h[i]*BA[i] for i in range(N))

        pi = minimal_polynomial(a,A)

        if pi.degree()==n and gcd(pi,pi.derivative())==1:
            fact = [g for g,_ in pi.factor()]
            for g in fact:
                if len(g.real_roots()) > 0:
                    E = NumberField(g, 'gamma')
                    gamma = E.gen()
                    R = PolynomialRing(E, 't')
                    t = R.gen()
                    B = FiniteDimensionalAlgebra(E,A.table())
                    pi_E = R(pi) # View pi in E[t]
                    Q = pi_E // (t-gamma)
                    a_in_B = B.sum(c*b for c,b in zip(get_coefficients(a,A),B.basis()))
                    x = Q(a_in_B)
                    return x,B
    return "Not found"


#--------------------------------------------------------------------
#  Step 2: Deduce a rank-one idempotent
#--------------------------------------------------------------------


def rank_one_idempotent_MnR(A):
    """
    Constructs a rank-one idempotent from a rank-one element.

    INPUT:
        - A: An algebra over Q, assumed to be isomorphic to Mn(Q).

    OUTPUT:
        - e: A rank-one idempotent element in B.
        - B: The algebra A ⊗_Q E, where E is a number field E = Q(gamma).
    
    ALGORITHM:
        First, a rank-one element 'x' is found using rank_one_MnR(A).
        The equation x*a*x = x is then solved for 'a'. Such an 'a' exists
        because the space x*B*x is one-dimensional.
        The element e := a*x is then a rank-one idempotent, since
        e^2 = (a*x)*(a*x) = a*(x*a*x) = a*x = e.
    """
    x,B = rank_one_MnR(A)
    a = solve_xax_eq_x(x,B)
    e = a*x
    return e,B


#--------------------------------------------------------------------
#  Step 3: Recover a rational zero divisor by approximation
#--------------------------------------------------------------------


def approx_number_field_element(E, x, precision=1e-10):
    """
    Approximates an element of a number field with a rational number.

    INPUT:
        - E: A real number field, defined as a quotient of Q[t] by an 
             irreducible polynomial g.
        - x: An element in E.
        - precision: The desired precision for the approximation.

    OUTPUT:
        - r: A rational number (element of QQ) that approximates the real
             value of x under a fixed embedding of E into R.
    
    ALGORITHM:
        The defining polynomial g of E is retrieved. An embedding of E into R
        is chosen by picking the first real root of g. The element x, represented
        as a polynomial in the generator of E, is evaluated at this real root.
        The resulting real number is then converted to a rational number with
        the specified precision.
    """
    g = E.defining_polynomial()
    real_roots = g.real_roots()
    if not real_roots:
        raise ValueError("The defining polynomial of E has no real roots.")
    
    # This chooses a specific real embedding of E.
    gamma_real = real_roots[0] 
    
    poly_rep = x.polynomial()
    real_val = poly_rep(gamma_real)
    
    # Convert the numerical value to a Sage rational number (QQ).
    num_digits = int(-log(precision, 10))
    r = QQ(real_val.n(digits=num_digits))
    return r


def heuristic_zero_divisor_from_real(A):
    """
    Heuristically finds a zero divisor in A by approximating an idempotent
    from a real extension field.

    INPUT:
        - A: An algebra over Q, assumed to be isomorphic to Mn(Q).

    OUTPUT:
        - a: An element in A that is likely a zero divisor.
    
    ALGORITHM:
        Let A = Q*a_1 ⊕ ... ⊕ Q*a_N.
        1. We construct a rank-one idempotent e = c_1*a_1 + ... + c_N*a_N, 
           where the coefficients c_i belong to a real number field E = Q(gamma).
           This is done by calling rank_one_idempotent_MnR.
        
        2. We compute high-precision rational approximations for each 
           coefficient c_i, resulting in an element y in A whose coefficients 
           are rational numbers with large numerators and denominators.
           y is close to e, but is unlikely to be a zero divisor.

        3. To find a "simpler" element, we need to find good rational 
           approximations for the coefficients of y, ideally with a common
           denominator. This step requires a lattice-based algorithm like LLL
           to find small integer relations.

        4. The resulting element 'a' in A, constructed with these new "simple"
           rational coefficients, is our candidate for a zero divisor. While we
           don't expect a^2 = a, we hope that 'a' is a non-trivial zero divisor.
    """

    e, B = rank_one_idempotent_MnR(A)
    E = B.base_ring()

    coords = get_coefficients(e, B)

    # Step 2: Get high-precision rational approximations of the algebraic coefficients.
    high_precision_approx = [approx_number_field_element(E, c) for c in coords]

    better_approx = [r.n().nearby_rational(max_error = 1e-10) for r in high_precision_approx]  # This should be the list [s_1/d, ..., s_N/d]

    # Step 4: Construct the candidate zero divisor from the LLL-found approximation.
    a = A.sum(r*b for r,b in zip(better_approx, A.basis()))

    return a

