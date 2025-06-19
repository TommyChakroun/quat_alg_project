load("utilities/utilities.sage")
load("utilities/algebra_type.sage")

#------------------------------------------------------------------------------------------
#
#            NORMALIZE QUADRATIC FORM ON Z_(2)
#
#------------------------------------------------------------------------------------------



def normalize_quadratic_form_over_Z_loc_2(Q, e):
    """
    Normalize a quadratic form Q over Z_(2).

    INPUT:
        - Q: quadratic form Q: Q^n -> Q assumed to have coefficients in Z_(2), viewed as Q: Z_(2)^n -> Z_(2) because the real implementation of Z_(2) in Sage is broken.
        - e: list of vectors [e0, ..., en-1] forming a basis of Z_(2)^n.

    OUTPUT:
        - list of vectors [f0, ..., fn-1] forming a normalized basis.
    """


    T = lambda x,y : Q(x+y)-Q(x)-Q(y)
    ord = lambda x: QQ.valuation(2)(x)

    n = len(e)

    if n==1:
        return e
    
    # Step 1 of John Voight's paper
    is_Q_zero = True
    min_val = +Infinity
    min_ij = None
    for i in range(n):
        for j in range(i, n):
            val = ord(T(e[i], e[j]))
            if T(e[i], e[j]) != 0 and val < min_val:
                is_Q_zero = False
                min_val = val
                min_ij = (i, j)

    if is_Q_zero:
        return e

    i, j = min_ij

    # Step 2 of John Voight's paper
    if i == j:
        f1 = e[i]
        # Step 3 of John Voight's paper
        rest = [e[k]- (T(f1, e[k]) / T(f1, f1)) * f1 for k in range(1,n)]
        # Step 5 of John Voight's paper
        return [f1] + normalize_quadratic_form_over_Z_loc_2(Q, rest)


    # Step 4 of John Voight's paper
    ord_val = ord(T(e[i], e[j]))
    f1 = 2**ord_val / T(e[i], e[j])*e[i]
    f2 = e[j]
    e[i] = e[0]
    e[j]=e[1]
    
    d = T(f1, f1) * T(f2, f2) - T(f1, f2)**2

    rest = []
    for k in range(2,n):
        ek = e[k]
        tf1f2 = T(f1, f2)
        tf1f1 = T(f1, f1)
        tf2f2 = T(f2, f2)
        tf1ek = T(f1, ek)
        tf2ek = T(f2, ek)
        tk = tf1f2 * tf2ek - tf2f2*tf1ek
        uk = tf1f2 * tf1ek  - tf1f1 * tf2ek
        fk = ek + (tk / d) * f1 + (uk / d) * f2
        rest.append(fk)

    # Step 5 of John Voight's paper
    return [f1, f2] + normalize_quadratic_form_over_Z_loc_2(Q, rest)



def normalize_symetric_matrix_over_Z_loc_2(S):
    """
    """
    R = S.base_ring()
    n = S.nrows()
    Q = QuadraticForm(R,2*S) # Q(x) = x^t S x

    e = identity_matrix(R,n)
    basis = [vector(R,e[i]) for i in range(n)]
    new_basis = normalize_quadratic_form_over_Z_loc_2(Q, basis)

    P = Matrix(R,new_basis).transpose()
    D = P.transpose()*S*P

    return D,P



