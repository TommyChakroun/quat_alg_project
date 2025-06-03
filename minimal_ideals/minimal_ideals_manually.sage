load("maximal_orders/maximal_orders_utilities.sage")
load("utilities/algebra_type.sage")
load("minimal_ideals/idempotent_to_ideals.sage")


#------------------------------------------------------------------------------------------
#
#            MINIMAL IDEALS OF A SEMI SIMPLE ALGEBRA
#
#------------------------------------------------------------------------------------------



## Utilities 

def minimal_polynomial(a,A):
    """
        -- A -- a finite dimensional algebra over a field F
        -- a -- an element of A
        -- basis_field -- a list of element of 
    """

    # Get basis of the algebra A
    basis = A.basis()
    n = len(basis)
    F = A.base_ring()

    # Multiplication by a as a matrix
    M = Matrix(F, n)
    for j, bj in enumerate(basis):
        prod = a * bj
        coeffs = get_coefficients(prod,A)
        M.set_column(j, vector(F, coeffs))

    # Minimal polynomial of the multiplication matrix
    return M.minimal_polynomial()


def split_from_idemp(A,list_idemp):
    """
    INPUT :
        -- A -- a commutative algebra
        -- list_idemp -- a list e1,..,es of orthogonal idempotent of A
    OUTPUT :
        -- Subalgebras -- a list of couple (Ai,basis_Ai) where Ai is an algebra given by structure constant 
                        representing the algebra Aei and basis_Ai is a list of element of A representing the
                        basis of Ai view in A 
    """
    basis_A = A.basis()
    F = A.base_ring()
    dim_A = dimension(A)
    Subalgebras = []
    for e in list_idemp:
        SetGeneratorsIdeal = [a*e for a in basis_A]

        ## Write vector of SetGeneratorsIdeal in each row of a matrix
        rows = []
        for x in SetGeneratorsIdeal:
            rows.append(get_coefficients(x,A))
        M = Matrix(F,rows)
        ## Extract free family of SetGeneratorsIdeal

        basis_Ai = []
        independent_vectors = M.row_space().basis_matrix()
        for v in independent_vectors:
            basis_Ai.append(sum(v[i] * basis_A[i] for i in range(dim_A)))
        
        table = structure_constants_subspace(A, basis_Ai)
        Ai = FiniteDimensionalAlgebra(F,table)
       
        Subalgebras.append((Ai,basis_Ai))
    return Subalgebras






## Idempotents


def idempotents_commutative_perso_split(Z,eps = 0.000000000001):
    """
    INPUT :
        -- Z -- a commutative semi simple algebra over a fiel F given by structure constant
                equivalently Z is isomorphic to a product K1 x .. x Kr of field extension Ki/F
        -- eps -- a precision
    OUTPUT :
        -- [e1,..,es] -- a list of non zero orthogonal idempotent of r
                        the list is probably well of size s=r the abstract number of field extension
                        i.e. the list is probably made of primite orhtogonal indempotent.
                        But not sure 100%   
    ALGORITHM:
        let n = dim Z and p the cardinal of the base field of Z.
        Choose nb_its an integer such that (1/p^{n-1})^nb_ite < eps;
    """
    n = dimension(Z)
    F = Z.base_ring()
    p = F.order()

    if n==1:
        return [Z.one()]

    nb_ite = int(- log(eps)/((n-1)*log(p)) +1)
    N = 0

    ## We try to find strict idempotents in Z


    for N in range(nb_ite):
        alpha = Z.random_element()
        pi =  minimal_polynomial(alpha,Z)
        fact = pi.factor()
        if len(fact)>2:
            h = []
            moduli = [P for P, _ in list(fact)] 

            s = len(fact)
            for i in range(s):
                residues = [1 if j == i else 0 for j in range(s)]
                h.append(crt(residues, moduli))
            list_idemp = [h[i](alpha) for i in range(s)]
            subalgebras = split_from_idemp(Z,list_idemp)
            ListIdempFinal = []

            for A,BA in subalgebras:
                idemp_in_A = idempotents_commutative_perso_split(A,eps)
                dim_A = dimension(A)
                for e in idemp_in_A:
                    coords = get_coefficients(e,A)
                    ListIdempFinal.append(sum(coords[i]*BA[i] for i in range(dim_A)))
            return ListIdempFinal

    ## If we didn't find, Z is probably just a field extension of F so we return the unit of Z
    
    return [Z.one()]


def idempotents_commutative_perso(Z,nb_ite = 100):
    """
    INPUT :
        -- Z -- a commutative semi simple algebra over a fiel F given by structure constant
                equivalently Z is isomorphic to a product K1 x .. x Kr of field extension Ki/F
    OUTPUT :
        -- [e1,..,es] -- a list of non zero orthogonal idempotent of r
                        the list is probably well of size s=r the abstract number of field extension
                        i.e. the list is probably made of primite orhtogonal indempotent.
                        But not sure 100%   
    """
    nb_max_factor = 0
    list_max_factor = []
    good_alpha = Z.zero()
    for N in range(nb_ite):
        alpha = Z.random_element()
        pi =  minimal_polynomial(alpha,Z)
        fact = pi.factor()
        if len(fact)>nb_max_factor:
            good_alpha = alpha
            list_max_factor = fact
            nb_max_factor = len(fact)
    
    h = []
    moduli = [P for P, _ in list_max_factor] 
    for i in range(nb_max_factor):
        residues = [1 if j == i else 0 for j in range(nb_max_factor)]
        h.append(crt(residues, moduli))

    return [h[i](good_alpha) for i in range(nb_max_factor)]


def central_idempotents_perso(A,nb_ite=100):
    """
    INPUT :
        -- A -- a semi simple algebra over F
    OUTPUT :
        -- [e1,..,es] -- a list of non zero orthogonal idempotent of r
                        the list is probably well of size s=r the abstract number of field extension
                        i.e. the list is probably made of primite orhtogonal indempotent.
                        But not sure 100%   

    """
    Z,lift = center_perso(A)
    return [lift(e) for e in idempotents_commutative_perso(Z,nb_ite)]








## Minimal Ideal Perso


def minimal_ideals_perso(A,nb_ite =100):
    return idempotents_to_ideals(A,central_idempotents_perso(A,nb_ite))









## Old functions


def minimal_ideals_commutative(A):
    """
    INPUT :
        -- A -- a finite dimensional commutative algebra over a finite field given by structure constant
    OUTPUT :
        -- MinIdealsList -- a list of list of element of A which are the basis of the minimal two sided ideal of A
    """
    if split_commutative_algebra(A) == "IsField":
        return [A.basis()]
    
    basis_I,basis_J = split_commutative_algebra(A)

    ## We compute recursively the minimal ideal of I and J.
    
    I = finite_dimensional_algebra_format(A, basis_I)
    J = finite_dimensional_algebra_format(A, basis_J)

    dim_I = I.dimension()
    dim_J = J.dimension()

    list_basis_ideal_I = minimal_ideals_commutative(I)
    list_basis_ideal_J = minimal_ideals_commutative(J)

    list_basis_ideal_A = []

    for b in list_basis_ideal_I:
        lift_b = []
        for e in b:
            cooefs = get_coefficients(e,I)
            lift_b.append(sum(coefs[i]*basis_I[i] for i in range(dim_I)))

    for b in list_basis_ideal_J:
        lift_b = []
        for e in b:
            cooefs = get_coefficients(e,J)
            lift_b.append(sum(coefs[i]*basis_J[i] for i in range(dim_J)))

    return list_basis_ideal_A


def minimal_ideals_manually(A):
    """
    INPUT : 
        -- A -- a finite dimensional algebra over a finite field given by structure constant
    OUTPUT :
        -- MinIdealsList -- a list of list of element of A which are the basis of the minimal two sided ideal of A
    """
    F = A.base_ring()

    basis_Z,Z = center(A)
    basis_A = A.basis()

    dim_Z = dimension(Z)
    dim_A = dimension(A)
    

    ListBasisMinIdealCenter = minimal_ideals_commutative(Z)
    MinIdealsList = []

    for basis_ZAi_in_Z in ListBasisMinIdealCenter:   # visit all basis of Z(Ai) where Z = Z(A1) ⊕ ...  ⊕ Z(Ar)
        ## View the basis of Z(Ai) in A 
        basis_ZAi=[]
        for e in basis_ZAi_in_Z :
            coords = get_coefficient(e,Z)
            basis_ZAi.append(sum(coords[i]*basis_Z[i] for i in range(dim_Z)))

        ## Construct the set of generators of ideal Ai = Z(Ai)A :  {bj*ak for bj a basis of Z(Ai) and ak a basis of A }
        dim_ZAi = len(basis_ZAi)
        SetGeneratorsIdeal = []

        for j in range(dim_ZAi):
            for k in range(dim_A):
                SetGeneratorsIdeal.append(basis_ZAi[j]*basis_A[k])

        ## Write vector of SetGeneratorsIdeal in each row of a matrix
        rows = []
        for x in SetGeneratorsIdeal:
            rows.append(get_coefficient(x,A))
        M = Matrix(F,rows)
        ## Extract free family of SetGeneratorsIdeal

        basis_Ai = []
        independent_vectors = M.row_space().basis_matrix()
        for v in independent_vectors:
            basis_Ai.append(sum(v[i] * basis_A[i] for i in range(dim_A)))
        
        ## Add the basis of ideal to the list

        MinIdealsList.append(basis_Ai)

    return MinIdealsList

