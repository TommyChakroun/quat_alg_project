load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("core/rank_one.sage")
load("maximal_orders/maximal_orders_utilities.sage")
load("maximal_orders/find_maximal_orders.sage")
load("explicit_embeding_real_number/explicit_embedding.sage")
load("core/identify_quaternion_algebra.sage")


#------------------------------------------------------------------------------------------
#
#            ZERO DIVISOR IN ALGEBRA ISOMORPHIC TO M4(Q)
#  https://www.sciencedirect.com/science/article/pii/S0747717107000089?ref=pdf_download&fr=RR-2&rr=94d3039f3a517bf3
#
#------------------------------------------------------------------------------------------


## Utilities

def is_central(x,A):
    BA = A.basis()
    for a in BA:
        if a*x != x*a:
            return False,a
    return True,BA[0]


def centralizer(a,A):
    """
    INPUT :
        -- A -- an algebra isomorphic to M4(Q)
        -- a -- an element of A
    OUTPUT:
        -- basis_Za -- a basis of Za = {x in A | xa=ax}
        -- Za -- a finite dimensional algebra given byh structure constant representing Za with the basis basis_Za
        -- lift -- function from Za to A the inlcusion.
    """
    dim_A = dimension(A)
    BA = list(A.basis())
    a = get_coefficients(a,A)

    T = structure_constants(A,BA)

    rows = []
    for k in range(dim_A):
        rows.append([sum(a[i]*T[j][i][k] - a[i]*T[i][j][k] for i in range(dim_A)) for j in range(dim_A)])
    
    M = Matrix(QQ,rows)

    X = M.right_kernel().basis()

    basis_Za = [sum(x[i]*BA[i] for i in range(dim_A)) for x in X]

    Za,lift = subalgebra_from_subspace(A,basis_Za)

    def push(x):
        BZa = Za.basis() 
        coords = coordinate(x,A,basis_Za)
        return sum(c*b for c,b in zip(coords,BZa))


    return Za,push,lift


def find_structured_basis(A, a):
    """return a basis of the form [e1,a*e1,..,e4,a*e4] assuming a has an irreducible minimal polynomial of degree 2"""
    V = VectorSpace(A.base_ring(), A.dimension())
    building_basis = [A.one(), a]
    for _ in range(3):
        W = V.subspace([b.vector() for b in building_basis])
        for e in A.basis():
            if e.vector() not in W:
                building_basis.extend([e, a * e])
                break
    return building_basis

def change_field(A,a):
    """
    INPUT:
        -- A -- an algebra of dimension 8 over Q
        -- a -- in A with minimal polynomial irreducible of degree 2, so Q(a) is a field
    OUTPUT :
        -- B -- the algerba A view with the field Q(a)
        -- psuh -- A to B
        -- lift -- B to A
    """
    pi = minimal_polynomial(a,A)
    E = NumberField(pi,'aa')
   
    aa = E.gen()

    structured_basis = find_structured_basis(A,a)

    C = []
    for j in range(4):
        rows = []
        for i in range(4):
            fi = structured_basis[2*i]
            fj = structured_basis[2*j]
            coords = coordinate(fi*fj,A,structured_basis)
            rows.append([coords[2*k]+coords[2*k+1]*aa for k in range(4)])
        M = Matrix(E,rows)
        C.append(M)
    B = FiniteDimensionalAlgebra(E,C)
    BB = B.basis()

    def push(a):
        coords = coordinate(a,A,structured_basis)
        return sum((coords[2*k]+coords[2*k+1]*aa)*BB[k] for k in range(4))

    def lift(b):
        coords = get_coefficients(b,B)
        return sum(coords[k].vector()[0]*structured_basis[2*k]+coords[k].vector()[1]*structured_basis[2*k+1] for k in range(4))

    
    return B,push,lift
            





## Quadratic element : a in A such that pi_a has degree 2.

def quadratic_element(A):
    """
    INPUT :
        -- A -- an algebra isomrophic to M4(Q)
    OUTPUT :
        -- a -- element of A such with a minimal polynomial of degree 2
                or a is a zero divisor in A
    """
    BA = list(A.basis())
    c = BA[0]
    if is_central(BA[0],A):
        c = BA[1]

    pi = minimal_polynomial(c,A)
    c3 = pi.coefficient(3)

    c = c + c3/4 * A.one()

    pi = minimal_polynomial(c,A)

    if pi.degree()==2:
        return c
    
    if pi.coefficient(1)==0:
        print("pi_c = T^4+T ^2 +1")
        return c**2

    _,y = is_central(c,A)

    x = c*y-y*c

    if right_rank(A,x)<16:
        return x

    c_prime = x*c*(x**(-1))

    a = -(c+c_prime)

    return a**2



## Zero divisor in A isom M4(Q)

def sub_quat_alg_M4Q(A):
    a = quadratic_element(A)

    ## If a directly give a zero divisor this is ok

    pi = minimal_polynomial(a,A)
    if not pi.is_irreducible():
        g = pi.factor()[0][0]
        return g(a)

    ## Now a is a quadratic element with minimal polynomial pi irreducible so Q(a) is a field

    Ca,push,lift = centralizer(a,A)  # should have dimension 8 over Q

    A2,push_2,lift_2 = change_field(Ca,push(a))

    E = A2.base_ring()
    i,j,a,b = quaternion_structure(A2)

    A3 = QuaternionAlgebra(E,a,b)

    def lift_final(x):
        """
        lift a x in (a,b | E) in our original algebra M4Q
        """
        t,s,v,w = get_coefficients(x,A3) 
        x_in_A2 = t*A2.one()+s*i+v*j+w*i*j
        x_in_Ca = lift_2(x_in_A2)
        x_in_A = lift(x_in_Ca)
        return x_in_A

    return A3,lift_final