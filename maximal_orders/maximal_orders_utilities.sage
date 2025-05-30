load("utilities/utilities.sage")
load("utilities/algebra_type.sage")



#------------------------------------------------------------------------------------------
#
#            MAXIMAL ORDERS UTILITIES
#
#------------------------------------------------------------------------------------------


## Modular matrix kernel

def kernel_mod(M,d):
    """
    INPUT :
        -- M -- a rectangular matrix of size m x n with coefficients in Z
        -- d -- a positive integer
    OUTPUT :
        -- Zbasis_sol -- a list of vector of Z^n which is a Z basis of the set {X in Z^n | MX = 0 mod p}
    """
    m = M.nrows()
    n = M.ncols()
    N = min(n,m)

    D,U,V = M.smith_form(transformation=True)

    Y = [Matrix(ZZ, n, 1, [d/gcd(d, D[i][i]) if j == i else 0 for j in range(n)]) for i in range(N)]
    for i in range(N,n):
        Y.append(Matrix(ZZ, n, 1, [1 if j == i else 0 for j in range(n)]))
    X = [V* Y[i] for i in range(n)]

    return X


## Latices and left order 

def lattice_LLL(B,Zbasis_I):
    """
    INPUT :
        -- B -- a finite dimensional algebra over Q
        -- Zbasis_I -- a list [e1,..,eN] which is a Q-basis of B and representing the lattice I := Ze1⊕... ⊕ZeN 
    OUTPUT : 
        -- Zbasis_J -- a list [f1,..,fN] which is a Q-basis of B and representing the lattice J := Zf1⊕... ⊕ZfN such that J = I and Zbasis_J is "simpler"
    """
    basis_B = list(B.basis())
    N = dimension(B)

    rows = []
    for e in Zbasis_I:
        rows.append(get_coefficients(e,B))
    A = Matrix(QQ,rows)
    M = A.LLL()
    Zbasis_J = [ sum(M[i][k]*basis_B[k] for k in range(N)) for i in range(N)]
    return Zbasis_J

def is_in_lattice(B,Zbasis_I,x):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_I -- a list [e1,..,eN] which is a Q-basis of B and representing the lattice I := Ze1⊕... ⊕ZeN 
        -- x -- an element of B
    OUTPUT :
        -- boolean -- True if x is in B, False otherwise.
    """
    coords = coordinate(x,B,Zbasis_I)
    for r in coords:
        if r not in ZZ:
            return False
    return True


def is_sub_lattice(B,Zbasis_I,Zbasis_J):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_I -- a list [e1,..,eN] which is a Q-basis of B and representing the lattice I := Ze1⊕... ⊕ZeN 
        -- Zbasis_I -- a list [f1,..,fN] which is a Q-basis of B and representing the lattice I := Zf1⊕... ⊕ZfN 
    OUTPUT :
        -- boolean -- True if I is a sublattice of J, False otehrwise
    """
    for e in Zbasis_I:
        if not  is_in_lattice(B,Zbasis_J,e):
            return False
    return True


def is_order(B,Zbasis_I):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_I -- a list [e1,..,eN] which is a Q-basis of B and representing the lattice I := Ze1⊕... ⊕ZeN
    OUTPUT :
        -- boolean -- True if I is an order of B
    """
    if not is_in_lattice(B,Zbasis_I,B.one()):
        return False

    for e in Zbasis_I:
        for f in Zbasis_I:
            x = f*e
            if not is_in_lattice(B,Zbasis_I,x):
                return False
    return True


def are_equals_lattices(B,Zbasis_I,Zbasis_J):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_I -- a list [e1,..,eN] which is a Q-basis of B and representing the lattice I := Ze1⊕... ⊕ZeN 
        -- Zbasis_I -- a list [f1,..,fN] which is a Q-basis of B and representing the lattice I := Zf1⊕... ⊕ZfN 
    OUTPUT :
        -- boolean -- True if I is a sublattice of J, False otehrwise
    """
    return is_sub_lattice(B,Zbasis_I,Zbasis_J) and is_sub_lattice(B,Zbasis_J,Zbasis_I)


def left_order(B,Zbasis_I,reduced = True):
    """
    INPUT :
        -- B -- a central simple algebra over Q given by structure constant of dimension N
        -- Zbasis_I -- a list [e1,..,eN]  which is a Q-basis of B representing the 
                       Z-lattice I = Ze1⊕... ⊕ZeN
    OUTPUT :
        -- Zbasis_O -- a list [f1,..,fN]  which is a Q-basis of B representing the 
                       Z-order O = Zf1⊕... ⊕ZfN such that O=O_L(I)
    """
    N = dimension(B)
    one_B = B.one()
    coords_one = coordinate(one_B, B, Zbasis_I)
    s = lcm([r.denominator() for r in coords_one])


    C = structure_constants(B,Zbasis_I)

    rows = []
    for j in range(N):
        for k in range(N):
            rowjk = []
            for i in range(N):
                rowjk.append(C[j][i][k]/s)
            rows.append(rowjk)
    M = Matrix(QQ,rows)
    d = lcm([a.denominator() for a in M.list()])
    A = d*M

    A = Matrix(ZZ,A)
    X = kernel_mod(A,d) # Z-basis of solution of MX = 0 mod d

    Zbasis_OLI = [(1/s) * sum(X[i][k,0] * Zbasis_I[k] for k in range(N)) for i in range(N)]

    if not reduced:
        return Zbasis_OLI
    return lattice_LLL(B,Zbasis_OLI)
        


## Reduced trace and discriminant

def reduced_trace(B,x):
    """
    INPUT :
        -- B -- a central simple algebra given by structure constant
        -- x -- an element of B
    OUTPUT :
        -- t -- the reduced trace of x = 1/n*trace (b ->xb) where n i such that dim B = n^2
    """
    basis_B = list(B.basis())
    N = dimension(B)
    n = int(sqrt(N))

    trd = 0
    for i in range(N):
        coords = get_coefficients(x*basis_B[i],B)
        trd = trd + coords[i]
    
    trd = trd/n
    return trd

def discriminant(B,Zbasis_I):
    """
    INPUT :
        -- B -- a finite dimensional central simple algebra over Q
        -- Zbasis_I -- a list [e1,..,eN] which is a Q-basis of B and representing the lattice I := Ze1⊕... ⊕ZeN 
    OUTPUT :
        -- d -- in Z the discriminant of I (modulo +-1).
    """
    N = dimension(B)
    e = Zbasis_I
    rows = []
    for i in range(N):
        row = []
        for j in range(N):
            row.append(reduced_trace(B,e[i]*e[j]))
        rows.append(row)
    M = Matrix(QQ,rows)
    return M.determinant()
            



## Quotient of algebra and kernel of Z module map



def finite_algebra_from_order(B,Zbasis_O,p):
    """
    INPUT :
        -- b -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list [e1,..,eN] which is a Q-basis of B and such that O := Ze1⊕... ⊕ZeN is a Z-order of B
    OUTPUT :
        -- A -- an instance of FiniteDimensional algebra representing the finite Fp-algebra A = Fp e1⊕... ⊕Fp eN with the multiplication from O
        -- pi -- a lambda function [0,..,N-1] -> A which represent pi : Z^n -> A the abelian group homomrohism O -> A : pi(i) = fi = ei in A
    """
    table = []
    for e in Zbasis_O:
        rows = []
        for f in Zbasis_O:
            coords = coordinate(e*f,B,Zbasis_O) # a priori coordinate in Q but in Z since O is an order
            coords_mod_p = [GF(p)(a) for a in coords]
            rows.append(coords_mod_p)
        M = Matrix(GF(p), rows)
        table.append(M)

    cat = AlgebrasWithBasis(GF(p)).Associative().FiniteDimensional()
    A = FiniteDimensionalAlgebra(GF(p),table,category=cat)
    f = A.basis() # fi = ei view in A
    pi = lambda i : f[i]
    return A,pi


def quotient_algebra_ideal(A,basis_I):
    """
    INPUT : 
        -- A -- a finite dimensional algebra given by structure constants (over a finite field in our case)
        -- I -- a list of element of A representing a basis of an ideal I
    OUTPUT :
        -- C -- a finite dimensional algebra given by structure constant represnting A/I
        -- pi -- function A ->C  representing the projection A ->A/I
    """

    I = A.ideal(basis_I)
    pi = A.quotient_map(I)
    C = pi.codomain()
    return C,pi


def kernel_Z_mod_map(C, phi,N):
    """
    INPUT:
        -- C -- a finite-dimensional algebra over a finite field, given by structure constants.
        -- phi -- a function from [0, ..., N-1] to C representing a Z-module homomorphism φ: Z^N → C.

    OUTPUT:
        -- basis_Ker --: a list of vectors in Z^N forming a Z-basis of the kernel of φ.
    """

    # Determine the dimension M of the codomain C
    M = C.dimension()
    p = C.base_ring().order()

    # Initialize an empty list to store the columns of the matrix representation of φ
    columns = []

    # For each basis element e_i of ℤ^N, compute φ(e_i) and express it as a vector in C
    for i in range(N):
        phi_ei = phi(i)
        coords = get_coefficients(phi_ei, C) # cooridnate in Fp 
        coords_int = [int(c) for c in coords]
        columns.append(coords_int)

    # Construct the matrix representation of φ over ℤ
    M_matrix = Matrix(ZZ, columns).transpose()

    # Extract the basis vectors of the kernel
    basis_Ker =  kernel_mod(M_matrix,p)

    return basis_Ker


