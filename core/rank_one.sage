load("utilities/utilities.sage")
load("utilities/algebra_type.sage")

#------------------------------------------------------------------------------------------
#
#            RANK ONE ELEMENT IN ALGEBRA ISOMORPHIC TO Mn(Q)
#
#------------------------------------------------------------------------------------------

##---List of functions---##

# right_rank(A,x)

def right_rank(A, x):
    """
    INPUT :
        -- A -- an algebra
        -- x -- an element of A
    OUTPUT :
        -- r -- the rank of A ->A,a ->ax i.e the dimension of Ax over F
    """
    F = A.base_ring()
    basis_A = A.basis()

    rows = []

    for e in basis_A:
        y = e * x
        coords = coordinate(y, A, basis_A)
        rows.append(coords)
    
    M = Matrix(F, rows)
    r = M.rank()
    return r

def min_rank(A):
    r_min = dimension(A)
    x_min = A.one()
    max_ite = 1000
    for i in range(max_ite):
        x = A.random_element()
        r = right_rank(A,x)
        if r !=0 and r<r_min:
            r_min=r
            x_min = x 
    return r_min,x_min

def rank_one(A):
    N = dimension(A)
    max_ite = 1000
    for i in range(max_ite):
        x = A.random_element()
        r = right_rank(A,x)
        if r**2 == N:
            return x
    return "not found"

    
        
