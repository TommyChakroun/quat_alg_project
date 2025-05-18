load("/quat_alg_project/utilities/utilities.sage")
load("/quat_alg_project/utilities/algebra_type.sage")

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