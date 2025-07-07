load("utilities/utilities.sage")
load("utilities/algebra_type.sage")



#------------------------------------------------------------------------------------------
#
#            EXPLICIT ISOMROPHISM PROBLEM FOR MATRIX RING
#               FIND phi : C -> M_N(F)
#------------------------------------------------------------------------------------------

##---List of functions---##

# matrix_ring_iso_from_algebra_iso(A,B,isom_dict,basis_A,basis_B)
# basis_of_translated(A, x)
# matrix_ring_iso_from_rank_one(A, x)



##--- C =  A ⊗ B^op and we know f : A ->B isomorphism ---##

def matrix_ring_iso_from_algebra_iso(A,B,isom_dict):
    """
    INPUT :
        -- A,B -- central simple algebra of dim n^2 over F
        -- isom_dict -- dictionarry with key 0,1,.., n²-1 and value in B representing an isomorphism f : A -> B , f(A.basis()[i]) =  isom_dict[i]
    OUTPUT : 
        -- isom_to_mat_dict -- dictionarry with key the couples (i,j) 0<=i<n²-1 ,0<=j<= n²-1 and value in Mn²(F) representing an isomorphism phi : A ⊗ B^op -> M_n^2(F) : 
            phi(ei⊗fj) = isom_to_mat_dict[(i,j)] where e0,..en²-1 = A.basis() and f0,..,fn²-1 = B.basis().
    """
    F = A.base_ring()
    N = dimension(A) # = n²
    basis_B = B.basis()
    
    isom_to_mat_dict = {}
    for i in range(N):
        for j in range(N):
            g = lambda z : isom_dict[i]*z*basis_B[j]
            M = matrix_in_basis(B, B, basis_B, basis_B, g) # matrix of g :B -> B, z ->f(ei)*z*fj in  basis_B at the source and the target
            isom_to_mat_dict[(i,j)] = M

    return isom_to_mat_dict


##--- C whatever and we know x in C of rank 1 ---#

def basis_of_translated(A, x):
    """
    Return a basis of Ax.
    INPUT :
        -- A -- an algebra
        -- x -- an element of A
    OUTPUT :
        -- basis_of_Ax -- a list of element of A which is a basis of Ax.
    """
    F = A.base_ring()
    basis_A = A.basis()
    n = dimension(A)

    ## construct the matrix of e0*x,..,en-1*x in a basis_A
    rows = []

    for e in basis_A:
        y = e * x
        coords = coordinate(y, A, basis_A)
        rows.append(coords)
    
    M = Matrix(F, rows)
    
    ## extract a basis of Ax
    basis_of_Ax = []
    independent_vectors = M.row_space().basis_matrix()
    for v in independent_vectors:
        basis_of_Ax.append(sum(v[i] * basis_A[i] for i in range(n)))

    return basis_of_Ax


def matrix_ring_iso_from_rank_one(A, x):
    """
    INPUT : 
        -- A -- a finite dimensional algebra assume central simple
        -- x -- an element of A such that dim(Ax)=n where dim(A)=n^2
    OUTPUT :
        -- isom_dict --  a dictionary with keys 0,1,...,n²-1 and value in Mn(F) which represent an isomorphism f : A -> Mn(F) given by f(A.basis()[i]) = isom_dict[i]
    """
    F = A.base_ring()
    N = dimension(A)
    n = right_rank(x, A)
    basis_Ax = basis_of_translated(A, x)
    basis_A = A.basis()

    if n**2 != N:
        return "Ax has not dimension sqrt(dim(A))", None
    
    isom_dict = {}

    for i, e in enumerate(basis_A):
        rows = []
        for z in basis_Ax:
            y = e * z  # y is still in Ax
            coords = coordinate(y, A, basis_Ax)
            rows.append(coords)
        M = Matrix(F, rows).transpose()
        isom_dict[i] = M

    return isom_dict




