load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/iso_splitting_algebra/rank_one/rank_one_MnQ.sage")
load("src/iso_splitting_algebra/explicit_iso_matrix_ring.sage")

#------------------------------------------------------------------------------------------
#
#           EXPLICIT ISOMROPHISM PROBLEM FOR QUATERNION ALGEBRA
#               FIND f : A -> B isomorphism
#------------------------------------------------------------------------------------------



##--- When know phi : A ⊗ B^op -> M_n²(F) ---#

def quat_alg_iso_from_matrix_ring_iso(A,B,isom_dict,cpt=False):
    """
    INPUT : 
        -- A,B -- central simple algebra
        -- isom_dict -- a dictionnary with keys (i,j) 1<=i<=n²-1, 1<=j<=n^-1 and values in M_n²(F) 
                        representing an isomorphism phi: A ⊗ B^op to M_n²(F) given by isom_dict[(i,j)] = phi(ei⊗fj) 
                        where (ei)=A.basis() (fj) = B.basis() 
    OUTPUT :
        -- new_isom_dict -- a dictionnary with keys 0,1,..,n²-1 and values in B
                            representing an isomorphism f: A ->B given by new_isom_dict[i] = f(ei) 
                            where (ei)=A.basis() 
    """

    F = A.base_ring()
    N = dimension(A) # = n²
    V = VectorSpace(F,N)

    basis_V = V.basis()
    basis_A = A.basis()
    basis_B = B.basis() 

    new_isom_dict = {}

    n=0
    go = True
    while go and n<1000:
        v = V.random_element()
        phi = lambda a : isom_dict[(basis_A.index(a),0)]*v  # well defined for a in basis_A
        psi = lambda b : isom_dict[(0,basis_B.index(b))]*v  # well defined for b in basis_B

        M = matrix_in_basis(A,V,basis_A,basis_V,phi)  # So we can compute their matrices
        P = matrix_in_basis(B,V,basis_B,basis_V,psi)

        if M.is_invertible() and P.is_invertible():
            go = False
            P_inv = P.inverse()
            for k in range(N):
                coord = P_inv * phi(basis_A[k]) #cooridnate of psi^-1(phi(basis_A[k])) of e in the basis of B
                new_isom_dict[k] = sum (coord[i]*basis_B[i] for i in range(N)) # = psi^-1(phi(basis_A[k]))
        n=n+1

    if go :
        return "Suitable v not found"

    if cpt == True :
        return new_isom_dict,n
    return new_isom_dict



##------ From rank one element in A ⊗ B^op -----------##


def explicit_iso_quat_alg(A,B):
    """
    Return an explicit isomorphism between A and B.
    """
    C = tensor(A,opposite(B))
    x = rank_one_MnQ(C)
    if x == "not found":
        return "rank one element in  A ⊗ B^op not found"

    isom_dict = matrix_ring_iso_from_rank_one(C, x)  # isomorphim phi : C -> M_4(Q) but given with key in 0,..,n²-1
    phi_dict = {}
    for i in range(4):
        for j in range(4):
            phi_dict[(i,j)] = isom_dict[4*i+j]  # construct the dictionnary of phi but with entry pairs (i,j)

    isom_A_to_B = quat_alg_iso_from_matrix_ring_iso(A,B,phi_dict)
    return isom_A_to_B