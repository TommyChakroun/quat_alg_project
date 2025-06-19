load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/quaternion_recognition/identify_standard_involution.sage")

#------------------------------------------------------------------------------------------
#
#                  IDENTIFY QUATERNION ALGEBRA
#
#------------------------------------------------------------------------------------------

##---List of functions---##

# is_isomorphic_to_quaternion_algebra(B)
# quaternion_structure(B)


def is_isomorphic_to_quaternion_algebra(B):
    """
    INPUT : 
        -- B -- a finite dimensional algebra over a field F
    OUTPUT : 
        -- boolean -- True if B is isomorphic to a quaternion algebra, False otherwise
    """
    if not has_standard_involution(B):
        return False

    C = B.table()   
    e = B.basis()
    n = B.dimension()

    if n != 4 :
        return False

    # S = matrice (trd(ei ej_bar)) de la forme quadratic nrd
    S = Matrix([[reduced_trace(B,e[i]*standard_involution(B,e[j])) for i in range(4)] for j in range(4)])

    return S.determinant() != 0


def quaternion_structure(B):
    """
    INPUT : 
        -- B -- a finite dimensional algebra over a field F which is isomrophic to a quaternion algebra
    OUTPUT : 
        -- i,j,a,b -- such that i,j in B and a,b in F* with 1,i,j,ij basis of B and i**=a, j**2=b, i*j=-j*i
    """
    if not is_isomorphic_to_quaternion_algebra(B):
        raise TabError("is not isomrophic to a quternion algebra")

    e = B.basis()
    F = B.base_ring()

    # S = matrice (trd(ei ej_bar)) matric de Gram de la forme quadratic 2*nrd
    S = Matrix(F,[[reduced_trace(B,e[i]*standard_involution(B,e[j])) for i in range(4)] for j in range(4)])

    nrd = QuadraticForm(S)

    P = nrd.rational_diagonal_form(return_matrix = True)[1]
    D = P.transpose()*S*P

    
    i = P[0][1]*e[0]+P[1][1]*e[1]+P[2][1]*e[2]+P[3][1]*e[3]
    j = P[0][2]*e[0]+P[1][2]*e[1]+P[2][2]*e[2]+P[3][2]*e[3]

    a = -D[1][1]/2 # -nrd(i)
    b = -D[2][2]/2 # -nrd(j)

    return i,j,a,b

