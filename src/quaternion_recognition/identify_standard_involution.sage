load("utilities/utilities.sage")
load("utilities/algebra_type.sage")

#------------------------------------------------------------------------------------------
#
#                   IDENTIFY STANDARD INVOLUTION
#
#------------------------------------------------------------------------------------------



def has_standard_involution(B):
    """
    INPUT : 
        --B-- finite dimensional algebra over F
    OUTPUT : 
        --boolean-- true if and only if A has a standard involution

    ALGORITHM : John Voight 2012 Identifying the Matrix Ring, Algorithm 3.12

    EXAMPLES :

    """
    if not B.is_associative():
        return False
    if not B.is_unitary():
        return False

    C = B.table()   
    e = B.basis()
    n = B.dimension()

    # we compute trd[i] = the reduced trace of ei if B has a standard involution
    t  = [0 for i in range(n)] 

    for i in range(1,n):
        t[i] = C[i][i][i] #coefficient of ei in ei**2
        ni = e[i]*e[i]-t[i]*e[i]  # potential reduced norm of ei
        if not is_in_field(ni,B):
            return False

    for i in range(1,n):
        for j in range(i+1,n):
            nij = (e[i]+e[j])**2-(t[i]+t[j])*(e[i]+e[j]) #potential reduced norm of ei+ej
            if not is_in_field(nij,B) :
                return False

    return True


def reduced_trace_si(B,x):
    """
    INPUT: 
        -- B -- an algebra over a field F, associative, with first basis element e0 neutral, assume to have a standard involution
        -- x -- an element of B
    OUTPUT: 
        -- trd(x) the reduced trace of x
    """

    if not has_standard_involution(B):
        raise TypeError("algebra has not a standard involution")

    C = B.table()   
    e = B.basis()
    n = B.dimension()

    t = [0 for i in range(n)] # we will compute t[i] = the reduced trace of ei
    t[0] = 2 # 1 of F
    for i in range(1,n):
        t[i] = C[i][i][i] #coefficient of ei in ei**2

    trd_x = 0
    for i in range(n):
        trd_x = trd_x + x.vector()[i]*t[i]
    return trd_x

    return trd


def standard_involution(B,x):
    """
    INPUT: 
        -- B -- an algebra over F, associative, with first basis element e0 neutral, assumed to have a standard involution
        -- x -- an element of B
    OUTPUT: 
        --  bar(x) the conjugate of x through the standard involution
    """
   
    if not has_standard_involution(B):
        raise TypeError("algebra has not a standard involution")

    e0 = B.basis()[0] # 1 of B

    return reduced_trace_si(B,x)-x


def reduced_norm(B,x):
    """
    INPUT: 
        -- B -- an algebra over F,associative, with first basis element e0 neutral,assumed to have a standard involution
        -- x -- an element of B
    OUTPUT: 
        -- nrd(x) the reduced norm of x
    """
    if not has_standard_involution(B):
        raise TypeError("algebra has not a standard involution")

    return x*standard_involution(B,x)
