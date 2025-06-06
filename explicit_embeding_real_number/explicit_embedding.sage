load("utilities/utilities.sage")
load("minimal_ideals/minimal_ideals_manually.sage")

#------------------------------------------------------------------------------------------
#
#            EXPLICIT EMBEDING  A ->  A⊗R -> Mn(R)
#
#------------------------------------------------------------------------------------------

"""
Let A be a central simple algebra over Q of dimension N = n^2
Suppose that A is isomrophic to Mn(Q).
Then there exist alpha in R such that if K = Q(alpha) we can compute an explicit 
isomorphism
"""


def splitting_element(A,nb_ite = 10):
    """
    see https://www.ams.org/journals/mcom/1990-55-192/S0025-5718-1990-1035925-1/S0025-5718-1990-1035925-1.pdf
    Lemma 2.4.

    INPUT : 
        -- A -- central simple algebra over Q of dimension N=n^2
    OUTPUT :
        -- a -- spliiting element of A, that is the minimal polynomial pi of a over Q 
                is separable (gcd(pi,pi')=1 i.e. pi has simple root in bar(Q)) 
                and deg(pi)=n

    REMARK : In particular a is a cycliq matrix view in Mn(Q)
            not sure if that is useful.
    """

    N = dimension(A)
    BA = list(A.basis())
    n = int(sqrt(N))
    H = 2*n*(n-1)
    for ite in range(nb_ite):
        t = [randint(1,H) for l in range(N)]
        a = sum( t[i]*BA[i] for i in range(N))
        pi = minimal_polynomial(a,A)
        if pi.degree()==n and gcd(pi,pi.derivative())==1:
            return a
    return "Not found"