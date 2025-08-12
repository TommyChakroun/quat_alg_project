# We need to import 'choice' for the random matrix generation
from random import choice

#------------------------------------------------------------------------------------------
#
#            RANDOM QUADRATIC SIMPLIFICATION
#
#------------------------------------------------------------------------------------------

# We do some experiment on the following.
# Starting from a symmetric matrix S, we output a "random" S' = diag(-c,-d,cd)
# which is equivalent to S over the rational numbers.
#
# Equivalence means there exists an invertible matrix M in GL(3,Q)
# such that M^T * S * M = S'.



def random_equivalent_quadratic(S, transformation=False):
    """
    Takes a 3x3 symmetric matrix S over QQ with a positive square determinant
    and returns a "random" equivalent diagonal matrix of the form diag(-c,-d,c*d).
    """
    d = S.det()
    if d == 0 or d < 0 or not is_square(d):
        raise ValueError("The determinant of S must be a non-zero positive square.")

    coeffs = [-2, -1, -1/2, 1/2, 1, 3/2, 2]
    M = matrix(QQ, 3, 3)
    while M.det() == 0:
        M = matrix(QQ, 3, 3, [choice(coeffs) for _ in range(9)])

    S1 = M.transpose() * S * M
    QF = QuadraticForm(S1)
    Q2, N = QF.rational_diagonal_form(return_matrix = True)
    S2 = Q2.matrix()

    d1, d2, d3 = S2.diagonal() 
    c = -d1
    d = -d2

    k = sqrt(S2.det())
    P = diagonal_matrix([1, 1, (c * d) / k])
    S_final = diagonal_matrix([-c, -d, c * d])

    if transformation:
        M_total = M * N * P
        return S_final, M_total
    else:
        return S_final