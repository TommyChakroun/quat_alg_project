load("src/maximal_orders/maximal_orders_utilities.sage")
load("database/generate_data_base.sage")

#------------------------------------------------------------------------------------------
#
#            A list of global variable easy to use for examples and testing
#
#------------------------------------------------------------------------------------------


#####       1. ISOMORPHISM OF QUATERNION ALGEBRA     ######


ListIsoQuatAlgDb = read_file(5,"database/db_iso_quat_alg.txt")
ListIsoM4QDb_temp = read_file(20,"database/db_iso_matrix_ring.txt")
ListIsoM4QDb = [S[:4]+[Matrix(QQ,S[i]) for i in range(4,20)] for S in ListIsoM4QDb_temp ]







#####         2. Maximal Order in M4Q                ######



## M4(Q) and its basis 

MatriceRingM4Q = MatrixSpace(QQ,4,4) 
BasisM4Q = list(MatriceRingM4Q.basis())

## A list of invertible 16*16 rationals matrices

RadomInvertibleMatrix_16x16_example_1 = Matrix(QQ, [
    [   2,   -2,   -1,    0,    1,    1, -1/2,   -1,    0,  1/2,  1/2,    0,    2,  1/2,    1,    0],
    [   1,    0,    2,  1/2,    0,    0,    1,   -1,    0,    2,    0,    0,    0,   -1,    2,    1],
    [  -1,   -1,  1/2,    1,    0,    0,   -1,   -1,   -1,   -1,    0,   -1,    0,    1,    2,   -1],
    [   1,    2,    0,    1,    0,    0,   -1,    0,    1,    2,    0,   -2,    0,   -1,    0,    1],
    [   0,    1,    0,    0,   -1,    0,    0,    0,   -2, -1/2,   -1,   -1,    0, -1/2,   -1,   -2],
    [-1/2,    0,    1,    2,    0,  1/2, -1/2,    2,    0,    0,    1,    0,    0,    0,  1/2,   -2],
    [ 1/2,    0,    2,    1,    0,    2,    0,    2,    0,    1,    0,    0,    2,   -1,    2,    0],
    [   0,   -2,    2,  1/2,  1/2,    1,    1,   -1, -1/2,    2,    0,    0,    0,    0,    1,   -1],
    [   1,    0,   -1,    2,    2, -1/2,  1/2,   -1, -1/2,   -1,    0,  1/2, -1/2,    0,    1,    0],
    [   0,    0,    2,   -2,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0],
    [   2,    0, -1/2,   -2,    0,    1,    1,    2,   -1,    0, -1/2,   -2,   -2,    1,    0,   -1],
    [   0,   -2,  1/2,    0,    0,    1,    0, -1/2,   -2,   -2,    0,    0,   -2,   -1,    0,    1],
    [   1,    0,    0,   -1,   -1,    0,    1,    0,   -1,   -1,    1,    1,    0,    1,   -1,    0],
    [  -2,    1,   -1,    1,    0,    2,   -1,   -1,    0,   -2,    0,    0,  1/2,    0,    0,    0],
    [-1/2,   -1,    1,    1,  1/2,    0,    1,   -1,    0,    0,    0,    1,    1,    0,    0,    0],
    [   0,   -2, -1/2,    0,    0,    0,   -2,  1/2,   -1,  1/2,    1,    2,   -2,    0,   -1,    0]
])


RadomInvertibleMatrix_16x16_example_2 = Matrix(QQ, [
    [  -3,    1,  -10,    0,   -2,   -3,    7,   -1,   -1,    0,    0,    2,   13,   -1,   -1,    0],
    [  -1,    9,    1,   -2,   -2,    1,   -1,   -1,    0,    1,    1,    2,   13,   -1,    0,   -1],
    [   3,    0,   10,    0,   -1,    2,   -1,   -1,    1,   -1,   -1,   -8,    7,    0,    1,  -39],
    [  -1,   -6,    0,   -1,    2,    0,   -1,    1,   -1,    0,    1,   -2,   -2,  -12,   19,    7],
    [   0,   -2,   10,    4,    0,   -3,    3,   -1,    1,    2,    1,    0,    0,    0,   -1,   -9],
    [   1,   -1,    1,   -6,   -3,    1,    1,   -3,   -2,    2,   -1,    2,    4,   11,   -1,    0],
    [   3,    0,   -1,    1,   -1,   -2,    1,    1,    0,    3,   -3,   -6,    0,   -1,   -1,    1],
    [-291,   -1,    1,   -1,    1,    0,   -7,   -1,    0,    5,    1,    2,    0,   -4,   -1,   -1],
    [  -1,    1,    0,    1,   -1,   -2,    0,    0,   29,   -1,   -2,    1,   -1,    1,    0,   -1],
    [   2,    1,    0,   -2,    2,   -1,    1,   -1,   88,    2,    3,   -1,   -1,   -2,    1,    1],
    [  -1,    0,    3,   -2,    3,   -4,    1,    3,   -1,   -5,   -1,   -1,    1,  -43,    1,    0],
    [  -6,    1,    0,    0,    4,  -22,    1,    1,   -3,    1,   -2,   -4,  -11,    1,   -1,    0],
    [   1,    0,   -1,    0,   13,   -3,    0,    1,    4,   -1,    0,    1,    0,    1,    3,    1],
    [   1,   -1,    0,    2,   -1,    1,    0,    1,    9,   -3,   -2,   -1,    4,  -13,   -1,    1],
    [   1,    1,   -1,    0,    1,    4,   -1,    0,    2,    1,    1,    0,    1,    2,    0,    0],
    [   0,    0,   -1,    3,   -1,    0,    1,   -1,    1,    0,  -14,    1,    2,    9,    2,    0]
])


RadomInvertibleMatrix_16x16_example_3 = Matrix(QQ, [
    [  1,  0, -1,  1,  1,  3,  5,  0, -3, -5,  4,  1,  1, -5,  1,  5],
    [  0,  1, -4,  3,  4,  3,  0, -3, -3,  3,  4,  3, -3, -1, -3,  5],
    [  0,  0,  1,  4, -4,  3,  1,  1,  5, -2, -4,  5, -1,  1,  3,  1],
    [  0,  0,  0,  1, -3,  3, -1, -2, -3, -1, -5, -2,  4,  1,  4, -1],
    [  0,  0,  0,  0,  1, -1,  3, -4, -1, -2, -5, -2, -4,  2, -5, -3],
    [  0,  0,  0,  0,  0,  1, -5,  0, -3,  2, -1, -3,  2, -3,  1, -2],
    [  0,  0,  0,  0,  0,  0,  1,  1,  3, -5, -3,  5,  1,  2,  0,  0],
    [  0,  0,  0,  0,  0,  0,  0,  1,  0, -1, -3, -5, -2, -5, -4,  1],
    [  0,  0,  0,  0,  0,  0,  0,  0,  1,  5, -3,  5,  1,  3,  5, -3],
    [  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0, -5,  3, -3, -3],
    [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  1, -1,  5],
    [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -5, -1,  5,  3],
    [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -5,  5, -4],
    [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  2],
    [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2],
    [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1]
])




## A list of random basis of M4(Q) in order to view 

def basis_from_invertible_matrix(P):
    return [sum(P[i][k]*BasisM4Q[k] for k in range(16)) for i in range(16) ] 

RandomBasisM4Q_example_1 = basis_from_invertible_matrix(RadomInvertibleMatrix_16x16_example_1)
RandomBasisM4Q_example_2 = basis_from_invertible_matrix(RadomInvertibleMatrix_16x16_example_2)
RandomBasisM4Q_example_3 = basis_from_invertible_matrix(RadomInvertibleMatrix_16x16_example_3)


## A list of random Z basis of order of M4(Q) (construct respectivelu as the left order of the basis above)

#BasisRandomOrder_M4Q_example_1 = left_order(MatriceRingM4Q,RandomBasisM4Q_example_1)
#BasisRandomOrder_M4Q_example_2 = left_order(MatriceRingM4Q,RandomBasisM4Q_example_2)
#BasisRandomOrder_M4Q_example_3 = left_order(MatriceRingM4Q,RandomBasisM4Q_example_3)



## The abstract version, given by structure constant, of M4(Q) "view" with the random basis above

def structure_constants_from_new_basis(new_basis):
    table = []
    for e in new_basis:
        rows = []
        for f in new_basis:
            x = f*e
            rows.append(coordinate(x,MatriceRingM4Q,new_basis))
        M = Matrix(QQ,rows)
        table.append(M)
    return table


MixedM4Q_example_1 = FiniteDimensionalAlgebra(QQ,structure_constants_from_new_basis(RandomBasisM4Q_example_1))
MixedM4Q_example_2 = FiniteDimensionalAlgebra(QQ,structure_constants_from_new_basis(RandomBasisM4Q_example_2))
MixedM4Q_example_3 = FiniteDimensionalAlgebra(QQ,structure_constants_from_new_basis(RandomBasisM4Q_example_3))