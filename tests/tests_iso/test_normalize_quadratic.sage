load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/quaternion_recognition/normalize_quadratic.sage")

#------------------------------------------------------------------------------------------
#
#            TEST OF NORMALIZE QUADRATIC
#
#------------------------------------------------------------------------------------------


R = QQ  

matrices = []

matrices.append(Matrix(R, [
    [1, 2],
    [2, 1]
]))

matrices.append(Matrix(R, [
    [1, 1],
    [1, 1]
]))

matrices.append(Matrix(R, [
    [0, 1],
    [1, 0]
]))

matrices.append(Matrix(R, [
    [2, 1],
    [1, 4]
]))


matrices.append(Matrix(R, [
    [1, 2, 3, 4],
    [2, 5, 6, 7],
    [3, 6, 8, 9],
    [4, 7, 9, 10]
]))

matrices.append(Matrix(R, [
    [2, 0, 0, 0],
    [0, 4, 1, 1],
    [0, 1, 6, 2],
    [0, 1, 2, 8]
]))

matrices.append(Matrix(R, [
    [3, -1, 0, 1],
    [-1, 3, -1, 2],
    [0, -1, 3, -2],
    [1, 2, -2, 3]
]))

matrices.append(Matrix(R, [
    [10, 5, 0, 2],
    [5, 11, 1, 3],
    [0, 1, 12, 4],
    [2, 3, 4, 13]
]))

matrices.append(Matrix(R, [
    [0, 1, 2, 3],
    [1, 0, 4, 5],
    [2, 4, 0, 6],
    [3, 5, 6, 0]
]))

matrices.append(Matrix(R, [
    [1, 0, 1, 0],
    [0, 2, 0, 1],
    [1, 0, 3, 0],
    [0, 1, 0, 4]
]))

matrices.append(Matrix(R, [
    [5, 3, 2, 1],
    [3, 6, 4, 2],
    [2, 4, 7, 5],
    [1, 2, 5, 8]
]))

matrices.append(Matrix(R, [
    [7, 1, 0, -1],
    [1, 9, 1, 0],
    [0, 1, 11, 1],
    [-1, 0, 1, 13]
]))

matrices.append(Matrix(R, [
    [2, 2, 2, 2],
    [2, 2, 2, 2],
    [2, 2, 2, 2],
    [2, 2, 2, 2]
]))

matrices.append(Matrix(R, [
    [1, 1, 1, 1],
    [1, 2, 2, 2],
    [1, 2, 3, 3],
    [1, 2, 3, 4]
]))


def all_tests():
    for S in matrices:
        D,P = normalize_symetric_matrix_over_Z_loc_2(S)
        print("S = ")
        print(S)
        print("D = ")
        print(D)
        print("P = ")
        print(P)
        print("P^t S P = D : " + str(P.transpose()*S*P ==D))
        print("---------------")