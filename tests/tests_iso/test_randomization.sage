load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/isomorphism/explicit_iso_quat_alg.sage")

#------------------------------------------------------------------------------------------
#
#            TEST OF THE RANDOMIZATION IN THE ALGO FOR FINDING AN
#               ISOMORPHISM f:A->B FROM phi A Bop ->M_n²(F)
#------------------------------------------------------------------------------------------

## We use only this example from the data base

# If A = (-252,-36 | Q)) and B =(-9,-7 | Q)
# The following table give an isomorphism phi : A ⊗ B^op -> M_n²(Q)
# We want to deduce an isomrophism f : A -> B with the randomized algorithm.

c = -252
d = -36
a = -9
b = -7
list_of_list_of_list = [
[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
[[0, -9, 0, 0], [1, 0, 0, 0], [0, 0, 0, 9], [0, 0, -1, 0]],
[[0, 0, -7, 0], [0, 0, 0, -7], [1, 0, 0, 0], [0, 1, 0, 0]],
[[0, 0, 0, -63], [0, 0, 7, 0], [0, -9, 0, 0], [1, 0, 0, 0]],
[[0, 0, 0, -126], [0, 0, -14, 0], [0, 18, 0, 0], [2, 0, 0, 0]],
[[0, 0, 126, 0], [0, 0, 0, -126], [18, 0, 0, 0], [0, -18, 0, 0]],
[[0, -126, 0, 0], [-14, 0, 0, 0], [0, 0, 0, -126], [0, 0, -14, 0]],
[[-126, 0, 0, 0], [0, 126, 0, 0], [0, 0, 126, 0], [0, 0, 0, -126]],
[[0, -18, 0, 0], [2, 0, 0, 0], [0, 0, 0, -18], [0, 0, 2, 0]],
[[-18, 0, 0, 0], [0, -18, 0, 0], [0, 0, 18, 0], [0, 0, 0, 18]],
[[0, 0, 0, 126], [0, 0, -14, 0], [0, -18, 0, 0], [2, 0, 0, 0]],
[[0, 0, -126, 0], [0, 0, 0, -126], [-18, 0, 0, 0], [0, -18, 0, 0]],
[[0, 0, -252, 0], [0, 0, 0, 252], [36, 0, 0, 0], [0, -36, 0, 0]],
[[0, 0, 0, -2268], [0, 0, -252, 0], [0, -324, 0, 0], [-36, 0, 0, 0]],
[[-252, 0, 0, 0], [0, 252, 0, 0], [0, 0, -252, 0], [0, 0, 0, 252]],
[[0, 2268, 0, 0], [252, 0, 0, 0], [0, 0, 0, -2268], [0, 0, -252, 0]]
]


A = QuaternionAlgebra(QQ,c,d)
B = QuaternionAlgebra(QQ,a,b)

## Good format :

list_of_matrix = [Matrix(QQ,e) for e in list_of_list_of_list]

dict_of_matrix = {}
for i in range(4):
    for j in range(4):
        dict_of_matrix[(i,j)]=list_of_matrix[4*i+j]

def test_5(iteration):
    wait_dict = {}
    for n in range(iteration):
        correct = True
        isom_dict_A_to_B,k = quat_alg_iso_from_matrix_ring_iso(A,B,dict_of_matrix,cpt=True)
        if k in wait_dict:
            wait_dict[k] += 1
        else:
            wait_dict[k] = 1
        if not is_linear_iso(isom_dict_A_to_B,A,B):
            correct = False
        if not is_algebra_homomorphism(isom_dict_A_to_B,A,B):
            correct = False
    print("The algorithm is correct:", correct)
    print(f"Among the {iteration} iterations")
    print("the number of cases where we needed to make N random choices to find v is given by:")
    for key, value in wait_dict.items():
        print(f"{key} : {value}")