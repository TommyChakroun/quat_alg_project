load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("core/explicit_iso_matrix_ring.sage")
load("database/database_utilities.sage")

#------------------------------------------------------------------------------------------
#
#            TEST OF ISOMROPHISM phi: C -> M_N(F)
#
#------------------------------------------------------------------------------------------


##--- C = A ⊗ B^op  and we know f : A -> B isomorphism --- ##

def test_3():
    x = read_file(20,"/quat_alg_project/database/db_iso_matrix_ring.txt")
    for block in x:
        c,d,a,b = block[0],block[1],block[2],block[3]
        list_of_mat = [Matrix(QQ,block[i]) for i in range(4,20)] # 16 matrices
        A = QuaternionAlgebra(QQ,c,d)
        B = QuaternionAlgebra(QQ,a,b)

        C = tensor(A, opposite(B))
        isom_dict={}

        for k in range(16):
            isom_dict[k]=list_of_mat[k]
        T = MatrixSpace(QQ,4,4)
        
        print("iso : " + str(is_linear_iso(isom_dict,C,T)))
        print("alg homo : " +str(is_algebra_homomorphism(isom_dict, C, T)))

##--- C whatever and we know x in C of rank 1 ---#