load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/isomorphism/explicit_iso_quat_alg.sage")

#------------------------------------------------------------------------------------------
#
#            TEST OF ISOMROPHISM OF QUATERNION ALGEBRAS f : A -B
#
#------------------------------------------------------------------------------------------


##--- When know phi : A ⊗ B^op -> M_n²(F) ---#

def test_4(output_file = "db_iso_quat_alg_from_tensor_iso"):
    x = read_file(20,"/quat_alg_project/database/db_iso_matrix_ring.txt")

    with open(output_file, "w") as f:
        iso_generated = 0
        for block in x:
            c,d,a,b = block[0],block[1],block[2],block[3]
            list_of_mat = [Matrix(QQ,block[i]) for i in range(4,20)] # 16 matrices
            A = QuaternionAlgebra(QQ,c,d)
            B = QuaternionAlgebra(QQ,a,b)

            C = tensor(A, opposite(B))
            isom_dict={}

            for i in range(4):
                for j in range(4):
                    isom_dict[(i,j)]=list_of_mat[4*i+j]

            new_isom_A_to_B = quat_alg_iso_from_matrix_ring_iso(A,B,isom_dict)

            R = [get_coefficients(new_isom_A_to_B[1],B),get_coefficients(new_isom_A_to_B[2],B),get_coefficients(new_isom_A_to_B[3],B)]

            # Format the quaternion algebra and isomorphism
            f.write(f"{c}\n")
            f.write(f"{d}\n")
            f.write(f"{a}\n")
            f.write(f"{b}\n")
            f.write(f"["+str(list(R[0][1:]))+","+str(list(R[1])[1:])+","+str(list(R[2])[1:]) +"] \n" )
            f.write("\n")  # Add empty line between examples
                        
            iso_generated += 1
                        
            # Print progress after every 100 examples
            if iso_generated % 100 == 0:
                print(f"Generated {iso_generated} examples")
    
        print(f"Successfully compute again {iso_generated} isomrophisms in {output_file}")