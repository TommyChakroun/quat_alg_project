load("utilities/utilities.sage")
load("database/database_utilities.sage")
load("core/explicit_iso_matrix_ring.sage")

#------------------------------------------------------------------------------------------
#
#            GENERATE DATA BASE OF ISOMORPHISM  phi : A ⊗ Bop -> Mn²(Q)
#
#------------------------------------------------------------------------------------------

##---List of functions---## 

# generate_matrix_ring_iso_db(output_file="db_iso_matrix_ring.txt")


##-- Convert each isomorphism f:A->B in an isomorphism phi : A ⊗ Bop -> Mn²(Q) --## 

def generate_matrix_ring_iso_db(output_file="db_iso_matrix_ring.txt"):
    """
    Generate the database of isomorphisms phi : A ⊗ Bop -> Mn²(Q).
    The data base is a file db_iso_matrix_ring.txt with blocks :
    c
    d
    a
    b
    [[..]],[..],..,[..]]
    [[..]],[..],..,[..]]
    .
    .
    [[..]],[..],..,[..]]
    [[..]],[..],..,[..]]

    which mean that for A = (a,b | F) and B = (c,d | F) we have an isomosphism 
    phi : A ⊗ Bop -> Mn²(Q)
    given by 
    phi(A.basis()[j],B.basis()[j]) = the matrix of the n²*i+j line 

    """
    F = QQ  
    list_blocks = read_file(5,"/quat_alg_project/database/db_iso_quat_alg.txt")
    
    with open(output_file, "w") as f:
        examples_generated = 0
        
        # Generate all possible parameter pairs
        for c,d,a,b,R in list_blocks:
            A = QuaternionAlgebra(F,c,d)
            B = QuaternionAlgebra(F,a,b)
            one_B,i,j,k = B.basis()

            isom_dict= {}
            isom_dict[0]=one_B
            isom_dict[1]=R[0][0]*i+R[0][1]*j+R[0][2]*k
            isom_dict[2]=R[1][0]*i+R[1][1]*j+R[1][2]*k
            isom_dict[3]=R[2][0]*i+R[2][1]*j+R[2][2]*k

            tensor_isom_dict = matrix_ring_iso_from_algebra_iso(A,B,isom_dict)
                    
            # Format the quaternion algebra and isomorphism
            f.write(f"{c}\n")
            f.write(f"{d}\n")
            f.write(f"{a}\n")
            f.write(f"{b}\n")
            
            for key in tensor_isom_dict:
                list_of_tuple = tensor_isom_dict[key].rows()
                list_of_list = []
                for tuple in list_of_tuple:
                    list_of_list.append(list(tuple))
                f.write(f"{list_of_list}\n")
            f.write("\n")  # Add empty line between examples

            examples_generated += 1
                    
            # Print progress after every 100 examples
            if examples_generated % 100 == 0:
                print(f"Generated {examples_generated} examples")
    
    print(f"Successfully generated {examples_generated} examples in {output_file}")
