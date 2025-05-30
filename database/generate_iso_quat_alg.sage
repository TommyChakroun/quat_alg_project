load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("database/database_utilities.sage")

#------------------------------------------------------------------------------------------
#
#            GENERATE DATA BASE OF ISOMORPHISM OF QUATERNION ALGEBRAS
#
#------------------------------------------------------------------------------------------


##---List of functions---## 

# random_iso_quat_alg(F,a,b)
# generate_quaternion_iso_db(output_file="quaternion_isomorphisms_formated.txt")

def random_iso_quat_alg(F,a,b):
    """
    Return c,d in F* and an isomorphism f : (c,d | F) -> (a,b |F). 
    INPUT :
        -- F -- a field
        -- a -- in F*
        -- b -- in F*
    OUTPUT :
        -- c -- in F*
        -- d -- in F*
        -- isom_dict -- dictionnary with key 0,1,2,3 and value in (a,b |F) representing f : (c,d | F) -> (a,b |F) i.e 
                        f(1) = isom_dict[0], f(i) = isom_dict[1], f(j) = isom_dict[2], f(k) = isom_dict[3]
        -- M -- a matrix 3*3 representing in each mine  isom_dict[1],isom_dict[2],isom_dict[3] in the basis 1,i,j,k of (a,b | F)
    """

    # We make two change of basis of (a,b | F): 
    # (1,i,j,k) ----random---> (1,e1,e2,e3) ----normalize nrd ----> (1,new_i,new_j,new_k)

    A = QuaternionAlgebra(F,a,b)
    one_A,i,j,k = A.basis()

    B,P = quat_alg_mixed_with_table(F,a,b,transformation = True) # (a,b | F) view with a random basis as a FiniteDimensionaAlgebra

    # We have with 1,i,j,k the natural basis of  (a,b | F) 
    # We view this new basis really in A
    e0 = one_A
    e1 = P[0][0]*i + P[0][1]*j + P[0][2]*k
    e2 = P[1][0]*i + P[1][1]*j + P[1][2]*k
    e3 = P[2][0]*i + P[2][1]*j + P[2][2]*k

    new_i,new_j,c,d = quaternion_structure(B) 
    new_k = new_i*new_j
    
    new_i_vect = list(new_i.vector())
    new_j_vect = list(new_j.vector()) 
    new_k_vect = list(new_k.vector())

    rows = [new_i_vect[1:],new_j_vect[1:],new_k_vect[1:]]
    S = Matrix(F,rows)

    isom_dict = {}
    isom_dict[0] = one_A
    isom_dict[1] = new_i_vect[0]*e0 + new_i_vect[1]*e1 + new_i_vect[2]*e2 + new_i_vect[3]*e3
    isom_dict[2] = new_j_vect[0]*e0 + new_j_vect[1]*e1 + new_j_vect[2]*e2 + new_j_vect[3]*e3
    isom_dict[3] = new_k_vect[0]*e0 + new_k_vect[1]*e1 + new_k_vect[2]*e2 + new_k_vect[3]*e3

    return c,d,isom_dict,S*P


def generate_quaternion_iso_db(output_file="db_iso_quat_alg.txt"):
    """
    Generate database of quaternion algebra isomorphisms for all pairs (a,b) in [-10,10]*[-10,10]
    excluding zero values.
    The data base is a file db_iso_quat_alg.txt with blocks :
    c
    d
    a
    b
    (x1,y1,z1)
    (x2,y2,z2)
    (x3,y3,z3)

    which mean we have an isomorphism f : (c,d|F) ->(a,b | F) given by 
    f(1)=1
    f(i)=x1*i+y1*j+z1*k 
    f(j)=x2*i+y2*j+z2*k
    f(k)=x3*i+y3*j+z3*k
    """
    F = QQ  # Field of rational numbers
    
    # Generate all valid parameter values in the range [-10, 10] excluding 0
    values = list(range(-10, 0)) + list(range(1, 11))
    
    with open(output_file, "w") as f:
        examples_generated = 0
        
        # Generate all possible parameter pairs
        for a in values:
            for b in values:
                try:
                    # Generate a random isomorphism using your function
                    c, d, isom_dict,R= random_iso_quat_alg(F, a, b)
                    
                    # Format the quaternion algebra and isomorphism
                    f.write(f"{c}\n")
                    f.write(f"{d}\n")
                    f.write(f"{a}\n")
                    f.write(f"{b}\n")
                    f.write(f"["+str(list(R[0]))+","+str(list(R[1]))+","+str(list(R[2])) +"] \n" )
                    f.write("\n")  # Add empty line between examples
                    
                    examples_generated += 1
                    
                    # Print progress after every 100 examples
                    if examples_generated % 100 == 0:
                        print(f"Generated {examples_generated} examples")
                        
                except Exception as e:
                    print(f"Error generating example with parameters a={a}, b={b}: {e}")
                    continue
    
    print(f"Successfully generated {examples_generated} examples in {output_file}")