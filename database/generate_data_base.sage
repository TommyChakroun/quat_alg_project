load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/quaternion_recognition/identify_quaternion_algebra.sage")
load("src/iso_splitting_algebra/explicit_iso_matrix_ring.sage")

#------------------------------------------------------------------------------------------
#
#            GENERATE DATA BASE OF ISOMORPHISM OF QUATERNION ALGEBRAS 
#            and ISOMORPHISM phi : A ⊗ Bop -> M4(Q)
#
#------------------------------------------------------------------------------------------


##---List of functions---## 

# random_iso_quat_alg(F,a,b)
# generate_quaternion_iso_db(output_file="quaternion_isomorphisms_formated.txt")


## --- UTILITIES --- #

import ast

def read_line(string):
    """
    Parse a line:
    - If it's an integer, return int.
    - If it's a list of lists, parse manually and return as list of list of QQ.
    """
    s = string.strip()

    # Case: integer
    if s[0] != "[":
        return QQ(s)

    import re
    try:
        # Match all lists inside the outer list using regex
        row_strings = re.findall(r'\[([^\[\]]+)\]', s)
        matrix = []
        for row_str in row_strings:
            entries = [QQ(entry.strip()) for entry in row_str.split(',')]
            matrix.append(entries)
        return matrix
    except Exception as e:
        raise ValueError(f"Could not parse matrix string: {e}")


def read_one_block(lines):
    """
    Reads a block starting at start_index.
    The block alternates: int, list_of_lists, ..., and may end in either one.
    Returns (block_data, next_index).
    """
    return [read_line(line) for line in lines if line.strip()]


def read_file(size_block,filename):
    """
    Return list_of_blocks the list of each block
    """
    list_blocks = []
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]  # skip empty lines
    
    i = 0
    while i + size_block+1 < len(lines):
        block = []
        for k in range(size_block):
            block.append(lines[i+k])
        list_blocks.append(read_one_block(block))
        i = i+size_block
    return list_blocks


## --- Generate data base of isomorphism f : A -> B --- ##

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


def generate_quaternion_iso_db(output_file="database/db_iso_quat_alg.txt"):
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
    values = list(range(-50, 0)) + list(range(1, 51))
    
    with open(output_file, "w") as f:
        examples_generated = 0
        
        # Generate all possible parameter pairs
        for N in range(400):
            a = values[randint(0,98)]
            b = values[randint(0,98)]
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



## --- Generate data base of isomorphism phi : A ⊗ Bop -> M4(Q) from the date base of f : A -> B  --- ##

def generate_matrix_ring_iso_db(output_file="database/db_iso_matrix_ring.txt"):
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
    list_blocks = read_file(5,"database/db_iso_quat_alg.txt")
    
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



## --- Regenerate Data Base --- ##

def regenerate_data_base():
    generate_quaternion_iso_db(output_file="database/db_iso_quat_alg.txt")
    generate_matrix_ring_iso_db(output_file="database/db_iso_matrix_ring.txt")
    print("Please load again load('load_all_file'.sage) to construct the global variables ListIsoQuatAlgeDb and ListIsoM4QDb")
    return None