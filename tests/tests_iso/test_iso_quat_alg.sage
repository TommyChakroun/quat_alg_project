load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/isomorphism/explicit_iso_quat_alg.sage")
load("src/rank_one/rank_one_MnQ.sage")
load("src/isomorphism/explicit_iso_matrix_ring.sage")

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




import time

def test_19(a, b, c, d, p, n,nb_tries=100):
    """
    a,b,c,d integers representing the quat alg A = (a,b | Q) and B = (c,d | Q)
    p a prime
    n >0 an integer

    This test furnishes a series of prints and timers detailing all the process starting from
    A = (a,b | Q) and B = (c,d | Q) two quaternion algebras (assumed isomorphic)
    to the computation of an explicit isomorphism.
    """

    A = QuaternionAlgebra(a, b)
    B = QuaternionAlgebra(c, d)

    print("\n========================================")
    print("EXPLANATION")
    print("========================================\n")
    print(f"We consider two quaternion algebras over Q:")
    print(f"A = ({a}, {b} | Q),   B = ({c}, {d} | Q)")
    print("We are going to search for an explicit isomorphism A -> B.\n")
    print("Basic invariants:")
    print(" - Ramified primes of A:", A.ramified_primes())
    print(" - Ramified primes of B:", B.ramified_primes())
    print(" - Are A and B isomorphic?", A.is_isomorphic(B))

    print("\n========================================")
    print("STEP 1 - Compute maximal order in A")
    print("========================================\n")
    print("We compute a Z-basis of a maximal order in A.")
    t0 = time.time()
    BOA = A.maximal_order().basis()
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))

    print("========================================")
    print("STEP 2 - Compute maximal order in B")
    print("========================================\n")
    print("We compute a Z-basis of a maximal order in B.")
    t0 = time.time()
    BOB = B.maximal_order().basis()
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))

    print("========================================")
    print("STEP 3 - Tensor product A ⊗ B^op")
    print("========================================\n")
    print("We form the algebra C = A ⊗ B^op.")
    t0 = time.time()
    C_temp = tensor(A, opposite(B))
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))

    print("========================================")
    print("STEP 4 - Maximal order in A ⊗ B^op")
    print("========================================\n")
    print("We use the known maximal orders of A and B to deduce one in A ⊗ B^op.")
    t0 = time.time()
    BOC = max_order_tensor_quat_alg(A, B, C=C_temp, Zbasis_O1=BOA, Zbasis_O2=BOB)
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))

    print("========================================")
    print("STEP 5 - Format C with new basis")
    print("========================================\n")
    print("We rewrite C in terms of the basis of its maximal order.")
    t0 = time.time()
    C = finite_dimensional_algebra_format(C_temp, BOC)
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))

    print("========================================")
    print("STEP 6 - Compute F_p algebra")
    print("========================================\n")
    print(f"We compute a reduction modulo p={p} of the maximal order in C.")
    t0 = time.time()
    Op, _ = finite_algebra_from_order(C, C.basis(), p)
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))

    print("========================================")
    print("STEP 7 - Heuristic zero divisor search")
    print("========================================\n")
    print(f"We start {nb_tries} heuristic tries, working with p = {p} and precision p^(2^{n}).")
    t0 = time.time()
    dico = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    last_zero_div = None
    nb_reducible = 0
    for _ in range(nb_tries):
        e = heuristic_zero_divisor_MnZ(C, primes=[p], n_lifting_steps=n, Op=Op)
        r = right_rank(e, C) // 4
        dico[r] += 1
        if r==1 or r==3:
            last_zero_div = e
        if r==2 and last_zero_div==None:
            last_zero_div = e

        pi = minimal_polynomial(e,C)
        if not pi.is_irreducible():
            nb_reducible = nb_reducible+1
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds".format(t1 - t0))
    print(f"\nOn {nb_tries} heuristic tries we found:")
    for rank in [0,1,2,3,4]:
        print(f" - {dico[rank]} element(s) of rank {rank}")
    print(f" - {nb_reducible} element(s) with reducible minimal polynomial")
    if dico[1] == 0 and dico[2] == 0 and dico[3] == 0:
        print("We haven't found any zero divisor.")
        return None

    print("\n========================================")
    print("STEP 8 - Deduce rank one element in C")
    print("========================================\n")
    print(f"We use the last zero divisor found, with if possible rank 1 or 3, which has rank {right_rank(last_zero_div, C)//4}.")
    t0 = time.time()
    x = rank_one_MnQ(C, zero_divisor=last_zero_div)
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))

    print("========================================")
    print("STEP 9 - Write element in original basis")
    print("========================================\n")
    print("We express the rank one element in the original basis of C.")
    t0 = time.time()
    coords = x.vector()
    x = sum(coords[i] * BOC[i] for i in range(16))
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))

    print("========================================")
    print("STEP 10 - Isomorphism C -> M4(Q)")
    print("========================================\n")
    print("We compute an isomorphism from C to M4(Q) using this rank one element.")
    t0 = time.time()
    isom_list = matrix_ring_iso_from_rank_one(C_temp, x)
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))

    print("========================================")
    print("STEP 11 - Recover isomorphism A -> B")
    print("========================================\n")
    print("We recover an isomorphism A -> B from the matrix ring isomorphism.")
    t0 = time.time()
    isom_dict = {(i, j): isom_list[4 * i + j] for i in range(4) for j in range(4)}
    isom_A_B = quat_alg_iso_from_matrix_ring_iso(A, B, isom_dict)
    t1 = time.time()
    print("Done.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))


    print("\n========================================")
    print("STEP 11 - Check")
    print("========================================\n")
    print("Check that A -> B is indeed an algebra isomorphism.")
    t0 = time.time()
    assert is_linear_iso(isom_A_B, A, B), "The map is not Q-linear."
    assert is_algebra_homomorphism(isom_A_B, A, B), "The map is not a ring homomorphism."
    t1 = time.time()
    print("Check passed.")
    print("Time: {:.2f} seconds\n".format(t1 - t0))



    print("\n========================================")
    print("STEP 12 - Conclusion")
    print("========================================\n")
    print("The isomorphism A -> B is given by:")
    print("i ->", isom_A_B[1])
    print("j ->", isom_A_B[2])
    print("k ->", isom_A_B[3])

    return None




