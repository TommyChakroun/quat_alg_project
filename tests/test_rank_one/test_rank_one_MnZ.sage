load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/rank_one/rank_one_MnQ.sage")

#------------------------------------------------------------------------------------------
#
#            TEST RANK ONE IN Mn(Z)
#
#------------------------------------------------------------------------------------------


def test_17(a,b,c,d):
    """
    Tests the zero divisor algorithm on a tensor product of quaternion algebras.
    """
    # --- 1. Setup the Algebra ---
    A = QuaternionAlgebra(QQ,a,b)
    B = QuaternionAlgebra(QQ,c,d)
    C = tensor(A,opposite(B))
    BO = max_order(C)
    C = finite_dimensional_algebra_format(C,BO)

    # --- 2. Run the Zero Divisor Algorithm ---
    # Note: Corrected to run on the proper algebra 'C'.
    p = 5
    n_lifting_steps = 4
    try:
        e = zero_divisor_MnZ(C, p=p, n_lifting_steps=n_lifting_steps)
    except Exception as err:
        print(f"\n--- Test Failed for ({a}, {b}) ⊗ ({c}, {d}) ---")
        print(f"  ERROR: The function zero_divisor_MnZ failed with: {err}")
        return None

    # --- 3. Analyze the Results ---
    rank_e = right_rank(e,C)
    x = e**2-e
    rank_x = right_rank(x,C)
    coords_x = x.vector()

    # --- 4. Simple Printer ---
    # This section prints a concise summary of the test outcome.
    print(f"\n--- Test Result for Quat({a}, {b}) ⊗ Quat({c}, {d})^{{op}} ---")

    if rank_e < 16:
        print("✅ Success: Found a zero divisor.")
        print(f"  - Rank of element 'e': {rank_e} (out of 16)")

        if rank_x == 0:
            print("  - Element 'e' is a perfect idempotent (e^2 - e = 0).")
        else:
            print("  - Element 'e' is NOT a perfect idempotent.")
            print(f"  - Rank of error 'x = e^2 - e': {rank_x}")

            # Verify if the error term is small in the p-adic sense
            modulus_K = p**(2**n_lifting_steps)
            is_padically_small = all(coord % modulus_K == 0 for coord in coords_x)
            print(f"  - Error 'x' is p-adically small (coords % {modulus_K} == 0): {is_padically_small}")
    else:
        print("❌ Failure: Did not find a zero divisor.")
        print(f"  - Rank of element 'e': {rank_e} (full rank)")

    return None



def test_stat(A,ite,p,n):
    cpt = 0
    Op,_ = finite_algebra_from_order(A, A.basis(), p)
    dico = {0:0,4:0,8:0,12:0,16:0}
    for _ in range(ite):
        e = heuristic_zero_divisor_MnZ(A,primes=[p],n_lifting_steps = n,Op=Op)
        r = right_rank(e,A)
        dico[r] = dico[r]+1
    return dico

import time

def big_test_stat(A, ite, primes, n_values):
    """
    Repeats the test_stat function with different values of p and n,
    and prints the execution time for each combination.

    Args:
        A: The algebra object to be tested.
        ite (int): The number of iterations for each call to test_stat.
        primes (list): A list of prime numbers to be used as 'p'.
        n_values (list): A list of integer values for 'n' (n_lifting_steps).
    """
    print("Starting the comprehensive test.")
    print("-" * 40)

    for p in primes:
        print(f"Working with prime p = {p}")
        print("-" * 40)
        
        # It's generally more efficient to compute Op once for each prime
        try:
            Op, _ = finite_algebra_from_order(A, A.basis(), p)
        except Exception as e:
            print(f"Could not create finite algebra for p = {p}. Error: {e}")
            continue

        for n in n_values:
            print(f"  Testing with n = {n}")
            
            start_time = time.time()
            
            # Since test_stat is provided, we will call it.
            # We pass the pre-computed Op to avoid redundant calculations.
            # We modify the call to test_stat to accept Op.
            
            def test_stat_modified(A, ite, p, n, Op):
                cpt = 0
                for _ in range(ite):
                    e = heuristic_zero_divisor_MnZ(A, [p], n_lifting_steps=n, Op=Op)
                    r = right_rank(e, A)
                    if r < 16:
                        cpt = cpt + 1
                return cpt, ite, cpt / ite

    
            cpt, total_ite, ratio = test_stat_modified(A, ite, p, n, Op)
                
            end_time = time.time()
            elapsed_time = end_time - start_time
                
            print(f"    - test_stat run time: {elapsed_time:.4f} seconds")
            print(f"    - Results: Found {cpt} elements with rank < 16 out of {total_ite} iterations.")
            print(f"    - Ratio: {ratio}")

            
            print("-" * 20)
            
    print("Comprehensive test finished.")


def test_18(A,k,l = [-3,-2,-1,0,1,2,3]):
    dim_A = dimension(A)
    BA = A.basis()
    n = len(l)
    cpt = 0
    for _ in range(k):
        x = A.sum(l[randint(0,n-1)]*BA[i] for i in range(dim_A))
        pi = minimal_polynomial(x,A)
        if not pi.is_irreducible():
            cpt = cpt +1
    return cpt



