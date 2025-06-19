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