load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("src/isomorphism/explicit_iso_quat_alg_equations.sage")

#------------------------------------------------------------------------------------------
#
#            TEST OF SOLVING QUADRATICS
#
#------------------------------------------------------------------------------------------


A0 = QuaternionAlgebra(QQ,-1, -1594331)
B0 = QuaternionAlgebra(QQ,-20726307, -22877028415373)

A1 = QuaternionAlgebra(QQ,-1, -1594331)
B1 = QuaternionAlgebra(QQ,-7233482772, -31495065018040921275)

A2 = QuaternionAlgebra(QQ,-1, -1594331)
B2 = QuaternionAlgebra(QQ,-247775045735, -226277280128968267100)

A3 = QuaternionAlgebra(QQ,-1, -1594331)
B3 = QuaternionAlgebra(QQ,-43920971776735, -164154462096737670211116900)


import time
def test_20(N, p):
    """
    Tests and compares isomorphism computation between two quaternion algebras
    using both SageMath and Magma. Assumes helper functions and imports
    are provided in the execution context.
    """

    A = QuaternionAlgebra(QQ, [p], [1/2])
    i, j, k = A.gens()

    alpha = -i.reduced_norm()
    beta = -j.reduced_norm()

    mu = randint(1, N) * i + randint(1, N) * j + randint(1, N) * k
    nu = quadratic_complement(A, mu) 

    a_B = -mu.reduced_norm()
    b_B = -nu.reduced_norm()

    if b_B.denominator() != 1:
        d = b_B.denominator()
        b_B = b_B * d^2
    
    a_B = ZZ(a_B)
    b_B = ZZ(b_B)

    B = QuaternionAlgebra(QQ, a_B, b_B)

    print("======================================================")
    print("                      SET UP")
    print("======================================================")
    print(f"A = QuaternionAlgebra(QQ, {alpha}, {beta})")
    print(f"B = QuaternionAlgebra(QQ, {a_B}, {b_B})")
    print(f"\nSanity check: A.is_isomorphic(B) returns: {A.is_isomorphic(B)}")
    print("------------------------------------------------------\n")


    print("======================================================")
    print("     COMPUTE ISOMORPHISM f: A -> B WITH SAGEMATH")
    print("======================================================")
    
    # Time the execution using time.time()
    t0_sage = time.time()
    v, iso_sage = iso_quat_alg(A, B) # Calling your function
    t1_sage = time.time()
    
    print(f"\nSage computation took: {t1_sage - t0_sage:.4f} seconds.")
    print("------------------------------------------------------\n")

    # --- Step 5: Compute with Magma ---
    print("======================================================")
    print("      COMPUTE ISOMORPHISM f: A -> B WITH MAGMA via norm equation")
    print("======================================================")
    
    magma_code = f"""
    Q := Rationals();
    A := QuaternionAlgebra<Q | {alpha}, {beta}>;
    B := QuaternionAlgebra<Q | {a_B}, {b_B}>;
    t0 := Cputime();
    is_iso, iso := IsIsomorphic(A, B : Isomorphism:=true);
    time_taken := Cputime(t0);
    printf "Magma internal computation took: %o seconds (CPU time).\\n", time_taken;
    printf "--- Magma Isomorphism Map ---\\n";
    iso;
    """
    
    t0_magma = time.time()
    magma_output = magma.eval(magma_code)
    t1_magma = time.time()
        
    print(f"Total time for Magma call: {t1_magma - t0_magma:.4f} seconds.")

    print("------------------------------------------------------\n")

    print("======================================================")
    print("  COMPUTE ISOMORPHISM f: A -> B WITH MAGMA quadratic of dim 4")
    print("======================================================")
    
    magma_code = f"""
    Q := Rationals();
    a := {a_B};
    b := {b_B};
    alpha := {alpha};
    beta := {beta};
    A := QuaternionAlgebra<Q | alpha, beta>;
    B := QuaternionAlgebra<Q | a, b>;
    D := DiagonalMatrix(Q, [a, b, -a*b, -alpha]);
    sol := Basis(IsotropicSubspace(D))[1];
    mu := sol[1]/sol[4]*B.1 + sol[2]/sol[4]*B.2 + sol[3]/sol[4]*B.3;
    nu := B ! [0, sol[1], sol[2], sol[3]];
    gamma := -Norm(nu);
    C := Conic([Q|1, -alpha, -beta/gamma]);
    v, point := HasRationalPoint(C);
    """
    
    t0_magma = time.time()
    magma_output = magma.eval(magma_code)
    t1_magma = time.time()
        
    print(f"Total time for Magma call: {t1_magma - t0_magma:.4f} seconds.")

    print("------------------------------------------------------\n")

    # --- Step 6: Conclusion ---
    print("======================================================")
    print("                     CONCLUSION")
    print("======================================================")
    print("Test complete.")
    print("======================================================")




import time

def iso_quat_alg_printers(A, B):
    """
    Finds an explicit isomorphism between two quaternion algebras A and B.

    ALGORITHM:
    1. Handle the case where A and B are matrix rings.
    2. To embed A = (alpha,beta | K) into B = (a,b | K), we first find an element
       `mu` in B such that `mu^2 = alpha`. This is the image of `i_A`. This step
       requires finding a rational point on the conic `a*Y^2 + b*Z^2 - a*b*W^2 + alpha = 0`.
    3. We then find an element `nu_0` in B that anticommutes with `mu`.
    4. We find the image of `j_A`, which will be of the form `j' = (x + y*mu)*nu_0`.
       We solve for `x, y` such that `(j')^2 = beta`. This reduces to solving the
       norm equation `x^2 - alpha*y^2 = beta / nu_0^2`, which is done by finding a point
       on the conic `X^2 - alpha*Y^2 - (beta/nu_0^2)*Z^2 = 0`.
    5. The isomorphism is then defined by `i_A -> mu` and `j_A -> j'`.

    INPUT:
        - ``A``, ``B`` -- Two quaternion algebras over the same number field K.

    OUTPUT:
        - A tuple ``(is_isomorphic, isom_map)`` where:
            - ``is_isomorphic`` is `True` if an isomorphism exists, `False` otherwise.
            - ``isom_map`` is a list of four elements in B which are the images
              of the basis elements 1, i, j, k of A.
    """
    if A.base_ring()!=QQ or B.base_ring()!=QQ:
        raise ValueError("Algebras must be over rational field")
    ram_A = A.ramified_primes()
    ram_B = B.ramified_primes()
    if ram_A != ram_B:
        return False,[]
    if ram_A ==[]:
        raise ValueError("A,B are split, use the other function")

    alpha, beta = A.invariants()
    a, b = B.invariants()
    i_A, j_A, k_A = A.gens()
    i_B, j_B, k_B = B.gens()

    print("=====================================")
    print("PRESENTATION")
    print("=====================================")
    print(f"A = ({alpha},{beta}|Q)")
    print(f"B = ({a},{b}|Q)")
    print()
    print(f"OUR NOTATION: ALPHA = {alpha}, BETA = {beta}, A = {a}, B = {b}")
    print("A AND B ARE ISOMORPHIC.")
    print("WE START COMPUTING AN ISOMORPHISM FROM A TO B.")
    print()
    print()

    print("=====================================")
    print("SECTION 1: PREPARATION, FACTORIZATION")
    print("=====================================")
    start_time = time.time()
    Falpha = factor(alpha)
    print(f"FACTORIZATION OF ALPHA TOOK: {time.time() - start_time:.6f} SECONDS")
    
    start_time = time.time()
    Fbeta = factor(beta)
    print(f"FACTORIZATION OF BETA TOOK: {time.time() - start_time:.6f} SECONDS")

    start_time = time.time()
    Fa = factor(a)
    print(f"FACTORIZATION OF A TOOK: {time.time() - start_time:.6f} SECONDS")

    start_time = time.time()
    Fb = factor(b)
    print(f"FACTORIZATION OF B TOOK: {time.time() - start_time:.6f} SECONDS")
    print()
    print()

    print("=====================================")
    print("SECTION 2: SEEK MU IN B SUCH THAT MU^2 = ALPHA")
    print("=====================================")
    print("I.E SOLVE IN Q, a*X^2 + b*Y^2 - a*b*Z^2 - alpha*W^2 = 0")
    start_time = time.time()
    sol = diagonal_qfsolve([a,b,-a*b,-alpha],factors =[Fa,Fb,Fa*Fb,Falpha])
    print(f"SOLVING IS DONE, IT TOOK: {time.time() - start_time:.6f} SECONDS")
    mu = sol[0]/sol[3]*i_B+sol[1]/sol[3]*j_B+sol[2]/sol[3]*k_B
    print()
    print()

    print("=====================================")
    print("SECTION 3: COMPUTE A QUATERNION COMPLEMENT TO MU")
    print("=====================================")
    start_time = time.time()
    nu = quaternionic_complement(B, mu)
    gamma = -nu.reduced_norm() 
    r,s = gamma.numerator(),gamma.denominator() 
    print(f"COMPUTATION TOOK: {time.time() - start_time:.6f} SECONDS")
    print(f"NU = {nu}")
    print(f"GAMMA = -NU.REDUCED_NORM() = {gamma}")
    print()
    print()

    print("=====================================")
    print("SECTION 4: FACTORIZE GAMMA")
    print("=====================================")
    start_time = time.time()
    Fr = factor(abs(r))
    Fs = factor(abs(s))
    print(f"FACTORIZATION OF GAMMA TOOK: {time.time() - start_time:.6f} SECONDS")
    print()
    print()

    print("=====================================")
    print("SECTION 5: MODIFY 1, MU, NU, MU*NU INTO THE FINAL BASIS OF B")
    print("=====================================")
    print("WE MUST SOLVE X^2 - alpha*Y^2 - (beta/gamma)*Z^2 = 0")
    start_time = time.time()
    sol2 = diagonal_qfsolve([r,-alpha*r,-beta*s],factors = [Fr,Falpha*Fr,Fbeta*Fs])
    print(f"SOLVING TOOK: {time.time() - start_time:.6f} SECONDS")
    x = sol2[0] / sol2[2]
    y = sol2[1] / sol2[2]
    print()
    print()

    isom_i = mu
    isom_j = (x + y*mu) * nu
    isom_map = [B.one(), isom_i, isom_j, isom_i * isom_j]

    return True, isom_map




def quaternionic_complement_special(B, mu,LLL = True):
    """
    Finds an element nu in B that anticommutes with a given pure quaternion mu.

    Given mu in B with mu^2 in K, the set of elements nu satisfying
    mu*nu + nu*mu = 0 forms a 2-dimensional K-subspace of the pure quaternions in B.
    This function finds one non-zero element in that space.

    INPUT:
        - ``B`` -- A quaternion algebra (a,b | K).
        - ``mu`` -- A pure quaternion in B (i.e., trace is 0).

    OUTPUT:
        - ``nu`` -- A non-zero pure quaternion in B such that nu*mu = -mu*nu.
    """
    if mu.reduced_trace() != 0:
        raise ValueError("Input element mu must be a pure quaternion.")

    (i, j, k), (a, b) = B.gens(), B.invariants()
    x0,x1,x2,x3 = mu.coefficient_tuple()


    if x1 != 0:
        v1 = (b*x2,-a*x1,0)
        v2 = (b*x3,0,x1)
    elif x2 != 0:
        v1 = (b*x2,-a*x1,0)
        v2 = (0,a*x3,x2)
    elif x3!=0 :
        v1 = (1,0,0)
        v2 = (0,1,0)

    if not LLL:
        return v1[0]*i+v1[1]*j+v1[2]*k

    M = Matrix(QQ,[v1,v2])
    d = M.denominator()
    dM = d*M
    N = M.LLL()

    nu1 = N[0][0]*i + N[0][1]*j + N[0][2]*k
    nu2 = N[1][0]*i + N[1][1]*j + N[1][2]*k

    return nu1


def size_of_gamma(A, B,LLL = True):
    if A.base_ring()!=QQ or B.base_ring()!=QQ:
        raise ValueError("Algebras must be over rational field")
    ram_A = A.ramified_primes()
    ram_B = B.ramified_primes()
    if ram_A != ram_B:
        return False,[]
    if ram_A ==[]:
        raise ValueError("A,B are split, use the other function")

    alpha, beta = A.invariants()
    a, b = B.invariants()
    i_A, j_A, k_A = A.gens()
    i_B, j_B, k_B = B.gens()

    Falpha = factor(abs(alpha))
    Fbeta = factor(abs(beta))
    Fa = factor(abs(a))
    Fb = factor(abs(b))
   
 
    sol = diagonal_qfsolve([a,b,-a*b,-alpha],factors =[Fa,Fb,Fa*Fb,Falpha])
    mu = sol[0]/sol[3]*i_B+sol[1]/sol[3]*j_B+sol[2]/sol[3]*k_B


    start_time = time.time()
    nu = quaternionic_complement_special(B, mu,LLL=LLL)
    gamma = -nu.reduced_norm() 
    r,s = gamma.numerator(),gamma.denominator()
   
    return max(r.nbits(),s.nbits())







def mean_gamma(t,LLL=True):
    """
    Reads a data file for a given bit length 't', runs the timer on each
    isomorphic pair, and returns the mean computation time.

    Args:
        t (int): The bit length corresponding to the data file.

    Returns:
        float: The mean time for the benchmark, or -1.0 if an error occurs.
    """
    db_directory = "database"
    filename = os.path.join(db_directory, f"{t}_bits_iso.txt")

    if not os.path.exists(filename):
        print(f"Error: Data file not found at '{filename}'")
        return -1.0

    gamma_sizes = []

    print("===========================")
    
    try:
        with open(filename, 'r') as f:
            # Read all lines from the file at once
            lines = f.read().splitlines()

        # The file is in blocks of 4 numbers, sometimes separated by a blank line.
        # We can iterate through the lines in steps of 4, ignoring blank lines.
        # A robust way is to filter out empty lines first.
        
        numeric_lines = [line for line in lines if line.strip()]
        
        if len(numeric_lines) % 4 != 0:
            print(f"Warning: The file {filename} has an incomplete block of numbers.")
        
        # Process in chunks of 4
        for i in range(0, len(numeric_lines), 4):
            # Read the block of four integer coefficients
            a_val = ZZ(numeric_lines[i])
            b_val = ZZ(numeric_lines[i+1])
            c_val = ZZ(numeric_lines[i+2])
            d_val = ZZ(numeric_lines[i+3])

            # Recreate the quaternion algebras. Remember to make the coefficients
            # negative, as per the generation logic for division algebras.
            A = QuaternionAlgebra(QQ, -a_val, -b_val)
            B = QuaternionAlgebra(QQ, -c_val, -d_val)

            # Measure the time for the isomorphism from B to A, as requested
            res = size_of_gamma(B, A,LLL = LLL)
            print(res)
            gamma_sizes.append(res)
            

    except (IOError, IndexError, ValueError) as e:
        print(f"An error occurred while processing {filename}: {e}")
        return -1.0


        
    return sum(gamma_sizes)/len(gamma_sizes)
