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



