load("utilities/utilities.sage")
load("maximal_orders/maximal_orders_utilities.sage")
load("maximal_orders/minimal_ideals_from_magma.sage")
load("maximal_orders/minimal_ideals_from_sage.sage")


#------------------------------------------------------------------------------------------
#
#            FIND A MAXIMAL ORDER CONTAINING O IN CENTRAL SIMPLE ALGEBRA
#
#------------------------------------------------------------------------------------------

import time


def strictly_bigger_order_local(B, Zbasis_O, p):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of elements of B representing the lattice O := Ze1⊕...⊕ZeN; assumed to be a Z-order.
        -- primes_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of elements of B representing the lattice Gamma = Zf1⊕...⊕ZfN, a strictly bigger order than O at p if possible; "Maximal loc" otherwise.
    """
    N = dimension(B)

    A, pi_1 = finite_algebra_from_order(B, Zbasis_O, p)
    basis_rad_A = A.radical_basis()
    C, pi_2 = quotient_algebra_ideal(A, basis_rad_A)

    phi = lambda i: pi_2(pi_1(i))

    Zbasis_I_vect = kernel_Z_mod_map(C, phi, N)
    Zbasis_I = [sum(Zbasis_I_vect[k][i] * Zbasis_O[i] for i in range(N)) for k in range(N)]
    Zbasis_I = lattice_LLL(B,Zbasis_I)

    Zbasis_Gamma = left_order(B, Zbasis_I)

    if not are_equals_lattices(B, Zbasis_O, Zbasis_Gamma):
        return Zbasis_Gamma

    # Test minimal non-zero ideals in A/Rad A
        
    MinIdealsList = minimal_ideals_magma(C)

    for basis_K in MinIdealsList:
        D, pi_3 = quotient_algebra_ideal(C, basis_K)

        psi = lambda i: pi_3(pi_2(pi_1(i)))

        Zbasis_I_vect = kernel_Z_mod_map(D, psi, N)
        Zbasis_I = [sum(Zbasis_I_vect[k][i] * Zbasis_O[i] for i in range(N)) for k in range(N)]
        Zbasis_I = lattice_LLL(B,Zbasis_I)
      
        Zbasis_Gamma = left_order(B, Zbasis_I)

        if not are_equals_lattices(B, Zbasis_O, Zbasis_Gamma):
            return Zbasis_Gamma
    return "Maximal loc"


def strictly_bigger_order(B, Zbasis_O, primes_disc=None):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of elements of B representing the lattice O := Ze1⊕...⊕ZeN; assumed to be a Z-order.
        -- primes_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of elements of B representing the lattice Gamma = Zf1⊕...⊕ZfN, a strictly bigger order than O if possible; "Maximal" otherwise.
    """
    
    if primes_disc is None:
        disc = discriminant(B, Zbasis_O)
        factor_disc = factor(disc)
        primes_disc = [p for p, _ in factor_disc]

    for p in primes_disc.copy():
        Zbasis_Gamma = strictly_bigger_order_local(B, Zbasis_O, p)
        if Zbasis_Gamma != "Maximal loc":
            return Zbasis_Gamma
        primes_disc.remove(p)

    return "Maximal"


def strictly_bigger_order_printers(B, Zbasis_O, primes_disc=None):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of elements of B representing the lattice O := Ze1⊕...⊕ZeN; assumed to be a Z-order.
        -- primes_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of elements of B representing the lattice Gamma = Zf1⊕...⊕ZfN, a strictly bigger order than O if possible; "Maximal" otherwise.
    """
    start_total = time.time()
    
    N = dimension(B)
    start_disc = time.time()
    disc = discriminant(B, Zbasis_O)
    time_compute_disc = time.time() - start_disc

    t0=time.time()
    if primes_disc is None:
        factor_disc = factor(disc)
        primes_disc = [p for p, _ in factor_disc]
    time_factor_disc = time.time() -t0
        
    print("discriminant = " + str(factor(disc)))
    print(f"[Timing] Discriminant computed in {time_compute_disc} s")
    print(f"[Timing] Factorization discriminant computed in {time_factor_disc} s")
    print(f"[Info] Prime divisors of discriminant: {primes_disc}")

    for p in primes_disc:
        print(f"[Info] Working at prime p = {p}")
        
        t0 = time.time()
        A, pi_1 = finite_algebra_from_order(B, Zbasis_O, p)
        print(f"[Timing] Finite algebra and projection computed in {time.time() - t0:.4f} s")

        t0 = time.time()
        basis_rad_A = A.radical_basis()
        C, pi_2 = quotient_algebra_ideal(A, basis_rad_A)
        print(f"[Timing] Radical and quotient A/Rad(A) computed in {time.time() - t0:.4f} s")

        t0 = time.time()
        phi = lambda i: pi_2(pi_1(i))
        Zbasis_I_vect = kernel_Z_mod_map(C, phi, N)
        print(f"[Timing] Zero ideal kernel computed in {time.time() - t0:.4f} s")

        Zbasis_I = [sum(Zbasis_I_vect[k][i] * Zbasis_O[i] for i in range(N)) for k in range(N)]
        Zbasis_I = lattice_LLL(B,Zbasis_I)

        t0 = time.time()
        Zbasis_Gamma = left_order(B, Zbasis_I)
        print(f"[Timing] Left order computed in {time.time() - t0:.4f} s")

        if not are_equals_lattices(B, Zbasis_O, Zbasis_Gamma):
            print(f"[Result] Found strictly bigger order at prime {p} (zero ideal case).")
            print(f"[Total Time] {time.time() - start_total:.4f} s")
            print()
            return Zbasis_Gamma

        # Test minimal non-zero ideals in A/Rad A
        t0 = time.time()
        MinIdealsList = minimal_ideals_magma(C)
        print(f"[Timing] Minimal ideals list at p={p} computed in {time.time() - t0:.4f} s")
        print(f"[INFO] Number of minimal ideals : {len(MinIdealsList)}")
        
        for basis_K in MinIdealsList:
            t0 = time.time()
            D, pi_3 = quotient_algebra_ideal(C, basis_K)
            psi = lambda i: pi_3(pi_2(pi_1(i)))
            Zbasis_I_vect = kernel_Z_mod_map(D, psi, N)
            print(f"[Timing] Zero ideal kernel computed in {time.time() - t0:.4f} s")
            Zbasis_I = [sum(Zbasis_I_vect[k][i] * Zbasis_O[i] for i in range(N)) for k in range(N)]
            Zbasis_I = lattice_LLL(B,Zbasis_I)
            t0 = time.time()
            Zbasis_Gamma = left_order(B, Zbasis_I)
            print(f"[Timing] Left order computed in {time.time() - t0:.4f} s")
            if not are_equals_lattices(B, Zbasis_O, Zbasis_Gamma):
                print(f"[Result] Found strictly bigger order at prime {p} (minimal ideal case).")
                print(f"[Total Time] {time.time() - start_total:.4f} s")
                print()
                return Zbasis_Gamma
        print(f"[Timing] All minimal ideals processed at p={p} in {time.time() - t0:.4f} s")

    print(f"[Result] Order is already maximal.")
    print(f"[Total Time] {time.time() - start_total:.4f} s")
    return "Maximal"


def is_maximal_order(B,Zbasis_O,primes_disc = None):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of element of B representing the lattice O := Ze1⊕... ⊕ZeN assume it is an Z-order.
        -- prime_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- boolean -- True if O is a maximal order, False otherwise
    """
    return strictly_bigger_order(B,Zbasis_O,primes_disc) == "Maximal"




def max_order_containing_order(B,Zbasis_O,primes_disc = None):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of element of B representing the lattice O := Ze1⊕... ⊕ZeN assume it is an Z-order.
        -- prime_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of element of B representing the lattice Gamma = Zf1⊕... ⊕ZfN which is a maximal order containing O 
    """
    Zbasis_max_O = Zbasis_O
    if primes_disc is None :
        disc = discriminant(B, Zbasis_O)
        factor_disc = factor(disc)
        primes_disc = [p for p, _ in factor_disc]
    go = True
    while go:
        previous = Zbasis_max_O
        Zbasis_max_O = strictly_bigger_order(B,Zbasis_max_O,primes_disc)
        if Zbasis_max_O == "Maximal":
            go = False
    return previous


def max_order_containing_order_printers(B,Zbasis_O,primes_disc = None):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of element of B representing the lattice O := Ze1⊕... ⊕ZeN assume it is an Z-order.
        -- prime_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of element of B representing the lattice Gamma = Zf1⊕... ⊕ZfN which is a maximal order containing O 
    """
    Zbasis_max_O = Zbasis_O
    go = True
    while go:
        previous = Zbasis_max_O
        Zbasis_max_O = strictly_bigger_order_printers(B,Zbasis_max_O)
        if Zbasis_max_O == "Maximal":
            go = False
    return previous




def max_order(B):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of element of B representing the lattice Gamma = Zf1⊕... ⊕ZfN which is a maximal order of B

    Rq : Start from the lattice I := Ze1⊕... ⊕ZeN where (ei) basis of B, compute the left order O of I and then compute a maximal order containing O.
    """
    Zbasis_I = list(B.basis())
    Zbasis_O = left_order(B,Zbasis_I)
    return max_order_containing_order(B,Zbasis_O)

def max_order_printers(B):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of element of B representing the lattice Gamma = Zf1⊕... ⊕ZfN which is a maximal order of B

    Rq : Start from the lattice I := Ze1⊕... ⊕ZeN where (ei) basis of B, compute the left order O of I and then compute a maximal order containing O.
    """
    Zbasis_I = list(B.basis())
    Zbasis_O = left_order(B,Zbasis_I)
    return max_order_containing_order_printers(B,Zbasis_O)




def max_order_tensor_quat_alg(A,B,Zbasis_O1=None,Zbasis_O2=None):
    """
    INPUT :
        -- A -- quaternion algebra (such that A.basis() = 1,i,j,k)
        -- B -- queternion algebra (such that B.basis() = 1,i,j,k)
        -- Zbasis_O1 -- Z basis of a maximal order in A
        -- Zbasis_O2 -- Z basis of a maximal order in B
    OUTPUT :
        -- Zbasis_Gamma -- a Z basis of a maximal order in C = tensor(A,opposit(B))
    REMARK :
        a1,a2,a3,a4 = 1,i,j,k in A
        b1,b2,b3,b4 = 1,i,j,k in B

        O1 = Ze1 ⊕ Ze2 ⊕ Ze3 ⊕ Ze4
        O2 = Zf1 ⊕ Zf2 ⊕ Zf3 ⊕ Zf4

        O = Z(e1 ⊗ f1) ⊕ Z(e1 ⊗ f2) ⊕ Z(e1 ⊗ f3) ⊕ Z(e1 ⊗ f4) ⊕ 
            Z(e2 ⊗ f1) ⊕ Z(e2 ⊗ f2) ⊕ Z(e2 ⊗ f3) ⊕ Z(e2 ⊗ f4) ⊕
            Z(e3 ⊗ f1) ⊕ Z(e3 ⊗ f2) ⊕ Z(e3 ⊗ f3) ⊕ Z(e3 ⊗ f4) ⊕
            Z(e4 ⊗ f1) ⊕ Z(e4 ⊗ f2) ⊕ Z(e4 ⊗ f3) ⊕ Z(e4 ⊗ f4)

        each ei are express in a1,a2,a3,a4
        each fj are express in b1,b2,b3,b4

        and by definition the basis of C is 
        c1,..,c16 with
        c_{4*k+l} = ak ⊗ bl
        so :

        In C: ei ⊗ fj = (sum_k sk*ak) ⊗ (sum_l rl*bl) = 
                        = sum_k sum_l sk*rl*(ak ⊗ bl) = 
                        = sum_k sum_l sk*rl*c_{4*k+l}

    """
    C = tensor(A,opposite(B))
    BC = C.basis()

    if Zbasis_O1 is None:
        Zbasis_O1 = max_order(A)

    if Zbasis_O2 is None:
        Zbasis_O2 = max_order(B)

    Zbasis_O = []
    for e in Zbasis_O1:
        for f in Zbasis_O2 :
            s = get_coefficients(e,A)
            r = get_coefficients(f,B)
            e_tensor_f = sum( sum (s[k]*r[l]*BC[4*k+l]) for k in range(4) for l in range(4))
            Zbasis_O.append(e_tensor_f)

    Zbasis_Gamma = max_order_containing_order(B,Zbasis_O)

    return Zbasis_Gamma