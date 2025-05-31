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


def strictly_bigger_order(B, Zbasis_O, primes_disc=None):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of elements of B representing the lattice O := Ze1⊕...⊕ZeN; assumed to be a Z-order.
        -- primes_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of elements of B representing the lattice Gamma = Zf1⊕...⊕ZfN, a strictly bigger order than O if possible; "Maximal" otherwise.
    """
    N = dimension(B)
    disc = discriminant(B, Zbasis_O)
 
    if primes_disc is None:
        factor_disc = factor(disc)
        primes_disc = [p for p, _ in factor_disc]


    for p in primes_disc:
        A, pi_1 = finite_algebra_from_order(B, Zbasis_O, p)
        basis_rad_A = A.radical_basis()
        C, pi_2 = quotient_algebra_ideal(A, basis_rad_A)

        phi = lambda i: pi_2(pi_1(i))

        Zbasis_I_vect = kernel_Z_mod_map(C, phi, N)
        Zbasis_I = [sum(int(Zbasis_I_vect[k][i, 0]) * Zbasis_O[i] for i in range(N)) for k in range(N)]
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
            Zbasis_I = [sum(Zbasis_I_vect[k][i, 0] * Zbasis_O[i] for i in range(N)) for k in range(N)]
            Zbasis_I = lattice_LLL(B,Zbasis_I)
      
            Zbasis_Gamma = left_order(B, Zbasis_I)

            if not are_equals_lattices(B, Zbasis_O, Zbasis_Gamma):
                return Zbasis_Gamma

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

        
        Zbasis_I = [sum(int(Zbasis_I_vect[k][i, 0]) * Zbasis_O[i] for i in range(N)) for k in range(N)]
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
            Zbasis_I = [sum(Zbasis_I_vect[k][i, 0] * Zbasis_O[i] for i in range(N)) for k in range(N)]
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
    go = True
    while go:
        previous = Zbasis_max_O
        Zbasis_max_O = strictly_bigger_order(B,Zbasis_max_O)
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