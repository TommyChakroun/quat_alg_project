load("utilities/utilities.sage")
load("src/iso_splitting_algebra/maximal_orders/maximal_orders_utilities.sage")
load("src/iso_splitting_algebra/minimal_ideals/minimal_ideals_manually.sage")


#------------------------------------------------------------------------------------------
#
#            FIND A MAXIMAL ORDER CONTAINING O IN CENTRAL SIMPLE ALGEBRA
#
#------------------------------------------------------------------------------------------

import time



### --------- Main algo : strictly bigger order at a prime p  ---------- ###


def strictly_bigger_order_local(B, Zbasis_O, p,lattice_format = "LLL"):
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

    if lattice_format == "LLL":
        Zbasis_I = lattice_LLL(B,Zbasis_I)
    else :
        Zbasis_I = lattice_hermite(B,Zbasis_I)

    Zbasis_Gamma = left_order(B, Zbasis_I,lattice_format)

    if not are_equals_lattices(B, Zbasis_O, Zbasis_Gamma):
        return Zbasis_Gamma

    # Test minimal non-zero ideals in A/Rad A
        
    MinIdealsList = minimal_ideals_perso(C)

    for basis_K in MinIdealsList:
        D, pi_3 = quotient_algebra_ideal(C, basis_K)

        psi = lambda i: pi_3(pi_2(pi_1(i)))

        Zbasis_I_vect = kernel_Z_mod_map(D, psi, N)
        Zbasis_I = [sum(Zbasis_I_vect[k][i] * Zbasis_O[i] for i in range(N)) for k in range(N)]
        Zbasis_I = lattice_LLL(B,Zbasis_I)
      
        Zbasis_Gamma = left_order(B, Zbasis_I,lattice_format)

        if not are_equals_lattices(B, Zbasis_O, Zbasis_Gamma):
            return Zbasis_Gamma
    return "Maximal loc"

def strictly_bigger_order_local_printers(B, Zbasis_O, p,lattice_format = "LLL"):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of elements of B representing the lattice O := Ze1⊕...⊕ZeN; assumed to be a Z-order.
        -- primes_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of elements of B representing the lattice Gamma = Zf1⊕...⊕ZfN, a strictly bigger order than O at p if possible; "Maximal loc" otherwise.
    """
    start_total = time.time()

    print(f"[Info] Working at prime p = {p}")

    N = dimension(B)
        
    t0 = time.time()
    A, pi_1 = finite_algebra_from_order(B, Zbasis_O, p)
    print(f"[Timing] Finite algebra and projection computed in {time.time() - t0:.4f} s")

    t0 = time.time()
    basis_rad_A = A.radical_basis()
    C, pi_2 = quotient_algebra_ideal(A, basis_rad_A)
    phi = lambda i: pi_2(pi_1(i))
    print(f"[Timing] Radical and quotient A/Rad(A) computed in {time.time() - t0:.4f} s")

    
    t0 = time.time()
    Zbasis_I_vect = kernel_Z_mod_map(C, phi, N)
    Zbasis_I = [sum(Zbasis_I_vect[k][i] * Zbasis_O[i] for i in range(N)) for k in range(N)]
    if lattice_format == "LLL":
        Zbasis_I = lattice_LLL(B,Zbasis_I)
    else :
        Zbasis_I = lattice_hermite(B,Zbasis_I)
    print(f"[Timing] Zero ideal kernel computed in {time.time() - t0:.4f} s")

    t0 = time.time()
    Zbasis_Gamma = left_order(B, Zbasis_I,lattice_format)
    print(f"[Timing] Left order computed in {time.time() - t0:.4f} s")

    if not are_equals_lattices(B, Zbasis_O, Zbasis_Gamma):
        print(f"[Result] Found strictly bigger order at prime {p} (zero ideal case).")
        print(f"[Total Time] {time.time() - start_total:.4f} s")
        print()
        return Zbasis_Gamma

    # Test minimal non-zero ideals in A/Rad A
    t0 = time.time()
    MinIdealsList = minimal_ideals_perso(C)
    print(f"[Timing] Minimal ideals list at p={p} computed in {time.time() - t0:.4f} s")
    print(f"[INFO] Number of minimal ideals : {len(MinIdealsList)}")
        
    for basis_K in MinIdealsList:
        t0 = time.time()
        D, pi_3 = quotient_algebra_ideal(C, basis_K)
        psi = lambda i: pi_3(pi_2(pi_1(i)))
        print(f"[Timing] Quotient by ideal computed in {time.time() - t0:.4f} s")


        t0 = time.time()
        Zbasis_I_vect = kernel_Z_mod_map(D, psi, N)
        Zbasis_I = [sum(Zbasis_I_vect[k][i] * Zbasis_O[i] for i in range(N)) for k in range(N)]
        Zbasis_I = lattice_LLL(B,Zbasis_I)
        print(f"[Timing] Zero ideal kernel computed in {time.time() - t0:.4f} s")

        t0 = time.time()
        Zbasis_Gamma = left_order(B, Zbasis_I,lattice_format)
        print(f"[Timing] Left order computed in {time.time() - t0:.4f} s")

        if not are_equals_lattices(B, Zbasis_O, Zbasis_Gamma):
            print(f"[Result] Found strictly bigger order at prime {p} (minimal ideal case).")
            print(f"[Total Time] {time.time() - start_total:.4f} s")
            print()
            return Zbasis_Gamma

    print(f"[Result] Order already maximal at prime {p}.")
    print(f"[Total Time] {time.time() - start_total:.4f} s")
    return "Maximal loc"



### --------- The followings algortihme are easily deduced  ---------- ###




## Strictly bigger order global


def strictly_bigger_order(B, Zbasis_O, primes_disc=None, printers = False,lattice_format = "LLL"):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of elements of B representing the lattice O := Ze1⊕...⊕ZeN; assumed to be a Z-order.
        -- primes_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of elements of B representing the lattice Gamma = Zf1⊕...⊕ZfN, a strictly bigger order than O if possible; "Maximal" otherwise.
    """
    
    t0 = time.time()
    disc = discriminant(B, Zbasis_O)
    t1 = time.time()
    factor_disc = MyFactor(disc,primes=primes_disc)
    primes_disc = [p for p, _ in factor_disc]
    t2 = time.time()

    if printers :
        print(f"[INFO] discriminant :{disc}={factor_disc}")
        print(f"[INFO] list of primes divisor of the discriminant : {primes_disc}")
        print(f"[TIMING] Compute discriminant : {t1 - t0:.4f} s")
        print(f"[Timing] Factor the discriminant : {t2 - t1:.4f}s")

    for p in primes_disc:
        if printers:
            Zbasis_Gamma = strictly_bigger_order_local_printers(B, Zbasis_O, p,lattice_format)
        else :
            Zbasis_Gamma = strictly_bigger_order_local(B, Zbasis_O, p,lattice_format)

        if Zbasis_Gamma != "Maximal loc":
            return Zbasis_Gamma

        primes_disc.remove(p)

    return "Maximal"



## Is maximal Order global



def is_maximal_order(B,Zbasis_O,primes_disc = None,lattice_format = "LLL"):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of element of B representing the lattice O := Ze1⊕... ⊕ZeN assume it is an Z-order.
        -- prime_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- boolean -- True if O is a maximal order, False otherwise
    """
    return strictly_bigger_order(B,Zbasis_O,primes_disc,lattice_format) == "Maximal"





## Maximal Order at a prime p containing a given order 


def max_order_containing_order_local(B,Zbasis_O,p,printers = False,lattice_format = "LLL"):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of element of B representing the lattice O := Ze1⊕... ⊕ZeN assume it is an Z-order.
        -- p -- a prime 
        -- prime_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of element of B representing the lattice Gamma = Zf1⊕... ⊕ZfN which is a maximal order at the prime p containing O 
    """
    Zbasis_max_O = Zbasis_O
    go = True
    if printers:
        print()
        print(f" --- Compute an overorder which maximal at the prime {p}... --- ")
        print()
    while go:
        previous = Zbasis_max_O
        if printers:
            Zbasis_max_O = strictly_bigger_order_local_printers(B,Zbasis_max_O,p,lattice_format)
        else :
            Zbasis_max_O = strictly_bigger_order_local(B,Zbasis_max_O,p,lattice_format)
        if Zbasis_max_O == "Maximal loc":
            go = False
    if printers:
        print()
        print(f"Overorder which maximal at the prime {p} succesfully compute.")
    return previous







## Maximal Order global containing a given order : itterative and parallel method


def max_order_containing_order(B,Zbasis_O,primes_disc = None,parallel = False, printers = False,lattice_format = "LLL"):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
        -- Zbasis_O -- a list e1,..,eN of element of B representing the lattice O := Ze1⊕... ⊕ZeN assume it is an Z-order.
        -- prime_disc -- the list of prime factors of disc(O)
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of element of B representing the lattice Gamma = Zf1⊕... ⊕ZfN which is a maximal order containing O 
    """

    
    disc = discriminant(B, Zbasis_O)
    if printers:
        print(" --- Start computing a maximal over order --- ")
        print()
        print(f"starting discriminant  : {disc}")
        print()
        print("start factoring the discriminant...")

    t0 = time.time()
    factor_disc = MyFactor(disc,primes=primes_disc)
    primes_disc = [p for p, _ in factor_disc]
    t1=time.time()

    if printers:
        print(f"done in {t1 - t0:.4f} s")
        print(f"disc =  {factor_disc}")

    


    if parallel:
        if printers:
            print("We use the parallel algorithm.")
            print()
        list_Zbasis_max_loc = []
        for p in primes_disc :
            list_Zbasis_max_loc.append(max_order_containing_order_local(B,Zbasis_O,p,printers,lattice_format))

        if printers :
            print()
            print(" --- Compute the sum of each locally maximal orders... ---")
        t0 = time.time()
        rows = []
        for Zbasis_L in list_Zbasis_max_loc:
            for e in Zbasis_L:
                rows.append(get_coefficients(e,B))
        M = Matrix(QQ,rows)

        vectors = Z_mod_span_by_rows(M)
        basis_B = B.basis()
        N = dimension(B)
        Zbasis_L = [ sum (x[i]*basis_B[i] for i in range(N)) for x in vectors]

        if printers :
            print(f"Sum of each locally maximal orders successfully compute in {time.time() - t0:.4f} s.")
            print()
        
        return Zbasis_L
    else :
        if printers:
            print("We use the iterrative algorithm.")
            print()
        Zbasis_max_O = Zbasis_O
        go = True
        while go:
            previous = Zbasis_max_O
            Zbasis_max_O = strictly_bigger_order(B,Zbasis_max_O,primes_disc,printers,lattice_format)
            if Zbasis_max_O == "Maximal":
                go = False
        return previous







## Maximal Order in an algebra : start from the left order of the canonical lattice of  B


def max_order(B,parallel = False, printers = False,lattice_format = "LLL"):
    """
    INPUT :
        -- B -- a central simple algebra over Q of dimension N
    OUTPUT :
        -- Zbasis_Gamma -- a list f1,..,fN of element of B representing the lattice Gamma = Zf1⊕... ⊕ZfN which is a maximal order of B

    Rq : Start from the lattice I := Ze1⊕... ⊕ZeN where (ei) basis of B, compute the left order O of I and then compute a maximal order containing O.
    """
    Zbasis_I = list(B.basis())
    Zbasis_O = left_order(B,Zbasis_I,lattice_format)
    return max_order_containing_order(B,Zbasis_O,parallel=parallel,printers=printers,lattice_format=lattice_format)








## Maximal Order in A⊗Bop where A and B are quaternion algebra


def max_order_tensor_quat_alg(A,B,C = None,Zbasis_O1=None,Zbasis_O2=None,printers=False):
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
    if C == None:
        C = tensor(A,opposite(B))
    BC = C.basis()

    if Zbasis_O1 is None:
        Zbasis_O1 = list(A.maximal_order(A).basis())  # for quaternion algebra the inner sage function maximal_order is way faster than our general function

    if Zbasis_O2 is None:
        Zbasis_O2 = list(B.maximal_order(B).basis())

    Zbasis_O = []
    for e in Zbasis_O1:
        s = get_coefficients(e, A)
        for f in Zbasis_O2:
            r = get_coefficients(f, B)
            e_tensor_f = sum(sum(s[k]*r[l]*BC[4*k+l] for k in range(4)) for l in range(4))
            Zbasis_O.append(e_tensor_f)


    primes_disc = A.ramified_primes()+B.ramified_primes()
    Zbasis_Gamma = max_order_containing_order(C,Zbasis_O,primes_disc = primes_disc,printers=printers)

    return Zbasis_Gamma



