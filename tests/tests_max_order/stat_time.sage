load("src/isomorphism/explicit_iso_quat_alg_equations.sage")
load("src/maximal_orders/find_maximal_orders.sage")


A = QuaternionAlgebra(QQ,-1,-7)
i,j,k = A.gens()


import time
def mean_time(t):
    timing = []
    while len(timing)<10:
        M = int(2^(t/2))
        mu = randint(-M,M)*i+ randint(-M,M)*j+ randint(-M,M)*k
        nu = quaternionic_complement(A,mu)
        a = -mu.reduced_norm()
        b = -nu.reduced_norm()
        B = QuaternionAlgebra(QQ,a,b)

        mu0 = randint(-M,M)*i+ randint(-M,M)*j+ randint(-M,M)*k
        nu0 = quaternionic_complement(A,mu0)
        c = -mu0.reduced_norm()
        d = -nu0.reduced_norm()
        C = QuaternionAlgebra(QQ,c,d)

        print(ZZ(a).nbits(),ZZ(b).nbits(),ZZ(c).nbits(),ZZ(d).nbits())

        # --- Pre-computation with a 10-second timeout ---
        try:
            alarm(10)  # Set the alarm
            BOB = B.maximal_order().basis()
            BOC = C.maximal_order().basis()
            cancel_alarm()  # Cancel the alarm if successful (Correct function name)
        except AlarmInterrupt:
            # If the alarm goes off, print a message and try with new algebras
            print("    -> Pre-computation timed out, trying new algebras...")
            continue  # Skip to the next iteration of the while loop

    
        t0 =time.time()
        BL = max_order_tensor_quat_alg(B,C,Zbasis_O1=BOB,Zbasis_O2=BOC)
        t1 = time.time()

        timing.append(t1-t0)
        print(f"    -> Sample {len(timing)}/10 collected.")

    return sum(timing)/len(timing)





dico = {
    10:1.9403368473052978,
    20:2.487557530403137,
    30:2.6674673318862916,
    40:2.866751527786255,
    50:3.155496525764465,
    60:3.4114906072616575,
    70:3.144124484062195,
    80:3.0105661869049074,
    90:3.5098113775253297,
    100:2.9791826725006105,
    110:3.334433102607727,
    120:4.4708130121231076,
    125:4.347476243972778,
    130:3.4875132560729982,
    135:3.3489542245864867,
    140:3.524424338340759,
    145:3.217971420288086,
    150:5.019127273559571
}
