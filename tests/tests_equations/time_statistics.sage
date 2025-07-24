load("src/isomorphism/explicit_iso_quat_alg_equations.sage")
import numpy


def iso_quat_alg_timer(A, B):
    """
    Executes the isomorphism algorithm and returns a breakdown of CPU times.

    This function is a modified version of the main algorithm, with timing
    calls inserted to measure the performance of distinct computational steps.

    INPUT:
        - ``A``, ``B`` -- Two isomorphic quaternion algebras over QQ.

    OUTPUT:
        - A tuple ``(time_fact_coeff, time_fact_gamma, time_solve)`` containing
          the CPU time in seconds for each part of the algorithm.
    """
    if not A.is_isomorphic(B):
        print(f"A : {A}")
        print(f"B : {B}")
        return None

    total_time_start = cputime()

    # --- Initial Setup ---
    alpha, beta = A.invariants()
    a, b = B.invariants()
    i_B, j_B, k_B = B.gens()

    # --- Timing Block 1: Factorization of initial coefficients ---
    time_fact_coeff_start = cputime()
    Falpha = factor(alpha)
    Fbeta = factor(beta)
    Fa = factor(a)
    Fb = factor(b)
    time_fact_coeff_end = cputime()
    time_fact_coeff = time_fact_coeff_end - time_fact_coeff_start

    # --- Step 2: Find mu in B such that mu^2 = alpha ---
    # This involves solving a norm equation on a conic.
    sol = diagonal_qfsolve([a, b, -a*b, -alpha], factors=[Fa, Fb, Fa * Fb, Falpha])
    mu = sol[0]/sol[3] * i_B + sol[1]/sol[3] * j_B + sol[2]/sol[3] * k_B

    # --- Step 3: Find an element nu that anticommutes with mu ---
    nu = quaternionic_complement(B, mu)
    gamma = -nu.reduced_norm()
    r,s = gamma.numerator(),gamma.denominator()


    #Integer size 
    alpha_size = ZZ(alpha).nbits()
    beta_size = ZZ(beta).nbits()
    a_size = ZZ(a).nbits()
    b_size = ZZ(b).nbits()
    gamma_size = max(r.nbits(),s.nbits())
    print(gamma_size)


    # --- Timing Block 2: Factorization of gamma ---
    time_fact_gamma_start = cputime()
    Fgamma = None  # Initialize Fgamma in case of timeout

    try:
        alarm(10)  # Set a 10-second alarm
        Fr = factor(abs(r))
        Fs = factor(abs(s))
        cancel_alarm()  # If factor() finishes in time, cancel the alarm
    except AlarmInterrupt:
        # This block executes if the alarm goes off
        print("Factorization timed out after 10 seconds. Algorithm information:")
        print(f"alpha ({ZZ(alpha).nbits()} bits): {alpha}")
        print(f"beta ({ZZ(beta).nbits()} bits): {beta}")
        print(f"a ({ZZ(a).nbits()} bits): {a}")
        print(f"b ({ZZ(b).nbits()} bits): {b}")
        print(f"r ({ZZ(r).nbits()} bits): {r}")
        print(f"s ({ZZ(s).nbits()} bits): {s}")
        return [gamma_size]

    time_fact_gamma_end = cputime()
    time_fact_gamma = time_fact_gamma_end - time_fact_gamma_start

    # --- Step 4: Find j' = (x + y*mu)*nu such that (j')^2 = beta ---
    # This involves solving a second norm equation.
    sol2 = diagonal_qfsolve([r,-alpha*r,-beta*s],factors = [Fr,Falpha*Fr,Fbeta*Fs])


    # --- Final Calculation ---
    total_time_end = cputime()
    total_time = total_time_end - total_time_start

    # time_solve is the total time minus the time spent on factorization.
    time_solve = total_time - time_fact_coeff - time_fact_gamma

    
    #return time_fact_coeff, time_fact_gamma, time_solve,alpha_size,beta_size,a_size,b_size,gamma_size
    return total_time,gamma_size




# ==============================================================================
#  FUNCTION TO READ A DATA FILE AND COMPUTE MEAN TIME
# ==============================================================================

def mean_time(t):
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

    timings = []
    gamma_sizes = []
    
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
            res = iso_quat_alg_timer(B, A)
            if len(res)>1:
                time_taken,gamma_size = res
                timings.append(time_taken)
                gamma_sizes.append(gamma_size)
            else:
                gamma_sizes.append(res[0])

    except (IOError, IndexError, ValueError) as e:
        print(f"An error occurred while processing {filename}: {e}")
        return -1.0


        
    return sum(timings) / len(timings), sum(gamma_sizes)/len(gamma_sizes),100-len(timings)




dict_time_July_16 = {
    10:(0.00772742000000008, 36.58, 0),
    20:(0.010123919999999877, 83.78, 0),
    30:(0.01599833000000004, 134.61, 0),
    40:(0.05296614999999992, 181.01, 0),
    50:(0.3633695300000005, 232.79, 0),
    60:(1.1000505000000027, 280.3854166666667, 4),
    70:(1.8428254835164877, 307.6043956043956, 9),
    75:(2.1186574805194818, 320.85714285714283, 23),
    80:(2.888688142857128, 342.5079365079365, 37),
    85:(3.969353145454546, 373.3272727272727, 45),
}


dict_time_July_19 = {
    10: (0.008852219999999952, 33.33, 0),
    20: (0.011558359999999998, 71.78, 0),
    30: (0.01602070999999994, 110.8, 0),
    40: (0.03844669000000005, 149.42, 0),
    50: (0.18481579000000004, 189.22, 0),
    60: (0.8118399690721647, 232.03, 3),
    70: (1.6559694239130436, 262.09, 8),
    75: (2.133714182926835, 283.63, 18),
    80: (2.558688893939387, 306.23, 34),
    85: (3.840828423728835, 337.8, 41),
    90: (3.789440850000006, 351.97, 40),
    95: (6.614299259999988, 358.52, 50),
    100: (4.924485827586131, 387.38, 71)
}