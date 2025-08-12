load("src/iso_solving_quadratics/solving_quadratics.sage")



# ==============================================================================
#  FUNCTION TO GENERATE A SINGLE ISOMORPHIC PAIR
# ==============================================================================

def random_isomorphic_quat_alg(t):
    """
    Generates coefficients for two isomorphic quaternion algebras.
    B = (-a, -b | Q) and A = (-c, -d | Q)
    The bit lengths will be approximately (t, t, 2t, 3t).

    Args:
        t (int): The target bit length for the first two coefficients.

    Returns:
        A tuple (a, b, c, d) of positive integers.
    """
    # 1. Create the first algebra B = (a', b' | Q)
    # We use negative coefficients to ensure it's a division algebra.
    a_prime = -ZZ.random_element(2**(t-1), 2**t)
    b_prime = -ZZ.random_element(2**(t-1), 2**t)
    B = QuaternionAlgebra(QQ, a_prime, b_prime)
    i, j, k = B.gens()

    # 2. Create a "random" pure quaternion mu
    # Ensure mu is not zero.
    while True:
        coeffs = [randint(-5, 5) for _ in range(3)]
        if any(c != 0 for c in coeffs):
            mu = coeffs[0]*i + coeffs[1]*j + coeffs[2]*k
            break
            
    # 3. Find a small anti-commuting partner nu
    nu = quaternionic_complement(B, mu)

    # 4. The new coefficients are the reduced norms.
    # We take the negative to match the convention A = (-c, -d | Q).
    # The .numerator() call ensures we have an integer.
    c_prime = -mu.reduced_norm().numerator()
    d_prime = -nu.reduced_norm().numerator()
    
    # Return the absolute values for the database file
    return -a_prime, -b_prime, -c_prime, -d_prime

# ==============================================================================
#  MAIN FUNCTION TO GENERATE DATABASE FILES
# ==============================================================================

def generate_files():
    """
    Creates a directory named 'database' and populates it with files
    containing coefficients of isomorphic quaternion algebras.
    """
    # Define the bit lengths for the files and the number of entries per file
    bit_lengths = [10, 20, 30, 40, 50, 60, 70, 75, 80, 85, 90, 95, 100]
    num_entries_per_file = 100
    db_directory = "database"


    print("Starting database generation...")

    # Loop through each specified bit length
    for t in bit_lengths:
        filename = os.path.join(db_directory, f"{t}_bits_iso.txt")
        print(f"  -> Generating file: {filename} ...")

        try:
            with open(filename, 'w') as f:
                # Generate the specified number of entries for the file
                for i in range(num_entries_per_file):
                    # Fetch the four coefficients
                    a, b, c, d = random_isomorphic_quat_alg(t)
                    
                    # Write the block of 4 lines
                    f.write(f"{a}\n")
                    f.write(f"{b}\n")
                    f.write(f"{c}\n")
                    f.write(f"{d}\n")
                    
                    # Add a blank line between blocks for readability
                    if i < num_entries_per_file - 1:
                        f.write("\n")
        except Exception as e:
            print(f"An error occurred while writing to {filename}: {e}")
            
    print("Database generation complete!")

