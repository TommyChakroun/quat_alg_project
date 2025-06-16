load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("maximal_orders/find_maximal_orders.sage")
load("explicit_embeding_real_number/explicit_embedding.sage")
load("maximal_orders/maximal_orders_utilities.sage")


#------------------------------------------------------------------------------------------
#
#            RANK ONE ELEMENT IN ALGEBRA ISOMORPHIC TO Mn(Q)
#
#------------------------------------------------------------------------------------------

##---List of functions---##

# right_rank(A,x)

def right_rank(A, x):
    """
    INPUT :
        -- A -- an algebra
        -- x -- an element of A
    OUTPUT :
        -- r -- the rank of A ->A,a ->ax i.e the dimension of Ax over F
    """
    F = A.base_ring()
    basis_A = list(A.basis())

    rows = []

    for e in basis_A:
        y = e * x
        coords = coordinate(y, A, basis_A)
        rows.append(coords)
    
    M = Matrix(F, rows)
    r = M.rank()
    return r



def rank_one(A, Zbasis_O_max=None):
    """
    INPUT :
        -- A -- an algebra over Q given by structure constants, assumed to be isomorphic to Mn(Q).
        -- Zbasis_O_max -- a Z-basis of a maximal order in A, if provided.
    OUTPUT :
        -- a -- an element of A of "rank one," meaning the dimension of Aa over Q is n.
               Note: The original function's return statement has been modified to reflect this docstring.
    """

    print("===================================================================")
    print("                      🚀 Starting rank_one function 🚀             ")
    print("===================================================================")
    print("\n")

    print("--------------- Input Parameters ---------------")
    print(f"Algebra A: {A}")
    print(f"Z-basis of maximal order (Zbasis_O_max): {Zbasis_O_max}")
    print("\n" * 2) # Extra space for readability

    dim_A = A.dimension()
    print(f"✨ Dimension of algebra A: {dim_A}")

    n = int(sqrt(dim_A))
    print(f"✨ Calculated n (sqrt of dim_A): {n}")
    print("\n" * 2)

    if Zbasis_O_max is None:
        print("--------------- Maximal Order Calculation ---------------")
        print("💡 Zbasis_O_max not provided, computing maximal order...")
        Zbasis_O = A.maximal_order()
        print(f"Computed maximal order Zbasis_O (first 5 elements):")
        for i, el in enumerate(Zbasis_O[:min(5, len(Zbasis_O))]):
            print(f"    - {el}")
    else:
        print("--------------- Using Provided Maximal Order ---------------")
        print("✅ Using provided Zbasis_O_max.")
        Zbasis_O = Zbasis_O_max
    
    # Ensure Zbasis_O is a list of elements, not a basis object itself
    if hasattr(Zbasis_O, 'list'):
        Zbasis_O = Zbasis_O.list()
    print("\n" * 2)

    print("--------------- Real Embedding Process ---------------")
    E, embedding_dict = explicit_real_embedding(A)
    print(f"Real embedding space E: {E}")
    # print(f"Embedding dictionary: {embedding_dict}") # This can be very large, uncomment with caution.
    print("\n" * 2)

    print("--------------- Pushing Lattice Embedding ---------------")
    Zbasis_O_in_EN = push_lattice_embedding(A, Zbasis_O, E, embedding_dict)
    print(f"Z-basis of maximal order in embedded space (first 5 elements):")
    for i, el in enumerate(Zbasis_O_in_EN[:min(5, len(Zbasis_O_in_EN))]):
        print(f"    - {el}")
    print("\n" * 2)

    print("--------------- Approximating Lattice ---------------")
    Zbasis_O_in_EN_approx = approx_lattice(E, Zbasis_O_in_EN, precision=1e-5)
    print(f"Approximated lattice basis (first 5 elements):")
    for i, el in enumerate(Zbasis_O_in_EN_approx[:min(5, len(Zbasis_O_in_EN_approx))]):
        print(f"    - {el}")
    print("\n" * 2)
 
    print("--------------- Constructing Matrix for LLL ---------------")
    # Convert Zbasis_O_in_EN_approx to a matrix in Sage.
    # It's usually a list of vectors, so we construct a matrix from its rows.
    A_matrix_for_LLL = Matrix(QQ, Zbasis_O_in_EN_approx)
    print(f"Matrix for LLL (first 5 rows):")
    for i, row in enumerate(A_matrix_for_LLL.rows()[:min(5, A_matrix_for_LLL.nrows())]):
        print(f"    - {row}")
    print("\n" * 2)

    print("--------------- Performing LLL Reduction ---------------")
    N, T = A_matrix_for_LLL.LLL(transformation = True)
    print("LLL reduced basis N (first 5 rows):")
    for i, row in enumerate(N.rows()[:min(5, N.nrows())]):
        print(f"    - {row}")
    print("\n")

    print("Transformation matrix T:")
    print(T)
    print(f"Determinant of T: {T.determinant()}")
    print("\n" * 2)

    print("--------------- Applying Transformation ---------------")
    Zbasis_O_new = [ sum(T[i][k]*Zbasis_O[k] for k in range(dim_A)) for i in range(dim_A) ]
    print("New Z-basis of maximal order (first 5 elements):")
    for i, el in enumerate(Zbasis_O_new[:min(5, len(Zbasis_O_new))]):
        print(f"    - {el}")
    print("\n" * 2)

    print("===================================================================")
    print("                          ✅ Results ✅                            ")
    print("===================================================================")
    print("\n")

    print("--------------- Lattice Equality Check ---------------")
    try:
        # Using your personal function for lattice equality
        if are_equals_lattices(A, Zbasis_O, Zbasis_O_new):
            print("👍 The original Z-basis and the new Z-basis span the same lattice (Z-module).")
        else:
            print("⚠️ The original Z-basis and the new Z-basis DO NOT span the same lattice.")
    except NameError:
        print("❌ Error: 'are_equals_lattices' function not found. Please define it.")
    except Exception as e:
        print(f"❌ An error occurred while checking lattice equality: {e}")
    print("\n" * 2)

    print("--------------- Discriminant of Zbasis_O_new ---------------")
    try:
        # Using your personal function for discriminant
        calculated_discriminant = discriminant(A, Zbasis_O_new)
        print(f"🔢 Discriminant of Zbasis_O_new: {calculated_discriminant}")
        
        if calculated_discriminant == 1:
            print("🎉 The discriminant of Zbasis_O_new is 1, as expected for a maximal order in Mn(Q).")
        else:
            print(f"⚠️ Warning: The discriminant of Zbasis_O_new is {calculated_discriminant}, not 1.")
            print("This might indicate that Zbasis_O_new is not a basis for a maximal order with discriminant 1.")
    except NameError:
        print("❌ Error: 'discriminant' function not found. Please define it.")
    except Exception as e:
        print(f"❌ An error occurred while calculating the discriminant: {e}")
    print("\n" * 2)

    for ite in range(1000):
        alpha = [randint(1,100) for k in range(dim_A)]
        x = sum(alpha[i]*Zbasis_O_new[i] for i in range(dim_A))
        
        r = right_rank(A, x)
        
        if r<dim_A:
            print(f"\n🎉 Found an element of rank strict less than 16 after {ite+1} iterations: x = {x}")
            print("===================================================================")
            print("                    🏁 rank_one function finished 🏁               ")
            print("===================================================================")
            print("\n")
            return x


    print("\n⚠️ Could not find a rank one element after 100 iterations.")
    print("===================================================================")
    print("                    🏁 rank_one function finished 🏁               ")
    print("===================================================================")
    print("\n")

    return None


