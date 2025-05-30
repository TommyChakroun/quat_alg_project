load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("maximal_orders/find_maximal_orders.sage")

#------------------------------------------------------------------------------------------
#
#            TEST FIND MAXIMAL ORDER
#
#------------------------------------------------------------------------------------------



def test_7():
    A = QuaternionAlgebra(QQ,-1,-1)
    one_A,i,j,k = A.basis()
    for a in [1/2,1/3,1/5]:
        for b in [1/2,1/3,1/5]:
            Zbasis_I = [a*one_A,b*i,a*j,k]
            I = A.ideal(Zbasis_I)
            Zbasis_L = I.left_order().basis()
            Zbasis_O = left_order(A,Zbasis_I)
            Zbasis_O_reduced = lattice_LLL(B,Zbasis_O)
            print("Zbasis_O = " + str(Zbasis_O))
            print("Zbasis_O_reduced = " + str(Zbasis_O_reduced))
            print("Zbasis_L = " + str(Zbasis_L))
            print(" O = L = Ored : "+ str(are_equals_lattices(A,Zbasis_O,Zbasis_L) and are_equals_lattices(A,Zbasis_O,Zbasis_O_reduced)))


def test_8():
    for a in [-1,1,2]:
        for b in [-1,1,-3]:
            print("------------------------------")
            print(" ")
            print(f"Quaternion Algebra ({a},{b}|Q)")
            print("We start with the ideal I of basis:")
            
            B = QuaternionAlgebra(QQ, a, b)
            Zbasis_I = list(B.basis())  # = [1, i, j, k]
            print(Zbasis_I)

            print("Now we start a sequence of bases of strictly increasing order, the first one is the left order of I.")
            
            i = 0
            Zbasis_O = left_order(B, Zbasis_I)
            O = Zbasis_O  # keep a placeholder to compute discriminant

            while not is_maximal_order(B, Zbasis_O):
                print(str(i) + " : " + str(Zbasis_O) + " with disc = " + str(discriminant(B,Zbasis_O)))
                Zbasis_O = strictly_bigger_order(B, Zbasis_O)
                i += 1
            print(str(i) + " : " + str(Zbasis_O) + " with disc = " + str(discriminant(B,Zbasis_O)))
            print()
            print("Finally we obtain the maximal order with basis:")
            print(Zbasis_O)
            print("disc = " + str(discriminant(B,Zbasis_O)))
            print(" ")
            print("----------------------------------")


def test_9():
    print("==========================================")
    print("   TEST OF MAXIMAL ORDERS IN MATRIX SPACES")
    print("==========================================\n")

    for n in range(2, 6):
        print("------------------------------------------")
        print(f"Matrix algebra M_{n}(Q)")
        print("------------------------------------------\n")

        A = MatrixSpace(QQ, n, n)
        print(f"Basis of M_{n}(Q):")
        Zbasis_O = list(A.basis())
        print(Zbasis_O)
        print()

        print("Checking if the standard basis gives a maximal order...")
        result = is_maximal_order(A, Zbasis_O)
        print("Result:", result)
        print()


def mixed_matrix_space(F,n):
    A = MatrixSpace(F,n,n)
    e = list(A.basis())
    P = random_invertible_matrix(F, n)
    P_inv = P.inverse()
    random_basis_A = [P*e[k]*P_inv for k in range(n**2)]
    B = finite_dimensional_algebra_format(A, random_basis_A)
    return B,random_basis_A


def test_10(n):
    print("==========================================")
    print("   TEST FOR MIXED MATRIX SPACES")
    print("==========================================\n")
    print("------------------------------------------")
    print(f"Matrix algebra M_{n}(Q)")
    print("------------------------------------------\n")

    B,random_basis = mixed_matrix_space(QQ,n)
    print(f"Basis of M_{n}(Q):")
    Zbasis_I = list(B.basis())
    print(random_basis)
    print()
    print("Compute left order  : ")
    Zbasis_O = left_order(B,Zbasis_I)
    print("Done. ")
    print("----------------------------------")
    print("Discriminant of the left order O :")
    print(discriminant(B,Zbasis_O))
    print("----------------------------------")
    print()
    print("==========================================")
    print(" START COMPUTING A MAXIMAL ORDER CONTAINING O ")
    print("==========================================\n")
    print()
    Zbasis_L = max_order_containing_order(B,Zbasis_O)
    print("Done.")
    print()
    print(" Check dicriminant L is +-1 :")
    print("discriminant = " + str(factor(discriminant(B,Zbasis_L))))
    print()


def test_11():
    random_basis = [
    Matrix(QQ, [
        [2, -2, -1, 0],
        [1, 1, -1/2, -1],
        [0, 1/2, 1/2, 0],
        [2, 1/2, 1, 0]
    ]),
    Matrix(QQ, [
        [1, 0, 2, 1/2],
        [0, 0, 1, -1],
        [0, 2, 0, 0],
        [0, -1, 2, 1]
    ]),
    Matrix(QQ, [
        [-1, -1, 1/2, 1],
        [0, 0, -1, -1],
        [-1, -1, 0, -1],
        [0, 1, 2, -1]
    ]),
    Matrix(QQ, [
        [1, 2, 0, 1],
        [0, 0, -1, 0],
        [1, 2, 0, -2],
        [0, -1, 0, 1]
    ]),
    Matrix(QQ, [
        [0, 1, 0, 0],
        [-1, 0, 0, 0],
        [-2, -1/2, -1, -1],
        [0, -1/2, -1, -2]
    ]),
    Matrix(QQ, [
        [-1/2, 0, 1, 2],
        [0, 1/2, -1/2, 2],
        [0, 0, 1, 0],
        [0, 0, 1/2, -2]
    ]),
    Matrix(QQ, [
        [1/2, 0, 2, 1],
        [0, 2, 0, 2],
        [0, 1, 0, 0],
        [2, -1, 2, 0]
    ]),
    Matrix(QQ, [
        [0, -2, 2, 1/2],
        [1/2, 1, 1, -1],
        [-1/2, 2, 0, 0],
        [0, 0, 1, -1]
    ]),
    Matrix(QQ, [
        [1, 0, -1, 2],
        [2, -1/2, 1/2, -1],
        [-1/2, -1, 0, 1/2],
        [-1/2, 0, 1, 0]
    ]),
    Matrix(QQ, [
        [0, 0, 2, -2],
        [0, 0, 1, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [2, 0, -1/2, -2],
        [0, 1, 1, 2],
        [-1, 0, -1/2, -2],
        [-2, 1, 0, -1]
    ]),
    Matrix(QQ, [
        [0, -2, 1/2, 0],
        [0, 1, 0, -1/2],
        [-2, -2, 0, 0],
        [-2, -1, 0, 1]
    ]),
    Matrix(QQ, [
        [1, 0, 0, -1],
        [-1, 0, 1, 0],
        [-1, -1, 1, 1],
        [0, 1, -1, 0]
    ]),
    Matrix(QQ, [
        [-2, 1, -1, 1],
        [0, 2, -1, -1],
        [0, -2, 0, 0],
        [1/2, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [-1/2, -1, 1, 1],
        [1/2, 0, 1, -1],
        [0, 0, 0, 1],
        [1, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [0, -2, -1/2, 0],
        [0, 0, -2, 1/2],
        [-1, 1/2, 1, 2],
        [-2, 0, -1, 0]
    ])  
    ]

    A = MatrixSpace(QQ,4,4)
    B = finite_dimensional_algebra_format(A, random_basis)

    print("==========================================  ")
    print("   STANDARD TEST M4(Q) WITH BASIS           ")
    print("==========================================\n")

    print(random_basis)
    print()
    Zbasis_I = list(B.basis())
    print("Compute left order  : ")
    Zbasis_O = left_order(B,Zbasis_I)
    print("Done. ")
    print("----------------------------------")
    print("Discriminant of the left order O :")
    print(discriminant(B,Zbasis_O))
    print("----------------------------------")
    print()
    print("==========================================")
    print(" START COMPUTING A MAXIMAL ORDER CONTAINING O ")
    print("==========================================\n")
    print()
    Zbasis_L = max_order_containing_order(B,Zbasis_O)
    print("Done.")
    print()
    print(" Check dicriminant L is +-1 :")
    print("discriminant = " + str(factor(discriminant(B,Zbasis_L))))
    print()



def test_12():
    Zbasis_O = [
    Matrix(QQ, [
        [-1, 0, 0, 0],
        [ 0, -1, 0, 0],
        [ 0,  0, -1, 0],
        [ 0,  0,  0, -1]
    ]),
    Matrix(QQ, [
        [0, 0, 0, 0],
        [-3075398820, 0, 0, 0],
        [0, 0, 0, 0],
        [-3075398820, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [0, 0, 0, 0],
        [3075398820, 0, 0, 0],
        [0, 0, 0, 0],
        [-3075398820, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [4613098230, 0, 0, 0],
        [0, -1537699410, 0, 0],
        [0, 0, -1537699410, 0],
        [0, 0, 0, -1537699410]
    ]),
    Matrix(QQ, [
        [-1537699410, 0, 0, 0],
        [0, -1537699410, 0, 0],
        [0, 0, -1537699410, 0],
        [0, 0, 0, 4613098230]
    ]),
    Matrix(QQ, [
        [1537699410, 0, 0, 0],
        [0, -4613098230, 0, 0],
        [0, 0, 1537699410, 0],
        [0, 0, 0, 1537699410]
    ]),
    Matrix(QQ, [
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, -6150797640, 0, 0],
        [0, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 6150797640],
        [0, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, -6150797640, 0]
    ]),
    Matrix(QQ, [
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, -6150797640, 0, 0]
    ]),
    Matrix(QQ, [
        [-768849705, 6150797640, 0, 0],
        [0, 2306549115, 0, 0],
        [0, 0, -768849705, 0],
        [0, 3075398820, 0, -768849705]
    ]),
    Matrix(QQ, [
        [0, 0, 0, 0],
        [0, 0, -6150797640, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [0, 0, 6150797640, 0],
        [0, 0, 3075398820, 0],
        [0, 0, 0, 0],
        [0, 0, -3075398820, 0]
    ]),
    Matrix(QQ, [
        [0, 0, 0, 0],
        [0, 0, 0, 6150797640],
        [0, 0, 0, 0],
        [0, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [6150797640, 0, 0, 0],
        [0, 0, 0, 0]
    ]),
    Matrix(QQ, [
        [768849705, 0, 0, -6150797640],
        [0, 768849705, 0, -3075398820],
        [0, 0, 768849705, 0],
        [0, 0, 0, -2306549115]
    ])
    ]

    B = MatrixSpace(QQ,4,4)


    print("================================================== ")
    print("   START WITH THE Z-ORDER O OF  M4(Q) WITH Z-BASIS  : " )
    print("=================================================\n")

    print()
    print(Zbasis_O)
    print()

    print("================================================== ")
    print("   DEFINE WELL AN ORDER  : " )
    print("=================================================\n")

    print(is_order(B,Zbasis_O))

    print()
    print("================================================== ")
    print("   DISCRIMINANT OF O  : " )
    print("=================================================\n")
    d = discriminant(B,Zbasis_O)
    print()
    print(d)
    print("equal")
    print(factor(d))
    print()
    print("==========================================")
    print(" START COMPUTING A MAXIMAL ORDER CONTAINING O ")
    print("==========================================\n")
    print()
    Zbasis_L = max_order_containing_order(B,Zbasis_O)
    print("Done.")
    print()

    print("==========================================")
    print(" THE MAXIMAL ORDER IS L WITH Z BASIS      ")
    print("==========================================\n")
    print(Zbasis_L)

    print()
    print("==========================================")
    print(" CHECK     ")
    print("==========================================\n")
    print("does L contain O : " +str(is_sub_lattice(B,Zbasis_O,Zbasis_L)))
    print("discriminant L : " + str(factor(discriminant(B,Zbasis_L))))
    print()
