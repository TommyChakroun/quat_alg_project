load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("maximal_orders/find_maximal_orders.sage")

#------------------------------------------------------------------------------------------
#
#            TEST FIND MAXIMAL ORDER in C =  A ⊗ B^op
#               when A,B quaternio algebras and we know O1 in A, O2 in B max order
#
#------------------------------------------------------------------------------------------




def test_max_order_tensor_quat_alg(A,B,Zbasis_O1=None,Zbasis_O2=None):
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
    print("=" * 80)
    print("STARTING test_max_order_tensor_quat_alg")
    print("=" * 80)

    print()
    print(f"Input A: {A}")
    print(f"Input B: {B}")
    print(f"Input Zbasis_O1: {Zbasis_O1}")
    print(f"Input Zbasis_O2: {Zbasis_O2}")
    print()
    
    print()
    print("Creating tensor product C = tensor(A, opposite(B))...")
    C = tensor(A, opposite(B))
    print(f"C = {C}")
    print()
    
    print()
    print("Getting basis of C...")
    BC = C.basis()
    print(f"BC = {BC}")
    print(f"Length of BC: {len(BC)}")
    print()
    
    print()
    if Zbasis_O1 is None:
        print("=" * 80)
        print("Zbasis_O1 is None, computing a maximal order in A...")
        print("=" * 80)
        print()
        Zbasis_O1 = max_order(A)
        print(f"Computed Zbasis_O1: {Zbasis_O1}")
    else:
        print(f"Using provided Zbasis_O1: {Zbasis_O1}")
    print()
    
    print()
    if Zbasis_O2 is None:
        print("=" * 80)
        print("Zbasis_O2 is None, computing a maximal order in B...")
        print("=" * 80)
        Zbasis_O2 = max_order(B)
        print(f"Computed Zbasis_O2: {Zbasis_O2}")
    else:
        print(f"Using provided Zbasis_O2: {Zbasis_O2}")
    print()

    print()
    print(f"Zbasis_O1 has {len(Zbasis_O1)} elements")
    print(f"Zbasis_O2 has {len(Zbasis_O2)} elements")
    print(f"Expected tensor basis size: {len(Zbasis_O1)} × {len(Zbasis_O2)} = {len(Zbasis_O1) * len(Zbasis_O2)}")
    print()
    
    print()
    print("=" * 80)
    print("Building a Z basis of O = O1⊗O2 in C")
    print("=" * 80)
    print()
    Zbasis_O = []
    
    for e in Zbasis_O1:
        s = get_coefficients(e, A)
        for f in Zbasis_O2:
            r = get_coefficients(f, B)
            e_tensor_f = sum(sum(s[k]*r[l]*BC[4*k+l] for k in range(4)) for l in range(4))
            Zbasis_O.append(e_tensor_f)

    print("Zbasis_O = ")
    print(Zbasis_O)


    print()
    print("=" * 80)
    print("Discriminant Check :")
    print("=" * 80)
    dO1 = discriminant(A, Zbasis_O1)
    dO2 = discriminant(B, Zbasis_O2)
    dO = discriminant(C, Zbasis_O)
    print(f"discriminant(A, O1) = {dO1}")
    print(f"discriminant(B, O2) = {dO2}")
    print(f"discriminant(C, O) = {dO}")
    print(f"Expected discriminant relation: dO = dO1**4* dO2**4 = {dO1**4 * dO2**4}")
    print(f"Discriminant check: {dO == dO1**4 * dO2**4}")
    print()


    print()
    print("=" * 80)
    print("Computing maximal order L containing O...")
    print("=" * 80)
    Zbasis_L = max_order_containing_order(C, Zbasis_O,printers = True)
    print()
    print(f"Zbasis_L = {Zbasis_L}")
    print(f"Zbasis_L has {len(Zbasis_L)} elements")

    print()
    print("=" * 80)
    print("Final discriminant check :")
    print("=" * 80)
    print()
    dL = discriminant(C, Zbasis_L)
    print(f"dGamma = {dL} = {factor(dL)}")
    if dL != 1 and dL != -1:
        print("The discriminant is not ±1 so A and B are not isomorphic and A⊗B^op is not isomorphic to M4(Q).")
        print("Hence a maximal order has not necessarily discriminant 1, so it is not a problem")

    print()
    print("=" * 80)
    print("COMPLETED test_max_order_tensor_quat_alg")
    print("=" * 80)
    
    return Zbasis_L


def test_15(a,b,c,d):
    A = QuaternionAlgebra(QQ,a,b)
    B = QuaternionAlgebra(QQ,c,d)
    test_max_order_tensor_quat_alg(A,B)
    return None