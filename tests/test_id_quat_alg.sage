load("utilities/utilities.sage")
load("utilities/algebra_type.sage")
load("core/identify_quaternion_algebra.sage")

#------------------------------------------------------------------------------------------
#
#            TEST OF IDENTIFY QUATERNION ALGEBRAS
#
#------------------------------------------------------------------------------------------

def test_1():
    B  = quat_alg_mixed_with_table(QQ,-1,-1)
    e0,e1,e2,e3 = B.basis()
    i,j,a,b = quaternion_structure(B)


    print("B :")
    print(B)
    print("B is isomorphic to quaternion algebra :" + str(is_isomorphic_to_quaternion_algebra(B)))
    print("quaternion structure (i,j,a,b) = " + str((i,j,a,b)))
    print("CHECK i**2=a, j**2=b, i*j=-j*i :")
    print(i**2==a*e0)
    print(j**2==b*e0)
    print(i*j==-j*i)


    B  = quat_alg_mixed_with_table(QQ,-2,-1)
    e0,e1,e2,e3 = B.basis()
    i,j,a,b = quaternion_structure(B)


    print("B :")
    print(B)
    print("B is isomorphic to quaternion algebra :" + str(is_isomorphic_to_quaternion_algebra(B)))
    print("quaternion structure (i,j,a,b) = " + str((i,j,a,b)))
    print("CHECK i**2=a, j**2=b, i*j=-j*i :")
    print(i**2==a*e0)
    print(j**2==b*e0)
    print(i*j==-j*i)


    B  = quat_alg_mixed_with_table(QQ,-5,-7)
    e0,e1,e2,e3 = B.basis()
    i,j,a,b = quaternion_structure(B)


    print("B :")
    print(B)
    print("B is isomorphic to quaternion algebra :" + str(is_isomorphic_to_quaternion_algebra(B)))
    print("quaternion structure (i,j,a,b) = " + str((i,j,a,b)))
    print("CHECK i**2=a, j**2=b, i*j=-j*i :")
    print(i**2==a*e0)
    print(j**2==b*e0)
    print(i*j==-j*i)


def test_2():
    x = read_file(5,"/quat_alg_project/database/db_iso_quat_alg.txt")
    for c,d,a,b,R in x:
        A = QuaternionAlgebra(QQ,a,b)
        B = QuaternionAlgebra(QQ,c,d)
        one_A,i,j,k = A.basis()
        isom_dict={}
        isom_dict[0] = one_A
        isom_dict[1] = R[0][0]*i+R[0][1]*j+R[0][2]*k
        isom_dict[2] = R[1][0]*i+R[1][1]*j+R[1][2]*k
        isom_dict[3] = R[2][0]*i+R[2][1]*j+R[2][2]*k
        print("iso : " + str(is_linear_iso(isom_dict,B,A)))
        print("alg homo : " +str(is_algebra_homomorphism(isom_dict, B, A)))
        print(" i**2 ==c : " +str(isom_dict[1]**2==c*one_A))