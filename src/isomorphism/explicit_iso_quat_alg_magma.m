IsoQuatAlg := function(A, B)

    Q := BaseRing(A);
    alpha, beta := -Norm(A.1),-Norm(A.2);
    a, b := -Norm(B.1),-Norm(B.2);

    D := DiagonalMatrix(Q, [a, b, -a*b, -alpha]);
   
    sol := Basis(IsotropicSubspace(D))[1];

    mu := (sol[1]/sol[4])*B.1 + (sol[2]/sol[4])*B.2 + (sol[3]/sol[4])*B.3;

    nu := QuaternionicComplement(mu);
    gamma := -Norm(nu);

    D2 := DiagonalMatrix(Q, [1, -alpha, -beta/gamma]);
    sol2 := Basis(IsotropicSubspace(D2))[1];
    
    x := sol2[1]/sol2[3];
    y := sol2[2]/sol2[3];
    
    return [mu, (x + y*mu)*nu];

end function;


Q := Rationals();


// TESTS

// A0
p0 := NextPrime(10);
A0 := QuaternionAlgebra(Q,-11,-1);
i,j,k := A0.1,A0.2,A0.3;
mu0 := 3*i + 7*j + 2*k;
nu0 := 2*i + 3*k;
a0 := -Norm(mu0);
b0 := -Norm(nu0);
B0 := QuaternionAlgebra(Q,-192,-143);

// A1
p1 := NextPrime(100);
A1 := QuaternionAlgebra(Q,-101,-3);
i,j,k := A1.1,A1.2,A1.3;
mu1 := 15*i + 83*j + 21*k;
nu1 := 21*i - 5*k;
a1 := -Norm(mu1);
b1 := -Norm(nu1);
B1 := QuaternionAlgebra(Q,-177015,-52116);

// A2
p2 := NextPrime(1000);
A2 := QuaternionAlgebra(Q,-1009,-11);
i,j,k := A2.1,A2.2,A2.3;
mu2 := 312*i + 101*j + 542*k;
nu2 := 2981*i - 156*k;
a2 := -Norm(mu2);
b2 := -Norm(nu2);
B2 := QuaternionAlgebra(Q,-3358818943, -9236443513);

// A3
p3 := NextPrime(2^10);
A3 := QuaternionAlgebra(Q,-1031,-1);
i,j,k := A3.1,A3.2,A3.3;
mu3 := 1052*i + 2771*j + 985*k;
nu3 := 985*i - 1052*k;
a3 := -Norm(mu3);
b3 := -Norm(nu3);
B3 := QuaternionAlgebra(Q,-2148992240, -2141313799);

// A4
p4 := NextPrime(2^16);
A4 := QuaternionAlgebra(Q,-65537,-3);
i,j,k := A4.1,A4.2,A4.3;
mu4 := 12512*i + 42188*j + 31854*k;
nu4 := 47781*i - 6256*k;
a4 := -Norm(mu4);
b4 := -Norm(nu4);
B4 := QuaternionAlgebra(Q,-209761888045436, -157317411422553);

// A5 
p5 := NextPrime(2^28);
A5 := QuaternionAlgebra(Q,-268435459,-1);
i,j,k := A5.1,A5.2,A5.3;
mu5 := 25156*i + 5446*j + 54*k;
nu5 := 27*i-12578*k;
a5 := -Norm(mu5);
b5 := -Norm(nu5);
B5 := QuaternionAlgebra(Q,-169873273887987584,-42468318464582167);
