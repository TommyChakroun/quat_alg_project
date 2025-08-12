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
