Q := Rationals();

// We define two isomorphic quaternion algebras over the rationals

a := -124545653451326554;
b := -545646534564633;

alpha := -73;
beta := -69531513;

A := QuaternionAlgebra(Q,a,b);
B :=QuaternionAlgebra(Q,alpha, beta);

 
K<w> := QuadraticField(alpha);

// Since K split B, K split also A
// Hence equivalently,
// (norm  equation on K(sqrt(a))/K )   X^2-aY^2=b has a solution X,Y in K       
// (conic equation on K)               X^2-a Y^2 - b Z ^2 = 0 has a non zero solution X,Y,Z in K.


// We are going to compare which magma function is better for this task
// Appearently solving the conic equation is faster than solving the norm equation because use reduction before to an easier norm equation.

C := Conic([K|1,-a,-b]);


P<x> := PolynomialRing(K);
L := ext<K | x^2 - a>;



                                                                                                                                                                            
