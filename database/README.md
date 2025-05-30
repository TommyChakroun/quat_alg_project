The data base db_iso_quat_alg... are made of block of the form

c
d
a
b
[[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]]

which mean that we have an isomorphism f : (c,d | Q) -> (a,b | Q) given by

f(1)=1
f(i)=x1*i+x2*j+x3*k
f(j)=y1*i+y2*j+y3*k
f(k)=z1*i+z2*j+z3*k


The data base db_iso_matrix_ring is made of block of the form 

c
d
a
b
[[..],..[..]] ( = M1  a matrix 4x4)
[[..],..[..]] ( = M2  a matrix 4x4)
..
[[..],..[..]] ( = M16 a matrix 4x4)

which mean that if A = (a,b|Q) and B = (c,d|Q)
we have an isomorphism phi : B ⊗ A^op -> M4(Q)
given by 
phi(ei⊗fj) = M_{4*i+j}
where 
e0,e1,e3,e4 = 1,i,j,k of B
f0,f1,f3,f4 = 1,i,j,k of A