

SetSeed(12346);
p:=RandomPrime(73);
d:=-(24^2+p*53^2+p*27^2);
N := p*p*d;
M:= DiagonalMatrix(Rationals(),[1,p,p,d]);