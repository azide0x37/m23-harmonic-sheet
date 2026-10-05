Fp := GF(7);
K<t> := FunctionField(Fp);
P<x> := PolynomialRing(K);
f := x^5 + t*x + 1;
F<z> := FunctionField(f);
G, r, S := GaloisGroup(F : Subfields := false);
print #G;
