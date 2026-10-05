Fp := GF(7);
K<t> := FunctionField(Fp);
P<x> := PolynomialRing(K);
f := x^5 + t*x + 1;
q := t^2 + t + 3;
NP := function(old)
    return old + 1;
end function;
G, r, S := GaloisGroup(f : Prime := q, NextPrime := NP);
print #G;
