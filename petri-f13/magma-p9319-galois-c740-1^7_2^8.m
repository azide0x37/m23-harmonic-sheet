p := 9319;
k := GF(p);
Kt<T> := FunctionField(k);
P<w> := PolynomialRing(Kt);
h := w^23 + (65*T^3 + 2851*T^2 + 2129*T + 3075)/(T^3 + 564*T^2 + 5046*T + 4384)*w^22
    + (2831*T^3 + 3142*T^2 + 6646*T + 6071)/(T^3 + 564*T^2 + 5046*T + 4384)*w^21
    + (5763*T^3 + 4924*T^2 + 8079*T + 447)/(T^3 + 564*T^2 + 5046*T + 4384)*w^20
    + (5261*T^3 + 5161*T^2 + 9013*T + 8484)/(T^3 + 564*T^2 + 5046*T + 4384)*w^19
    + (6257*T^3 + 4335*T^2 + 7218*T + 3328)/(T^3 + 564*T^2 + 5046*T + 4384)*w^18
    + (7669*T^3 + 3939*T^2 + 2117*T + 8050)/(T^3 + 564*T^2 + 5046*T + 4384)*w^17
    + (6561*T^3 + 3315*T^2 + 2027*T + 8365)/(T^3 + 564*T^2 + 5046*T + 4384)*w^16
    + (4670*T^3 + 1651*T^2 + 8256*T + 5626)/(T^3 + 564*T^2 + 5046*T + 4384)*w^15
    + (7140*T^3 + 8193*T^2 + 3562*T + 3629)/(T^3 + 564*T^2 + 5046*T + 4384)*w^14
    + (5965*T^3 + 4027*T^2 + 371*T + 2513)/(T^3 + 564*T^2 + 5046*T + 4384)*w^13
    + (3421*T^3 + 546*T^2 + 7710*T + 1810)/(T^3 + 564*T^2 + 5046*T + 4384)*w^12
    + (7237*T^3 + 1149*T^2 + 8136*T + 5990)/(T^3 + 564*T^2 + 5046*T + 4384)*w^11
    + (5961*T^3 + 6051*T^2 + 1480*T + 6074)/(T^3 + 564*T^2 + 5046*T + 4384)*w^10
    + (5203*T^3 + 6793*T^2 + 9056*T + 270)/(T^3 + 564*T^2 + 5046*T + 4384)*w^9 +
    (719*T^3 + 8006*T^2 + 2244*T + 5557)/(T^3 + 564*T^2 + 5046*T + 4384)*w^8 +
    (6891*T^3 + 6844*T^2 + 1140*T + 2235)/(T^3 + 564*T^2 + 5046*T + 4384)*w^7 +
    (2132*T^3 + 4976*T^2 + 3493*T + 8287)/(T^3 + 564*T^2 + 5046*T + 4384)*w^6 +
    (1063*T^3 + 5506*T^2 + 5619*T + 2574)/(T^3 + 564*T^2 + 5046*T + 4384)*w^5 +
    (4757*T^3 + 6938*T^2 + 5704*T + 7692)/(T^3 + 564*T^2 + 5046*T + 4384)*w^4 +
    (7956*T^3 + 5348*T^2 + 2028*T + 6987)/(T^3 + 564*T^2 + 5046*T + 4384)*w^3 +
    (1765*T^3 + 4928*T^2 + 1730*T + 110)/(T^3 + 564*T^2 + 5046*T + 4384)*w^2 +
    (2774*T^3 + 8013*T^2 + 7318*T + 2416)/(T^3 + 564*T^2 + 5046*T + 4384)*w +
    (1803*T^4 + 7839*T^3 + 3925*T^2 + 8995*T + 8465)/(T^4 + 564*T^3 + 5046*T^2 +
    4384*T);
assert Degree(h) eq 23;
assert IsMonic(h);
assert IsIrreducible(h);
assert IsSeparable(h);
print "EXACT_POLYNOMIAL_OK";

SetSeed(1);
SetVerbose("GaloisGroup", 5);
SetVerbose("Invariant", 3);
prime := T - k!740;
printf "PRIME %o\n", prime;
assert &and[Evaluate(Denominator(Coefficient(h,i)),k!740) ne 0 : i in [0..23]];
hc := Polynomial([Evaluate(Numerator(Coefficient(h,i)),k!740) /
                  Evaluate(Denominator(Coefficient(h,i)),k!740)
                  : i in [0..23]]);
assert Discriminant(hc) ne 0;
print "PRIME_PATTERN", Sort([Degree(v[1]) : v in Factorization(hc)]);
procedure CertifiedRun()
    time G, roots, data := GaloisGroup(h : Prime := prime, ShortOK := false);
    print "GROUP_ORDER", #G;
    print "GROUP_DEGREE", Degree(G);
    print "GROUP_TRANSITIVE", IsTransitive(G);
    print "GROUP_GENERATORS", Generators(G);
    print "CERTIFIED_GALOISGROUP_RETURNED";
end procedure;
CertifiedRun();
