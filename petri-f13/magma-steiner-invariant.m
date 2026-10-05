SetSeed(1);
G := Alt(23);
Hdb := TransitiveGroup(23, 5);
H := sub<G | [ G!h : h in Generators(Hdb) ]>;
assert #H eq 10200960;
assert IsMaximal(G, H);
classes := ConjugacyClasses(H);
involutions := [ c[3] : c in classes | Order(c[3]) eq 2 ];
assert #involutions eq 1;
tau := involutions[1];
B0 := { i : i in [1..23] | i^tau eq i };
assert #B0 eq 7;
blocks := B0^H;
assert #blocks eq 253;
four_counts := AssociativeArray();
for B in blocks do
    for S in Subsets(B, 4) do
        key := Sprint(S);
        if not IsDefined(four_counts, key) then
            four_counts[key] := 0;
        end if;
        four_counts[key] +:= 1;
    end for;
end for;
assert #Keys(four_counts) eq Binomial(23,4);
assert &and[ four_counts[key] eq 1 : key in Keys(four_counts) ];
D := Design< 4, 23 | Setseq(blocks) >;
A := AutomorphismGroup(D);
assert #A eq #H;
R := SLPolynomialRing(Integers(), 23);
J := &+[ (&+[ R.i : i in B ])^5 : B in blocks ];
assert &and[ IsInvariant(J, h) : h in Generators(H) ];
cost, make := GaloisGroupInvariant(G, H : DoCost := true, Worklevel := 5);
I := make();
assert &and[ IsInvariant(I, h) : h in Generators(H) ];
printf "STEINER seed=%o blocks=%o four_subsets=%o aut_order=%o\n",
       Setseq(B0), #blocks, #Keys(four_counts), #A;
printf "BUILTIN worklevel=5 cost=%o rank=%o\n", cost, Rank(Parent(I));
printf "EXPLICIT rank=%o formula=sum_253_(sum_heptad_x)^5\n", Rank(Parent(J));
print "BUILTIN_PRETTY_BEGIN";
print PrettyPrintInvariant(I);
print "BUILTIN_PRETTY_END";
print "STEINER_INVARIANT_CONSTRUCTED";
