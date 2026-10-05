SetSeed(1);
SetVerbose("GaloisGroup", 5);
SetVerbose("Invariant", 3);
G := Alt(23);
Hdb := TransitiveGroup(23, 5);
H := sub<G | [ G!h : h in Generators(Hdb) ]>;
assert #H eq 10200960;
assert IsMaximal(G, H);
printf "GROUPS A23_order=%o M23_order=%o index=%o\n", #G, #H, Index(G,H);
for wl in [0..8] do
    printf "WORKLEVEL_BEGIN %o\n", wl;
    try
        time inv := GaloisGroupInvariant(G, H : Worklevel := wl);
        if Type(inv) eq BoolElt then
            printf "WORKLEVEL_RESULT %o false\n", wl;
        else
            printf "WORKLEVEL_RESULT %o constructed rank=%o\n", wl, Rank(Parent(inv));
            print PrettyPrintInvariant(inv);
            assert &and[ IsInvariant(inv, h) : h in Generators(H) ];
        end if;
    catch err
        printf "WORKLEVEL_ERROR %o %o\n", wl, err`Object;
    end try;
end for;
print "WORKLEVEL_SCAN_COMPLETE";
