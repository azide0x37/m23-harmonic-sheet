#!/usr/bin/env python3
"""(V1): the exact model over L is p-integral and reduces, under y -> a and
sqrt(-23) -> 4, to the characteristic-13 canonical model in the echelon
coordinates at P_0 that canon13.py computed from the MON13 plane model.

Usage: python3 check_V1.py K
"""
import sys, json

P = 13
FASC = [4, 5, 1, -3, -3, -1, 1]          # f = x^6 - x^5 - 3x^4 - 3x^3 + x^2 + 5x + 4


def reduce_elt(r):
    """(m(beta) + sqrt(-23) n(beta))/d  mod p, as F_13[a]/(f mod 13) coordinates"""
    d = r["d"] % P
    assert d, "denominator divisible by 13"
    dinv = pow(d, P - 2, P)
    return [(dinv * ((m + 4 * n) % P)) % P for m, n in zip(r["m"], r["n"])]


def main(K):
    rec = {r["label"]: r for r in json.load(open(f"reconCanon-{K}.json"))["items"]}
    st = json.load(open(f"canon-{K}.json"))
    items = [(f"q{a}{b}", st["qhat"][f"{a},{b}"]) for a in range(4) for b in range(a, 4)]
    items += [("C" + "".join(map(str, m)), v) for m, v in zip(st["cubic_monomials"], st["Chat"])]
    items += [(f"Pm{j}", v) for j, v in enumerate(st["Pminus"])]
    bad = [lab for lab, v in items if reduce_elt(rec[lab]) != [c % P for c in v]]
    dens = [lab for lab, _ in items if rec[lab]["d"] % P == 0]
    print(f"coefficients checked: {len(items)}")
    print(f"denominators divisible by 13: {dens if dens else 'none'}")
    print(f"reduction mismatches: {bad if bad else 'none'}")
    assert not bad and not dens
    print("(V1) verified")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 400)
