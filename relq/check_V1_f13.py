#!/usr/bin/env python3
"""(V1), replayable from this repository: the exact model over L in model-hat.gp is
p-integral and reduces, under y -> a and sqrt(-23) -> 4 (the place p of L above 13 with
residue field F_{13^6} = F_13[a]/(f mod 13)), to the exact characteristic-13 canonical
model of ../petri-f13/exact-f13-map.json (reduced_model: quadric, cubic, P0, Pinf),
in the echelon coordinates at P_0, with P_+ -> P_0 and P_- -> P_inf.

The cubic of a canonical genus-4 curve is defined modulo x_i*Q and scaling; both cubics
are brought to the normal form of canon13.py (no monomial divisible by x_2^2, coefficient
1 on x_1^2 x_2) before comparison.  Points are compared projectively.
Usage: python3 check_V1_f13.py            (needs cypari2)
"""
import json, sys, cypari2

pari = cypari2.Pari()
pari('read("model-hat.gp")')
pari("ffa = ffgen(f*Mod(1,13), 'a)")
d = json.load(open("../petri-f13/exact-f13-map.json"))
fld = d["field"]
assert fld["characteristic"] == 13 and fld["degree"] == 6 and fld["s_residue"] == 4
assert fld["compact_polynomial_coefficients_descending_mod_13"] == [1, 12, 10, 10, 1, 5, 4]   # f mod 13
rm = d["reduced_model"]

def ffvec(v):            # [c5, ..., c0] high-to-low -> element of F_{13^6}
    return pari("sum(i=1,6, %s[i]*ffa^(6-i))" % str(list(v)))

pari('dens(t) = if(type(t) == "t_POL", lcm(vector(poldegree(t)+1, i, dens(polcoef(t,i-1)))), denominator(t))')
pari('red13(e) = {my(l = lift(lift(e)), d, num); if(l == 0, return(0*ffa)); d = dens(l); if(d % 13 == 0, error("denominator divisible by 13")); num = l*d; subst(subst(num, x, 4), y, ffa) / d}')
def red13(gp_expr):      # reduce an element of L (nested POLMOD) at p
    return pari("red13(%s)" % gp_expr)

bad = []
# ---- quadric
for i in range(4):
    for j in range(i, 4):
        lhs = red13("QUAD[%d,%d]" % (i+1, j+1))
        rhs = ffvec(rm["quadric"]["(%d, %d)" % (i, j)])
        if lhs != rhs: bad.append("q%d%d" % (i, j))
print("quadric coefficients agree mod p:", "all 10" if not [b for b in bad if b[0]=='q'] else bad)

# ---- cubic, in normal form modulo x_i Q
# monomial dictionaries: key = sorted index triple -> F_{13^6}
def key3(m):             # exponent vector (e0,e1,e2,e3) -> "(i, j, k)"
    idx = []
    for i, e in enumerate(m): idx += [i]*e
    return "(%d, %d, %d)" % tuple(idx)
cub_hat = {}
for t in range(20):
    m = [int(v) for v in pari("CUBM[%d]" % (t+1))]
    cub_hat[key3(m)] = red13("CUBC[%d]" % (t+1))
cub_f13 = {k: ffvec(v) for k, v in rm["cubic"].items()}
quad_f13 = {}
for i in range(4):
    for j in range(i, 4):
        quad_f13[(i, j)] = ffvec(rm["quadric"]["(%d, %d)" % (i, j)])
assert quad_f13[(2, 2)] == 1, "q22 != 1 in the reduced model"

def normal_form(cub):
    """kill the monomials x_2^2 x_i (i = 0..3) with x_i * Q (q22 = 1), then scale x_1^2 x_2 to 1."""
    c = dict(cub)
    zero = pari("0*ffa")
    def get(k): return c.get(k, zero)
    for i in (2, 0, 1, 3):      # x_2 Q first: it also carries x_2^2 x_0, x_2^2 x_1, x_2^2 x_3
        k = [0,0,0,0]; k[2] += 2; k[i] += 1; k = key3(k)
        coef = get(k)
        if coef == 0: continue
        # subtract coef * x_i * Q
        for (a, b), q in quad_f13.items():
            m = [0,0,0,0]; m[a] += 1; m[b] += 1; m[i] += 1
            kk = key3(m)
            c[kk] = get(kk) - coef * q
    for i in range(4):
        k = [0,0,0,0]; k[2] += 2; k[i] += 1
        assert get(key3(k)) == 0, "elimination failed"
    lead = get("(1, 1, 2)")
    assert lead != 0, "x_1^2 x_2 coefficient is zero"
    return {k: v / lead for k, v in c.items()}
nh = normal_form(cub_hat); nf = normal_form(cub_f13)
keys = set(nh) | set(nf)
zero = pari("0*ffa")
cb = [k for k in keys if nh.get(k, zero) != nf.get(k, zero)]
print("cubic coefficients agree mod p (normal form modulo x_i Q):", "all 20" if not cb else cb)
bad += cb

# ---- points
pm = [red13("PM[%d]" % (i+1)) for i in range(4)]
pinf = [ffvec(v) for v in rm["Pinf"]]
# both projective: scale pinf so that its x0 = 1 (pm has x0 = 1)
assert pinf[0] != 0
pinf = [v / pinf[0] for v in pinf]
pb = [i for i in range(4) if pm[i] != pinf[i]]
print("P_- reduces to P_inf (projectively):", "yes" if not pb else ("NO, coordinates %s differ" % pb))
bad += ["Pm%d" % i for i in pb]
p0 = [ffvec(v) for v in rm["P0"]]
print("P_+ = (1:0:0:0) = P_0:", p0[0] != 0 and all(v == 0 for v in p0[1:]))
print("(V1) verified" if not bad else "(V1) FAILED: %s" % bad)
sys.exit(0 if not bad else 1)
