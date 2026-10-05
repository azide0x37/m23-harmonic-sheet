#!/usr/bin/env python3
"""Eight-point test, stage 3: the third branch value b_3 of tau_1 and the exact
factorisation H_0(b_3, U) = q7(U) d8(U)^2 in L[U].

  (i)   reduce H_0 modulo the prime of L above 13 (y -> a in F_{13^6},
        sqrt(-23) -> 4): the U-discriminant locates b_3 mod 13 as the root of
        multiplicity >= 8, and H_0(b_3, U) factors as q7 d8^2 there;
  (ii)  Hensel-lift (b, q7, d8) in Z_13[y]/(f) to precision 13^N (quadratic
        Newton on an invertible 16x16 minor of the 23x16 Jacobian);
  (iii) recognise b_3, q7, d8 in L by LLL (lattice of dimension 13 + 6);
  (iv)  verify EXACTLY in L[U]: H_0(b_3,U) = q7 d8^2, d8 squarefree of degree
        8, gcd(q7, d8) = 1, gcd(d8, d/db H_0(b_3,U)) = 1, b_3 != 0 a 13-unit,
        q7 squarefree, c0 a 13-unit.
(iv) is the proof; (i)-(iii) only find the witness.   Usage: NPREC=N python3 stage3.py
"""
import cypari2, time, sys, os

N = int(os.environ.get("NPREC", "2000"))
pari = cypari2.Pari()
pari.allocatemem(6 * 10**9)
LOG = open("stage3.log", "w")
def log(*a):
    s = " ".join(str(t) for t in a)
    print(s); sys.stdout.flush(); LOG.write(s + "\n"); LOG.flush()
t0 = time.time()
def el(): return f"[{time.time()-t0:7.1f} s]"

pari('T = varhigher("T"); s = varhigher("s"); b = varhigher("b"); U = varhigher("U")')
pari('read("model-hat.gp"); default(seriesprecision, 60)')
pari(f'N = {N}')
pari('[H0junk, H, AU, UPM, e3junk] = read("stage2-H0.bin")')
pari('[DVEC, NVEC, c0] = read("stage1-DN.bin")')
pari('H0 = U^23 + sum(j = 0, 5, sum(k = 0, 23, if(H[j+1,k+1] != 0, H[j+1,k+1] * b^j * U^k, 0)))')
pari("ffa = ffgen(f*Mod(1,13), 'a)")
pari('read("stage3-lib.gp")')
log(f"== stage 3: 13-adic precision N = {N}; H_0 degrees (U, b) =", pari('[poldegree(H0, U), vecmax(vector(24, k, poldegree(polcoef(H0, k-1, U), b)))]'))

# (i)
r = pari('mod13_setup()')
log(el(), "mod 13: disc_U(H_0) degree", r[0], "val at b=0", r[1], "| multiplicity of b3bar:", r[2], "| factor (deg, mult):", r[3])
log(el(), "H_0(b3bar,U) mod 13 = q7 d8^2 with degrees (q7, d8) =", (r[4], r[5]), "| factors (deg, mult):", r[6])
assert int(r[4]) == 7 and int(r[5]) == 8

# (ii)
pari('hensel_setup()')
rk = pari('hensel_rows()')
log(el(), "Jacobian mod 13: rank", rk[0], "rows", rk[1])
assert int(rk[0]) == 16
for it in range(1, 40):
    v = int(pari('hensel_step()'))
    log(el(), f"Newton step {it}: residual valuation {v}")
    if v >= N: break
vfin = int(pari('Evalv(Eof(X))'))
log(el(), f"final residual valuation {vfin} (need {N})")
assert vfin >= N

# (iii)
pari('RB = recog(X[1]); b3 = RB[1]')
log(el(), "b_3 recognised: height 13^%.1f, margin 13^%.1f, digits %s" % (float(pari('RB[2]')), float(pari('RB[3]')), pari('#Str(lift(lift(b3)))')))
worst = float(pari('RB[3]'))
pari('q7 = U^7; d8 = U^8')
for i in range(7):
    pari(f'r = recog(X[2+{i}]); q7 += r[1]*U^{i}')
    m = float(pari('r[3]')); worst = min(worst, m)
    log(el(), f"  q7 coefficient {i}: height 13^{float(pari('r[2]')):.1f} margin 13^{m:.1f}")
for i in range(8):
    pari(f'r = recog(X[9+{i}]); d8 += r[1]*U^{i}')
    m = float(pari('r[3]')); worst = min(worst, m)
    log(el(), f"  d8 coefficient {i}: height 13^{float(pari('r[2]')):.1f} margin 13^{m:.1f}")
log(el(), f"worst recognition margin 13^{worst:.1f}")
if worst < 10:
    log("*** margin < 13^10: increase NPREC and rerun; exact verification skipped"); sys.exit(2)

# (iv) exact verification in L[U]
pari('HB = subst(H0, b, b3)')
c1 = bool(pari('HB - q7*d8^2 == 0')); log(el(), "[1] H_0(b_3, U) == q7 * d8^2 exactly in L[U]:", c1)
c2 = bool(pari('poldegree(gcd(d8, deriv(d8, U)), U) == 0 && poldegree(d8, U) == 8')); log(el(), "[2] d8 squarefree of degree 8:", c2)
c3 = bool(pari('poldegree(gcd(q7, d8), U) == 0')); log(el(), "[3] gcd(q7, d8) = 1:", c3)
pari('HBb = subst(deriv(H0, b), b, b3)')
c4 = bool(pari('poldegree(gcd(d8, HBb), U) == 0')); log(el(), "[4] gcd(d8, d/db H_0(b_3,U)) = 1 (the eight double roots are smooth points):", c4)
c5 = bool(pari('b3 != 0 && red13(b3) != 0')); log(el(), "[5] b_3 != 0 and b_3 is a 13-adic unit:", c5)
c6 = bool(pari('H[2,1] != 0 && vecmin(vector(23, k, H[1,k] == 0)) == 1 && H[6,1] == -1')); log(el(), "[6] H_0(0,U) = U^23, h_{1,0} != 0 (Eisenstein), h_{5,0} = -1:", c6)
c7 = bool(pari('poldegree(gcd(q7, deriv(q7, U)), U) == 0 && poldegree(q7, U) == 7')); log(el(), "[7] q7 squarefree of degree 7 (the fibre over b_3 has 15 distinct points):", c7)
c8 = bool(pari('red13(c0) != 0')); log(el(), "[8] c0 (tau_1 = c0 tau) is a 13-adic unit, so b_3/c0 is one too:", c8)
PASS = c1 and c2 and c3 and c4 and c5 and c6 and c7 and c8
log("==> EIGHT-POINT TEST:", "PASS -- the eight simple ramification points of tau lie over the single value b_3; tau has exactly three branch points" if PASS else "FAIL")
if not PASS: sys.exit(1)
for fn in ("eightpoint-witness.txt", "stage3-b3.bin"):
    if os.path.exists(fn): os.remove(fn)
pari(r'write("eightpoint-witness.txt", "\\ Eight-point witness. b3 in L = F(sqrt(-23)), f = y^6-y^5-3y^4-3y^3+y^2+5y+4, x = sqrt(-23); q7, d8 as coefficient vectors [c_0, ..., c_deg] (monic); H_0(b3, U) = q7(U) d8(U)^2 exactly in L[U].")')
pari('write("eightpoint-witness.txt", "b3 = ", b3)')
pari('write("eightpoint-witness.txt", "q7coeffs = ", vector(8, i, polcoef(q7, i-1)))')
pari('write("eightpoint-witness.txt", "d8coeffs = ", vector(9, i, polcoef(d8, i-1)))')
pari('writebin("stage3-b3.bin", [b3, vector(8, i, polcoef(q7, i-1)), vector(9, i, polcoef(d8, i-1))])')
log(el(), "written eightpoint-witness.txt / stage3-b3.bin")
