#!/usr/bin/env python3
"""13-adic Hensel lift of the plane model H(tau,U) -- singularity formulation.

Unknowns (143):
  c_{i,j}   48   coefficients of H = U^23 + sum c_{i,j} tau^j U^{23-i}, 1<=j<=floor(5i/23),
                 with the L-rational normalisation c_{6,1} = c_{5,1}
  N         40   monic node polynomial N(tau) = tau^40 + ... (tau-coordinates of the nodes)
  u0        40   u0(tau) in W[tau]/(N): the u-coordinate of the node over each root of N
  q7         7   monic degree-7 (simple roots of H(1,U))
  d8         8   monic degree-8 (double roots of H(1,U))
Equations (143):
  H(tau,u0) = H_U(tau,u0) = H_tau(tau,u0) = 0 in W[tau]/(N)          3 x 40
  H(1,U) = q7(U) d8(U)^2                                               23
A solution near the MON13 model is a plane curve of bidegree (5,23) with the
cusp v^23 = sigma^5 at infinity, 40 singular points, and 8 double roots over
tau = 1: genus <= 4 by the singularities and >= 4 by Riemann-Hurwitz, hence a
genus-4 cover branched exactly over {0,1,infinity} with type (23; 2^8 1^7; 23).
Usage: python3 lift13b.py K
"""
import sys, json, time, os
sys.path.insert(0, "/home/claude/m23/relq"); sys.path.insert(0, "/home/claude/m23/mon13")
import flint
from wk import *
from lift13 import (FQ, FQF, FQR, RQ, FREE, start_point, mat_inverse, sum_mul, log, REL, NU as NU_OLD)

NC = 48
IDX = {"c": (0, 48), "N": (48, 88), "u0": (88, 128), "q7": (128, 135), "d8": (135, 143)}
NU = 143

def unpack(x, R):
    H = {}
    for k, v in zip(FREE, x[0:48]): H[k] = v
    H[(6, 1)] = H[(5, 1)]
    N = list(x[48:88]) + [R.one]
    u0 = list(x[88:128])
    q7 = list(x[128:135]) + [R.one]
    d8 = list(x[135:143]) + [R.one]
    return H, N, u0, q7, d8

# ---- polynomials in tau over R, reduced mod N (monic): plain lists, ascending
def pmod(R, p, N):
    return pdivmod(R, p, N)[1]

def H_coeffs_tau(R, H):
    """c_i(tau) for i = 0..23 as tau-polynomials (lists); c_0 = 1"""
    cs = [[] for _ in range(24)]
    cs[0] = [R.one]
    for (i, j), c in H.items():
        while len(cs[i]) <= j: cs[i].append(R.zero)
        cs[i][j] = R.add(cs[i][j], c)
    return [ptrim(R, c) for c in cs]

def eval_at_u0(R, cs, u0, N, weights=None):
    """sum_i w_i c_i(tau) u0^{23-i} mod N; weights w_i default 1 (H itself)."""
    # Horner in U: H = ((c_0 U + c_1) U + c_2) U + ...
    acc = []
    for i in range(24):
        term = cs[i] if weights is None else pscal(R, weights[i], cs[i])
        acc = padd(R, pmod(R, pmul(R, acc, u0), N), term)
    return pmod(R, acc, N)

def equations(R, x):
    H, N, u0, q7, d8 = unpack(x, R)
    cs = H_coeffs_tau(R, H)
    # H(tau,u0) mod N
    e1 = eval_at_u0(R, cs, u0, N)
    # H_U: sum (23-i) c_i U^{22-i}: use Horner on shifted list
    csU = [pscal(R, R.scal(23 - i, R.one), cs[i]) for i in range(23)]   # coefficients of U^{22-i}
    accU = []
    for i in range(23):
        accU = padd(R, pmod(R, pmul(R, accU, u0), N), csU[i])
    e2 = pmod(R, accU, N)
    # H_tau: derivative of c_i(tau)
    csT = [pderiv(R, c) for c in cs]
    e3 = eval_at_u0(R, csT, u0, N)
    # H(1,U) - q7 d8^2
    one = R.one
    h1 = [peval(R, cs[23 - e], one) if 23 - e < 24 else R.zero for e in range(24)]   # coefficient of U^e
    h1 = ptrim(R, h1)
    rhs = pmul(R, q7, pmul(R, d8, d8))
    e4 = psub(R, h1, rhs)
    def pad(p, n): return list(p) + [R.zero] * (n - len(p))
    return pad(e1, 40) + pad(e2, 40) + pad(e3, 40) + pad(e4, 23)

# ---- mod-13 starting point ---------------------------------------------------
class FqQuot:
    """F_{13^6}[tau]/(N) with fq_default_poly arithmetic"""
    def __init__(self, Nbar):
        self.N = Nbar; self.zero = RQ([]); self.one = RQ([FQ.one()]); self.two = RQ([FQ(2)])
    def add(self, x, y): return x + y
    def sub(self, x, y): return x - y
    def neg(self, x): return -x
    def mul(self, x, y): return (x * y) % self.N
    def eq(self, x, y): return x == y
    def is_zero(self, x): return x.is_zero()
    def is_unit(self, x): return x.gcd(self.N).degree() == 0
    def inv(self, x):
        g, s, t = x.xgcd(self.N)
        assert g.degree() == 0
        return (s / g[0]) % self.N
    def scal(self, n, x): return x * FQ(n)

def start_point_b():
    H, Nbar, Cbar, x0_old, pts = start_point()
    A = FqQuot(Nbar)
    # H(tau,U) with tau the class of tau in F[tau]/(N): U-polynomial over A
    tauA = RQ([FQ.zero(), FQ.one()]) % Nbar
    tp = [A.one]
    for _ in range(5): tp.append(A.mul(tp[-1], tauA))
    poly = [A.zero] * 24; poly[23] = A.one
    for (i, j), c in H.items():
        poly[23 - i] = A.add(poly[23 - i], A.mul(RQ([c]), tp[j]))
    dpoly = pderiv(A, poly)
    # Euclid for gcd over A
    p, q = poly[:], dpoly[:]
    while True:
        if not A.is_unit(q[-1]): raise RuntimeError("non-unit leading coefficient in node gcd")
        _, r = pdivmod(A, p, q)
        if not r: break
        p, q = q, r
    assert len(q) == 2, f"gcd over the node algebra has degree {len(q)-1}, expected 1"
    u0 = A.neg(A.mul(q[0], A.inv(q[1])))            # root of the linear gcd
    # H(1,U) = q7 d8^2 over F_{13^6}
    h1 = RQ([FQ.zero()] * 24)
    cs = H_coeffs_tau(FQR, H)
    h1 = RQ([peval(FQR, cs[23 - e], FQ.one()) for e in range(24)])
    fac = h1.factor()
    d8 = RQ([FQ.one()]); q7 = RQ([FQ.one()])
    for f, e in fac[1]:
        if e == 2: d8 = d8 * f
        elif e == 1: q7 = q7 * f
        else: raise RuntimeError(f"multiplicity {e} over tau=1")
    assert d8.degree() == 8 and q7.degree() == 7 and h1 == fac[0] * q7 * d8 * d8 and fac[0] == FQ.one()
    x0 = [H[k] for k in FREE] + [Nbar[i] for i in range(40)] + [u0[i] for i in range(40)] \
         + [q7[i] for i in range(7)] + [d8[i] for i in range(8)]
    # sanity: equations vanish mod 13
    eq = equations(FQR, x0)
    assert all(v.is_zero() for v in eq), "start point does not satisfy the equations"
    return x0, H, Nbar, u0, q7, d8

def jacobian(R, x, dual):
    cols = []
    base = [dual.const(v) for v in x]
    for ci in range(NU):
        xd = base[:]; xd[ci] = dual.eps(x[ci], R.one)
        cols.append([e[1] for e in equations(dual, xd)])
    return [[cols[c][r] for c in range(NU)] for r in range(NU)]

class WkDual:
    def __init__(self, W):
        self.W = W; self.zero = (W.zero, W.zero); self.one = (W.one, W.zero); self.two = (W.two, W.zero)
    def const(self, v): return (v, self.W.zero)
    def eps(self, v, d): return (v, d)
    def add(self, x, y): return (x[0] + y[0], x[1] + y[1])
    def sub(self, x, y): return (x[0] - y[0], x[1] - y[1])
    def neg(self, x): return (-x[0], -x[1])
    def mul(self, x, y): return (self.W.mul(x[0], y[0]), self.W.mul(x[0], y[1]) + self.W.mul(x[1], y[0]))
    def eq(self, x, y): return x[0] == y[0] and x[1] == y[1]
    def is_zero(self, x): return x[0].is_zero() and x[1].is_zero()
    def is_unit(self, x): return self.W.is_unit(x[0])
    def inv(self, x):
        iv = self.W.inv(x[0]); return (iv, -self.W.mul(self.W.mul(x[1], iv), iv))
    def scal(self, n, x): return (self.W.scal(n, x[0]), self.W.scal(n, x[1]))

def rank_fq(M):
    basis = []
    for row in M:
        v = row[:]
        for pc, b in basis:
            if not v[pc].is_zero():
                f = v[pc]; v = [x - f * y for x, y in zip(v, b)]
        pc = next((i for i, t in enumerate(v) if not t.is_zero()), None)
        if pc is None: continue
        inv = 1 / v[pc]; basis.append((pc, [t * inv for t in v]))
    return len(basis)

def hensel(Ktarget, checkpoint):
    log("setup (singularity formulation)")
    x0, H, Nbar, u0, q7, d8 = start_point_b()
    log("mod-13 start point ok: 40 nodes with linear gcd, H(1,U) = q7 d8^2")
    J = jacobian(FQR, x0, FqDual(FQ))
    rk = rank_fq(J)
    log(f"jacobian mod 13: rank {rk} / {NU}")
    assert rk == NU, "singular Jacobian mod 13"
    Jinv = mat_inverse(FQR, J)
    log("jacobian inverted")
    W = Wk(Ktarget)
    import glob
    cands = [checkpoint] + sorted(glob.glob(f"{REL}/liftb-*.json"), key=lambda p: -json.load(open(p))["prec"])
    resumed = False
    for cp in cands:
        if os.path.exists(cp):
            st = json.load(open(cp))
            if st["prec"] <= Ktarget:
                x = [W.from_asc(v) for v in st["x"]]; n = st["prec"]; resumed = True
                log(f"resumed from {cp} at precision 13^{n}")
                break
    if not resumed:
        x = [W.from_fq(v) for v in x0]; n = 1
    pn = P ** n
    while n < Ktarget:
        t0 = time.time()
        res = equations(W, x)
        r = []
        for c in res:
            asc = W.asc(c)
            assert all(v % pn == 0 for v in asc), f"residual not zero mod 13^{n}"
            r.append(FQ(flint.fmpz_mod_poly_ctx(P)([(v // pn) % P for v in asc])))
        delta = [-sum_mul(FQR, Jinv[i], r) for i in range(NU)]
        x = [W.add(xi, W.scal(pn, W.from_fq(d))) for xi, d in zip(x, delta)]
        n += 1; pn *= P
        if n % 10 == 0 or n == Ktarget:
            log(f"precision 13^{n}  ({time.time()-t0:.2f}s/step)")
            json.dump({"prec": n, "K": Ktarget, "x": [W.asc(v) for v in x], "free": FREE, "idx": IDX},
                      open(checkpoint, "w"))
    res = equations(W, x)
    assert all(W.is_zero(c) for c in res), "final residual nonzero"
    log(f"lift verified: all {NU} equations hold mod 13^{Ktarget}")
    json.dump({"prec": n, "K": Ktarget, "x": [W.asc(v) for v in x], "free": FREE, "idx": IDX}, open(checkpoint, "w"))
    return W, x

if __name__ == "__main__":
    Kt = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    hensel(Kt, f"{REL}/liftb-{Kt}.json")
