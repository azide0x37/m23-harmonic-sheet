#!/usr/bin/env python3
"""The canonical model of the lifted cover in the intrinsic echelon coordinates at P_+.

Coordinates: x_j = omega-hat_j, the echelon basis of H^0(K) w.r.t. the adapted parameter
w-hat at P_+ scaled so that the Petri quadric has q_12 = q_22 = 1.  Then:
  Q-hat : the Petri quadric (q_00 = q_01 = 0, q_02 = -q_11, q_03 = -q_12 = -1),
  C-hat : the Petri cubic with no monomial divisible by x_2^2 and coefficient 1 on x_1^2 x_2,
  D-hat : quintic vanishing to order 23 at P_+ = (1:0:0:0), N-hat : quintic vanishing to
          order 23 at P_-, N-hat / D-hat = tau, normalised by N-hat(P_+) = 1,
  P_-   : coordinates of the other wild point.
All are L-rational and are written (in W_K) to canon-K.json.   Usage: python3 canon13.py K
"""
import sys, json, time
sys.path.insert(0, "/home/claude/m23/relq"); sys.path.insert(0, "/home/claude/m23/mon13")
from wk import *
from lift13 import FQ, FQR, FREE, log, REL, mat_inverse
from lift13b import unpack, rank_fq
from petri13 import nullspace_W
from petri13L import MON, branch_L, monomial_series, expand, petri_at, wpow

def series_pow_list(W, base_series, e):
    out = LSer.const(W, W.one, n=base_series.n + 50)
    for _ in range(e): out = out * base_series
    return out

def solve_overdet(W, rows, rhs, n):
    """solve the (overdetermined, consistent) system rows * z = rhs by unit-pivot elimination"""
    aug = [r[:] + [b] for r, b in zip(rows, rhs)]
    piv = []
    for c in range(n):
        i = next((i for i in range(len(aug)) if i not in [p for p, _ in piv] and W.is_unit(aug[i][c])), None)
        assert i is not None, f"no unit pivot in column {c}"
        inv = W.inv(aug[i][c]); aug[i] = [W.mul(v, inv) for v in aug[i]]
        for j in range(len(aug)):
            if j == i or W.is_zero(aug[j][c]): continue
            f = aug[j][c]; aug[j] = [W.sub(a, W.mul(f, b)) for a, b in zip(aug[j], aug[i])]
        piv.append((i, c))
    z = [W.zero] * n
    for i, c in piv: z[c] = aug[i][n]
    return z

def monomials(deg, nvars=4):
    out = []
    def rec(prefix, remaining, k):
        if k == nvars - 1:
            out.append(tuple(prefix + [remaining])); return
        for e in range(remaining, -1, -1): rec(prefix + [e], remaining - e, k + 1)
    rec([], deg, 0)
    return out

def eval_monomial(W, xs, mono):
    """product of series xs[i]^mono[i]"""
    s = None
    for xi, e in zip(xs, mono):
        if e == 0: continue
        t = series_pow_list(W, xi, e)
        s = t if s is None else s * t
    return s if s is not None else LSer.const(W, W.one, n=xs[0].n)

def main(K, prec=220):
    W = Wk(K)
    st = json.load(open(f"{REL}/liftb-{K}.json"))
    x = [W.from_asc(v) for v in st["x"]]
    H, N, u0, q7, d8 = unpack(x, W)
    log(f"loaded lift at precision 13^{st['prec']}")
    # ---- adjoints and echelon basis at P_+ (as in petri13L)
    def pmodN(p): return pdivmod(W, p, N)[1]
    u0pow = [[W.one]]
    for _ in range(21): u0pow.append(pmodN(pmul(W, u0pow[-1], u0)))
    node_rows = [[W.zero] * 88 for _ in range(40)]
    for col, (i, j) in enumerate(MON):
        v = pmodN(pmul(W, [W.zero] * i + [W.one], u0pow[j]))
        for r in range(40): node_rows[r][col] = v[r] if r < len(v) else W.zero
    u, tau, kappa = branch_L(W, H, True, prec)
    dtau = LSer(W, -24, [W.scal(-23, kappa)], prec + 100)
    Hu = dH_series(W, H, tau, u, u.absprec() + 100)
    factor = dtau * Hu.inv()
    mons = monomial_series(W, tau, u, factor, prec)
    cusp_rows = [[s.coeff(e) for s in mons] for e in range(-88, 0)]
    ker = nullspace_W(W, node_rows + cusp_rows, 88)
    assert len(ker) == 4
    om_plus = [expand(W, v, mons) for v in ker]
    Jp, q, Ainvp, ech, pairs = petri_at(W, om_plus)
    qd = dict(zip(pairs, q))
    mu = W.mul(qd[(1, 2)], W.inv(qd[(2, 2)]))
    log(f"echelon basis at P_+ done; mu mod 13 = {W.to_fq(mu)}")
    # ---- omega-hat_j(w-hat) = mu^{-j} f_j(mu w-hat)   (f_j = ech_j as series in w~)
    muinv = W.inv(mu)
    def rescale(f, j):
        g = f.scale_var(mu)            # f(mu w)
        return g.scale(wpow(W, muinv, j))
    xh = [rescale(ech[j], j) for j in range(4)]           # omega-hat_j at P_+, series in w-hat (dw-hat implicit)
    for j in range(4):
        for e in range(4): assert xh[j].coeff(e) == (W.one if e == j else W.zero)
    # tau in w-hat: tau = kappa w~^-23 = kappa mu^-23 w-hat^-23
    kappa_hat = W.mul(kappa, wpow(W, muinv, 23))
    tau_hat = LSer(W, -23, [kappa_hat], prec + 300)
    # ---- Petri quadric in the hat coordinates (recompute from 12-jets)
    Jp2, qh, _, _, _ = petri_at(W, xh)      # petri_at re-echelonises (identity here) and normalises free column
    qhd = dict(zip(pairs, qh))
    q22 = qhd[(2, 2)]; qhat = {jk: W.mul(v, W.inv(q22)) for jk, v in qhd.items()}
    assert qhat[(1, 2)] == W.one and qhat[(2, 2)] == W.one and W.mul(qhat[(1, 3)], W.one) == Jp2 == Jp
    assert W.is_zero(qhat[(0, 0)]) and W.is_zero(qhat[(0, 1)]) and qhat[(0, 2)] == W.neg(qhat[(1, 1)]) and qhat[(0, 3)] == W.neg(W.one)
    log("Petri quadric in hat coordinates: q12 = q22 = 1, J_+ = q13")
    # ---- cubic: kernel of the 19-jet map on cubic monomials, reduced modulo x_i * Q
    cub = monomials(3); quad = monomials(2)
    cub_series = [eval_monomial(W, xh, m) for m in cub]
    rows = [[s.coeff(e) for s in cub_series] for e in range(19)]
    kc = nullspace_W(W, rows, 20)
    assert len(kc) == 5, len(kc)
    def vec_of_form(coeffs_by_mono, monos):
        return [coeffs_by_mono.get(m, W.zero) for m in monos]
    Qform = {}
    for (j, k), v in qhat.items():
        m = [0, 0, 0, 0]; m[j] += 1; m[k] += 1; Qform[tuple(m)] = W.add(Qform.get(tuple(m), W.zero), v)
    # x_i * Q as cubic vectors
    xiQ = []
    for i in range(4):
        d = {}
        for m, v in Qform.items():
            mm = list(m); mm[i] += 1; d[tuple(mm)] = W.add(d.get(tuple(mm), W.zero), v)
        xiQ.append(vec_of_form(d, cub))
    # find combination of kc with zero coefficients on x_2^2 x_i
    bad = [cub.index(tuple(sorted_tuple)) for sorted_tuple in [(1,0,2,0),(0,1,2,0),(0,0,3,0),(0,0,2,1)]]
    # solve: sum a_t kc[t] has zeros at bad positions: 4 conditions, 5 unknowns
    Msolve = [[kc[t][b] for t in range(5)] for b in bad]
    sol = nullspace_W(W, Msolve, 5)
    assert len(sol) == 1
    Cvec = [W.zero] * 20
    for t in range(5):
        if W.is_zero(sol[0][t]): continue
        Cvec = [W.add(a, W.mul(sol[0][t], b)) for a, b in zip(Cvec, kc[t])]
    norm_idx = cub.index((0, 2, 1, 0))       # x_1^2 x_2
    assert W.is_unit(Cvec[norm_idx]), "coefficient of x1^2 x2 not a unit"
    ci = W.inv(Cvec[norm_idx]); Cvec = [W.mul(ci, v) for v in Cvec]
    # verify C vanishes on the curve to order >= 19 (it is in the kernel)
    Cser = None
    for v, s in zip(Cvec, cub_series):
        if W.is_zero(v): continue
        t = s.scale(v); Cser = t if Cser is None else Cser + t
    assert Cser.val() is None or Cser.val() >= 19
    log("Petri cubic normalised (no x2^2-divisible monomials, coefficient of x1^2 x2 = 1)")
    # ---- quintics modulo the ideal
    quint = monomials(5)
    ideal_vecs = []
    for m in cub:
        d = {}
        for qm, v in Qform.items():
            mm = tuple(a + b for a, b in zip(m, qm)); d[mm] = W.add(d.get(mm, W.zero), v)
        ideal_vecs.append(vec_of_form(d, quint))
    Cform = {m: v for m, v in zip(cub, Cvec) if not W.is_zero(v)}
    for m in quad:
        d = {}
        for cm, v in Cform.items():
            mm = tuple(a + b for a, b in zip(m, cm)); d[mm] = W.add(d.get(mm, W.zero), v)
        ideal_vecs.append(vec_of_form(d, quint))
    # RREF of the ideal subspace (rows) with unit pivots
    basis = []
    for r in ideal_vecs:
        v = r[:]
        for pc, b in basis:
            if not W.is_zero(v[pc]):
                f = v[pc]; v = [W.sub(a, W.mul(f, c)) for a, c in zip(v, b)]
        pc = next((i for i, t in enumerate(v) if W.is_unit(t)), None)
        if pc is None:
            assert all(W.is_zero(t) for t in v); continue
        iv = W.inv(v[pc]); v = [W.mul(t, iv) for t in v]
        basis = [(qc, [W.sub(a, W.mul(b[pc], c)) for a, c in zip(b, v)] if not W.is_zero(b[pc]) else b) for qc, b in basis]
        basis.append((pc, v))
    assert len(basis) == 29, len(basis)
    pivots = sorted(pc for pc, _ in basis)
    nonpiv = [i for i in range(56) if i not in set(pivots)]
    assert len(nonpiv) == 27
    def normal_form(vec):
        v = vec[:]
        for pc, b in basis:
            if not W.is_zero(v[pc]):
                f = v[pc]; v = [W.sub(a, W.mul(f, c)) for a, c in zip(v, b)]
        assert all(W.is_zero(v[pc]) for pc in pivots)
        return v
    log(f"quintic normal forms: 27 standard monomials")
    # series of the 27 standard quintic monomials at P_+ and at P_-
    quint_series_plus = {i: eval_monomial(W, xh, quint[i]) for i in nonpiv}
    # ---- expansions at P_-
    u2, tau2, kappa2 = branch_L(W, H, False, prec)
    dtau2 = LSer(W, 22, [W.scal(23, kappa2)], prec + 100)
    Hu2 = dH_series(W, H, tau2, u2, u2.absprec() + 100)
    factor2 = dtau2 * Hu2.inv()
    mons2 = monomial_series(W, tau2, u2, factor2, prec)
    om_minus = [expand(W, v, mons2) for v in ker]
    # ech_j at P_- : same combinations as at P_+
    ech_minus = []
    for j in range(4):
        s = None
        for i in range(4):
            c = Ainvp[j][i]
            if W.is_zero(c): continue
            t = om_minus[i].scale(c); s = t if s is None else s + t
        ech_minus.append(s)
    xh_minus = [ech_minus[j].scale(wpow(W, muinv, j + 1)) for j in range(4)]     # omega-hat_j at P_- (series in w~')
    Pminus = [s.coeff(0) for s in xh_minus]
    log(f"P_- coordinates mod 13: {[str(W.to_fq(c)) for c in Pminus]}")
    quint_series_minus = {i: eval_monomial(W, xh_minus, quint[i]) for i in nonpiv}
    # ---- D-hat: quintic (normal form) vanishing to order >= 23 at P_+
    rowsD = [[quint_series_plus[i].coeff(e) for i in nonpiv] for e in range(23)]
    kD = nullspace_W(W, rowsD, 27)
    assert len(kD) == 4, len(kD)
    rowsN = [[quint_series_minus[i].coeff(e) for i in nonpiv] for e in range(23)]
    kN = nullspace_W(W, rowsN, 27)
    assert len(kN) == 4, len(kN)
    # ---- multiplication by tau is an isomorphism H^0(5K-23P_+) -> H^0(5K-23P_-), so the
    # pair (N,D) with N/D = tau is determined only up to the 4-dimensional choice of D.
    # Fix D-hat to be an RREF basis vector with ord_{P_+} exactly 23 and D(P_-) != 0, and
    # put N-hat := tau D-hat.  (This is the "selected multiplier row" of the char-13 model.)
    def ser(vec, sd):
        s = None
        for idx, val in zip(nonpiv, vec):
            if W.is_zero(val): continue
            t = sd[idx].scale(val); s = t if s is None else s + t
        return s
    cand = [(t, ser(kD[t], quint_series_plus), ser(kD[t], quint_series_minus)) for t in range(4)]
    log("ord_{P_+} on the H^0(5K-23P_+) basis: " + str([sp.val() for _, sp, _ in cand]))
    t0 = next(t for t, sp, sm in cand if sp.val() == 23 and sm.val() == 0)
    Dvec = kD[t0]; Dp = cand[t0][1]; Dm = cand[t0][2]
    Nser = [ser(kN[t], quint_series_plus) for t in range(4)]
    tD = tau_hat * Dp
    bcoef = solve_overdet(W, [[n.coeff(e) for n in Nser] for e in range(16)],
                          [tD.coeff(e) for e in range(16)], 4)
    for e in range(40):
        lhs = W.zero
        for t in range(4): lhs = W.add(lhs, W.mul(bcoef[t], Nser[t].coeff(e)))
        assert lhs == tD.coeff(e), f"tau*D != N at order {e}"
    Nvec = [W.zero] * 27
    for t in range(4):
        Nvec = [W.add(a, W.mul(bcoef[t], b)) for a, b in zip(Nvec, kN[t])]
    Np = ser(Nvec, quint_series_plus); Nm = ser(Nvec, quint_series_minus)
    assert (Dp.val(), Dm.val(), Np.val(), Nm.val()) == (23, 0, 0, 23), \
        (Dp.val(), Dm.val(), Np.val(), Nm.val())
    log("N-hat = tau D-hat to order 40 at P_+;  ord_{P_+}(D,N) = (23,0), ord_{P_-}(D,N) = (0,23)")
    assert W.is_unit(Pminus[0])
    pinv = W.inv(Pminus[0]); Pminus = [W.mul(pinv, v) for v in Pminus]
    out = {"K": K, "quad_pairs": pairs, "qhat": {f"{j},{k}": W.asc(v) for (j, k), v in qhat.items()},
           "cubic_monomials": cub, "Chat": [W.asc(v) for v in Cvec],
           "quintic_monomials": quint, "standard_quintics": nonpiv,
           "Nhat": [W.asc(v) for v in Nvec], "Dhat": [W.asc(v) for v in Dvec],
           "kappa_hat": W.asc(kappa_hat), "mu": W.asc(mu),
           "Pminus": [W.asc(v) for v in Pminus], "J_plus": W.asc(Jp),
           "kappa_minus_tilde": W.asc(kappa2)}
    json.dump(out, open(f"{REL}/canon-{K}.json", "w"))
    log(f"written canon-{K}.json")

if __name__ == "__main__":
    main(int(sys.argv[1]))
