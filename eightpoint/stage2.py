#!/usr/bin/env python3
"""Eight-point test, stage 2: the exact plane relation H_0(tau_1, U_0) = 0.

  tau_1 = c0 * tau,   1/tau_1 = z^23 (1 + O(z^3))          (z = x1/x0 at P_+)
  U_0 in L(5 P_+),    U_0 = z^-5 (1 + O(z)),  U_0(P_-) = 0   (U_0 = A_U / x3^2)
  H_0  = U^23 + sum h_{jk} tau^j U^k,  23 j + 5 k <= 115.

The pole orders 23j + 5k of the monomials at P_+ are distinct except for
U^23 and tau^5 (both 115), so H_0 is determined by back-substitution on the
Laurent expansion at P_+; the 44 orders m in [0,115] that are not of the
form 23j + 5k are consistency checks.  A solution with H_0(tau_1,U_0) = O(z)
is a global function without poles, hence 0: H_0(tau_1, U_0) = 0 in L(X).
Usage: NPP2=160 python3 stage2.py
"""
import cypari2, time, sys, os

NPP2 = int(os.environ.get("NPP2", "160"))
pari = cypari2.Pari()
pari.allocatemem(6 * 10**9)
LOG = open("stage2.log", "w")
def log(*a):
    s = " ".join(str(t) for t in a)
    print(s); sys.stdout.flush(); LOG.write(s + "\n"); LOG.flush()

t0 = time.time()
def el(): return f"[{time.time()-t0:7.1f} s]"

pari('T = varhigher("T"); w = varhigher("w")')
pari('read("model-hat.gp"); read("verify_L.gp"); read("lser.gp"); default(seriesprecision, 60)')
pari('ls_init(f); ls_curve_init()')
pari('[DVEC, NVEC, c0] = read("stage1-DN.bin")')
log(f"== stage 2 (fast library): NPP2 = {NPP2}; g = {pari('g')}")

# ---------- 2a. expansion at P_+ ----------
pari(f'vp2 = ls_localexp_plus({NPP2})')
pari('E = ls_QC(vp2)')
log(el(), "P_+ expansion done; Q, C vanish to precision:", pari('[ls_val(E[1]) == [], ls_val(E[2]) == []]'),
    "absprec:", pari('[ls_absprec(E[1]), ls_absprec(E[2])]'))
log(el(), "echelon: ord(y - z^2), ord(w - z^3) >= 4:",
    pari('[ls_val(ls_sub(vp2[3], ls_shift(ls_const(Mod(1,g),0), 2))), ls_val(ls_sub(vp2[4], ls_shift(ls_const(Mod(1,g),0), 3)))]'))

# ---------- 2b. D, N, tau_1 ----------
pari(r'''
YPW = matrix(6, 6);
YPW[1,1] = ls_const(Mod(1,g), 0);
for(b = 1, 5, YPW[b+1, 1] = ls_mul(YPW[b, 1], vp2[3]));
for(c = 1, 5, YPW[1, c+1] = ls_mul(YPW[1, c], vp2[4]));
for(b = 1, 4, for(c = 1, 5 - b, YPW[b+1, c+1] = ls_mul(YPW[b+1, c], vp2[4])));
''')
pari(r'''quintser3(cf) = {my(s = ls_zero(LS_INF)); for(k = 1, #cf, if(cf[k] != 0, my(m = QUINTM[STDQ[k]+1]); s = ls_add(s, ls_scal(ls_shift(YPW[m[3]+1, m[4]+1], m[2]), Lw(cf[k]))))); s}''')
pari('DP2 = ls_strip(quintser3(DVEC)); NPL2 = ls_strip(quintser3(NVEC))')
log(el(), "D, N series: valuations", pari('[DP2[3], NPL2[3]]'), "(need 23, 0)")
pari('INV1 = ls_strip(ls_scal(ls_div(DP2, NPL2), 1/Lw(c0)))')
log(el(), "1/tau_1: valuation", pari('INV1[3]'), " leading coefficient 1:", pari('ls_coef(INV1, 23) == 1'),
    " a1, a2 = ", pari('[ls_coef(INV1, 24), ls_coef(INV1, 25)]'), " relative precision", pari('INV1[4]'))
pari('TAU1 = ls_inv(INV1)')
log(el(), "tau_1: valuation", pari('TAU1[3]'), "relative precision", pari('TAU1[4]'))

# ---------- 2c. E_3 and U_0 (nested arithmetic for the small part) ----------
pari(r'''
s = varhigher("s");
q02 = QUAD[1,3]; q11 = QUAD[2,2];
PS = [-(q11 + s + s^2), q02*s, q02*s^2, L0];
c6 = Cev(PS);
e3 = c6 / s^3; e3 = e3 / polcoef(e3, 3, s);
pB = vector(4, i, Mod(PS[i], e3));
gQ = vector(4, k, Qd(pB, k)); gC = vector(4, k, Cd(pB, k));
vB = [gQ[3]*gC[4] - gQ[4]*gC[3], 0, gQ[4]*gC[1] - gQ[1]*gC[4], gQ[1]*gC[3] - gQ[3]*gC[1]];
QM = [[i,j] | i <- [1..4]; j <- [i..4]];
''')
pari(r'''condA(ij) = {my(i = ij[1], j = ij[2]); [pB[i]*pB[j], pB[i]*vB[j] + pB[j]*vB[i]]}''')
pari(r'''coords3(e) = {my(l = lift(e)); vector(3, k, polcoef(l, k-1, s))}''')
pari(r'''
MA = matrix(7, 10);
for(c = 1, 10, my(wv = condA(QM[c]), a = coords3(wv[1]), bb = coords3(wv[2])); MA[1, c] = if(QM[c] == [1,1], L1, L0); for(k = 1, 3, MA[1+k, c] = a[k]; MA[4+k, c] = bb[k]));
KA = matker(MA);
Qvec = vector(10, c, QUAD[QM[c][1], QM[c][2]]);
X3vec = vector(10, c, if(QM[c] == [4,4], L1, L0));
AU = 0; for(t = 1, #KA, if(matrank(Mat([Qvec~, X3vec~, KA[,t]])) == 3, AU = KA[,t]~; break));
''')
log(el(), "E_3: sextic degree", pari('poldegree(c6, s)'), "valuation at 0", pari('valuation(c6, s)'),
    "| e3 squarefree", pari('poldegree(gcd(e3, deriv(e3, s)), s) == 0'), "| e3(0) != 0", pari('subst(e3, s, 0) != 0'),
    "| conic param on Q", pari('Qev(PS) == 0'), "| x3(P_-) != 0", pari('PM[4] != L0'))
log(el(), "tangent v orthogonal to grad Q, grad C:", pari('[sum(k=1,4, gQ[k]*vB[k]) == 0, sum(k=1,4, gC[k]*vB[k]) == 0]'),
    "| v nonzero on E_3:", pari('poldegree(gcd(gcd(gcd(lift(vB[1]), lift(vB[3])), lift(vB[4])), e3), s) == 0'))
log(el(), "quadrics through P_+ + 2E_3: dim", pari('#KA'), "(need 3) | A_U found:", pari('AU != 0'))
pari(r'''
AUw = vector(10, c, Lw(AU[c]));
''')
pari(r'''QT(i, j) = if(i == 1 && j == 1, ls_const(Mod(1,g),0), if(i == 1, vp2[j], if(i == 2, ls_shift(vp2[j], 1), YPW[(i==3)+(j==3)+1, (i==4)+(j==4)+1])))''')
pari(r'''
AUser = ls_zero(LS_INF); for(c = 1, 10, if(AUw[c] != 0, AUser = ls_add(AUser, ls_scal(QT(QM[c][1], QM[c][2]), AUw[c]))));
AUser = ls_strip(AUser);
U0raw = ls_div(AUser, YPW[1,3]);
lead = ls_coef(U0raw, -5);
UPM = (sum(c = 1, 10, AUw[c] * PMw[QM[c][1]] * PMw[QM[c][2]]) / PMw[4]^2) / lead;
U0 = ls_sub(ls_scal(U0raw, 1/lead), ls_const(UPM, 0));
''')
log(el(), "ord_{P_+}(A_U) (need 1):", pari('AUser[3]'), "| U_0 valuation", pari('U0[3]'), "leading coeff 1:", pari('ls_coef(U0,-5) == 1'),
    "| rel. precision", pari('U0[4]'), "| digits of U_0(P_-):", pari('#Str(lift(UPM))'))

if NPP2 < 140:
    log("NPP2 < 140: test run, stopping after 2c"); sys.exit(0)

# ---------- 2d. powers ----------
pari(r'''
TPW = vector(6); TPW[1] = ls_const(Mod(1,g),0); for(j = 2, 6, TPW[j] = ls_mul(TPW[j-1], TAU1));
UPW = vector(24); UPW[1] = ls_const(Mod(1,g),0); for(k = 2, 24, UPW[k] = ls_mul(UPW[k-1], U0));
''')
log(el(), "powers done; absprec of tau^5, U^23:", pari('[ls_absprec(TPW[6]), ls_absprec(UPW[24])]'))

# ---------- 2e. back-substitution ----------
pari(r'''
S = UPW[24];
H = matrix(6, 24, i, j, 0);
checks = List(); ndet = 0;
{
  for(m = 0, 115, my(mm = 115 - m, found = 0);
    for(j = 0, 5, my(r = mm - 23*j); if(r >= 0 && r % 5 == 0, my(k = r/5); if(k <= 23 && !(j == 0 && k == 23),
      found = 1; ndet++;
      my(hc = -ls_coef(S, -mm));
      H[j+1, k+1] = hc;
      if(hc != 0, S = ls_add(S, ls_scal(ls_mul(TPW[j+1], UPW[k+1]), hc))))));
    if(!found, listput(checks, [mm, ls_coef(S, -mm) == 0])));
}
''')
ok_gaps = bool(pari('vecmin(vector(#checks, i, checks[i][2])) == 1'))
log(el(), "coefficients determined:", pari('ndet'), "| gap orders:", pari('#checks'), "| all gap checks vanish:", ok_gaps)
log("gap orders (m, coefficient of z^-m vanishes):", pari('Vec(checks)'))
ok_rem = bool(pari('ls_val(S) == []'))
log(el(), "remainder S = H_0(tau_1,U_0): absolute precision", pari('ls_absprec(S)'), "| zero to that precision:", ok_rem)
ok_h = bool(pari('H[6,1] == -1 && vecmin(vector(23, k, H[1,k] == 0)) == 1 && H[2,1] != 0'))
log("h_{5,0} = -1:", pari('H[6,1] == -1'), "| H_0(0,U) = U^23 (h_{0,k} = 0 for k < 23):", pari('vecmin(vector(23, k, H[1,k] == 0)) == 1'),
    "| Eisenstein h_{1,0} != 0:", pari('H[2,1] != 0'))
if not (ok_gaps and ok_rem and ok_h and int(pari('ndet')) == 72 and int(pari('#checks')) == 44):
    log("*** STAGE 2 FAILED"); sys.exit(1)
log("largest coefficient digits (Lw numerator):", pari('vecmax(vector(6, j, vecmax(vector(24, k, if(H[j,k] == 0, 0, #Str(lift(H[j,k])))))))'))
# nested form for stage 3
for fn in ("stage2-H0.bin", "eightpoint-H0.txt"):
    if os.path.exists(fn): os.remove(fn)
pari(r'''
Hn = matrix(6, 24, i, j, if(H[i,j] == 0, 0, Lnest(H[i,j])));
UPMn = Lnest(UPM);
writebin("stage2-H0.bin", [0, Hn, AU, UPMn, 0]);
''')
pari(r'write("eightpoint-H0.txt", "\\ H_0(tau, U) = U^23 + sum_{j,k} Hcoef[j+1,k+1] tau^j U^k (rows j = 0..5, columns k = 0..23), coefficients in L = F(sqrt(-23)) as Mod(Mod(m(y),f) + Mod(n(y),f)*x, x^2+23); AU = the quadric A_U on x_i x_j (i<=j) with U_0 = A_U/x_3^2 before normalisation; U0PM = the subtracted value U_0(P_-).")')
pari('write("eightpoint-H0.txt", "f = ", f); write("eightpoint-H0.txt", "Hcoef = ", Hn); write("eightpoint-H0.txt", "AU = ", AU); write("eightpoint-H0.txt", "U0PM = ", UPMn)')
log(el(), "written stage2-H0.bin (coefficient matrix, nested form) and eightpoint-H0.txt")
