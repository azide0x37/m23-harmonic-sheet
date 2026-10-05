\\ ===================================================================
\\ Eight-point test, stage 3 library (loaded by stage3.py).
\\ (i)  reduction of H_0 modulo the prime of L above 13 (y -> a in F_{13^6},
\\      sqrt(-23) -> 4), the discriminant in U, and b_3 mod 13;
\\ (ii) Hensel lift of (b, q7, d8) in Z_13[y]/(f) with H_0(b,U) = q7 d8^2;
\\ (iii) LLL recognition in L;  (iv) the exact checks are in stage3.py.
\\ ===================================================================

\\ ---------- (i) reduction mod 13 ----------
dens(l) = if(type(l) == "t_POL", lcm(vector(poldegree(l)+1, i, dens(polcoef(l,i-1)))), denominator(l));
red13(e) =
{
  my(l = lift(lift(e)), d, num);
  if(l == 0, return(0*ffa));
  d = dens(l); num = l*d;
  if(d % 13 == 0, error("denominator divisible by 13"));
  subst(subst(num, x, 4), y, ffa) / d;
}
mod13_setup() =
{
  Hbar = sum(k = 0, 23, sum(j = 0, 5, red13(polcoef(polcoef(H0, k, U), j, b)) * b^j * U^k));
  Dbar = polresultant(Hbar, deriv(Hbar, U), U);
  FD = factor(Dbar);
  cand = [];
  for(i = 1, #FD~, if(poldegree(FD[i,1], b) == 1 && FD[i,2] >= 8 && polcoef(FD[i,1], 0, b) != 0,
    cand = concat(cand, [[-polcoef(FD[i,1], 0, b)/polcoef(FD[i,1], 1, b), FD[i,2]]])));
  if(#cand != 1, error("no unique root of multiplicity >= 8 mod 13"));
  b3bar = cand[1][1];
  FH = factor(subst(Hbar, b, b3bar));
  d8bar = 1; q7bar = 1;
  for(i = 1, #FH~, if(FH[i,2] == 2, d8bar *= FH[i,1], if(FH[i,2] == 1, q7bar *= FH[i,1], error("multiplicity > 2 over b_3 mod 13"))));
  [poldegree(Dbar, b), valuation(Dbar, b), cand[1][2], vector(#FD~, i, [poldegree(FD[i,1], b), FD[i,2]]),
   poldegree(q7bar, U), poldegree(d8bar, U), vector(#FH~, i, [poldegree(FH[i,1], U), FH[i,2]])];
}

\\ ---------- (ii) Hensel lift in W = Z_13[y]/(f) ----------
toW(e) = Mod(subst(lift(lift(e)) * (1 + O(13^N)), x, s13), f);
ffW(a) = Mod(subst(a.pol * 1, variable(ffa.pol), y) * (1 + O(13^N)), f);
Wmod13(wv) = subst(lift(wv) * Mod(1,13), y, ffa);
hensel_setup() =
{
  s13 = sqrt(-23 + O(13^N)); if(lift(s13) % 13 != 4, s13 = -s13);
  H0W = sum(k = 0, 23, sum(j = 0, 5, toW(polcoef(polcoef(H0, k, U), j, b)) * b^j * U^k));
  dH0W = deriv(H0W, b);
  X = vector(16); X[1] = ffW(b3bar);
  for(i = 0, 6, X[2+i] = ffW(polcoef(q7bar, i, U)));
  for(i = 0, 7, X[9+i] = ffW(polcoef(d8bar, i, U)));
  1;
}
q7of(X) = U^7 + sum(i = 0, 6, X[2+i]*U^i);
d8of(X) = U^8 + sum(i = 0, 7, X[9+i]*U^i);
Eof(X) = my(e = subst(H0W, b, X[1]) - q7of(X)*d8of(X)^2); vector(23, k, polcoef(e, k-1, U));
Jof(X) =
{
  my(q = q7of(X), d = d8of(X), dd = d^2, qd = 2*q*d, M = matrix(23, 16), col);
  col = subst(dH0W, b, X[1]); for(k = 1, 23, M[k,1] = polcoef(col, k-1, U));
  for(i = 0, 6, my(p = -U^i*dd); for(k = 1, 23, M[k, 2+i] = polcoef(p, k-1, U)));
  for(i = 0, 7, my(p = -U^i*qd); for(k = 1, 23, M[k, 9+i] = polcoef(p, k-1, U)));
  M;
}
Evalv(E) =
{
  my(m = N);
  for(k = 1, #E, my(l = lift(E[k])); if(l != 0,
    for(i = 0, poldegree(l), my(c = polcoef(l, i)); if(c != 0, m = min(m, valuation(c, 13))))));
  m;
}
hensel_rows() =
{
  my(Jbar = matrix(23, 16, i, j, Wmod13(Jof(X)[i,j])));
  rows = matindexrank(Jbar)[1];
  [matrank(Jbar), rows];
}
hensel_step() =
{
  my(E = Eof(X), v = Evalv(E));
  if(v >= N, return(v));
  my(J = Jof(X), Jr = matrix(16, 16, i, j, J[rows[i], j]), Er = vector(16, i, E[rows[i]])~);
  my(delta = matsolve(Jr, Er));
  X = vector(16, i, X[i] - delta[i]);
  v;
}

\\ ---------- (iii) LLL recognition ----------
recog(wv) =
{
  my(l = lift(wv), beta = vector(6, i, lift(polcoef(l, i-1))), pN = 13^N, sl = lift(s13), M, red, v, d, m, n, elt, nrm, h);
  M = matrix(19, 19);
  for(i = 1, 13, M[i, i] = 1);
  for(i = 1, 6, M[13 + i, 1] = pN * (beta[i] % pN));
  for(i = 1, 6, M[13 + i, 1 + i] = -pN);
  for(i = 1, 6, M[13 + i, 7 + i] = -pN * (sl % pN));
  for(i = 1, 6, M[13 + i, 13 + i] = pN * pN);
  red = M * qflll(M);
  v = red[, 1];
  if(vecmax(abs(vector(6, i, v[13+i]))) != 0, error("LLL: shortest vector violates the congruence"));
  d = v[1]; m = vector(6, i, v[1+i]); n = vector(6, i, v[7+i]);
  if(d < 0, d = -d; m = -m; n = -n);
  if(d == 0, error("LLL: zero denominator"));
  elt = Mod(Mod(Pol(vector(6, i, m[7-i]), y), f) + Mod(Pol(vector(6, i, n[7-i]), y), f)*x, x^2+23) / d;
  nrm = d^2 + sum(i=1,6, m[i]^2 + n[i]^2);
  h = log(nrm)/(2*log(13));
  [elt, h, 6*N/13 - h];
}
