\\ ===================================================================
\\ lser.gp -- fast exact Laurent series over L = F(sqrt(-23)) = Q[w]/(g).
\\
\\ An element of L is Mod(p(w), g) ("Lw form"); the paper's nested form
\\ Mod(Mod(m(y),f) + Mod(n(y),f)*x, x^2+23) is converted with Lw() / Lnest().
\\ A Laurent series is S = [A, d, v, n]:  value = z^v * (sum_i A[i+1] w^i) / d,
\\ A[i] in Z[z] of degree < n (integral!), d in Z, known modulo z^(v+n).
\\ Multiplication is 144 integer polynomial products (Kronecker in PARI)
\\ plus the integral reduction of w^12..w^22 modulo the monic g.
\\ ===================================================================
LS_INF = 10^9;

ls_init(f0) =
{
  my(RE = rnfequation(f0, x^2+23, 1));
  g = subst(RE[1], x, w);
  LS_k = RE[3];
  YW = Mod(subst(lift(RE[2]), x, w), g);         \\ the root y of f in Q[w]/g
  XW = Mod(w, g) - LS_k*YW;                       \\ sqrt(-23):  w = x + k y
  if(subst(f0, y, YW) != 0, error("ls_init: y wrong"));
  if(XW^2 + 23 != 0, error("ls_init: x wrong"));
  \\ reduction rows: w^i = sum_m RED[i-11][m+1] w^m for i = 12..22
  RED = vector(11); my(cur = vector(12, i, -polcoef(g, i-1)));   \\ w^12
  RED[1] = cur;
  for(i = 13, 22, my(nxt = vector(12));
    \\ w * cur = sum cur[m+1] w^(m+1); w^12 -> RED[1]
    for(m = 0, 10, nxt[m+2] += cur[m+1]);
    nxt += cur[12] * RED[1];
    cur = nxt; RED[i-11] = cur);
  \\ the 12x12 matrix of (m_0..m_5, n_0..n_5) -> coefficients in w^0..w^11 of m(YW) + XW n(YW)
  my(M = matrix(12, 12));
  for(i = 0, 5, my(c = lift(YW^i)); for(r = 0, 11, M[r+1, i+1] = polcoef(c, r, w)));
  for(i = 0, 5, my(c = lift(XW*YW^i)); for(r = 0, 11, M[r+1, 7+i] = polcoef(c, r, w)));
  LS_M = M; LS_Minv = M^(-1);
  1;
}

\\ nested form -> Lw form
Lw(e) = my(l = lift(lift(e))); Mod(subst(subst(l, x, XW), y, YW), g);
\\ Lw form -> nested form
Lnest(e) =
{
  my(c = vector(12, r, polcoef(lift(e), r-1, w))~, s = LS_Minv * c);
  Mod(Mod(Pol(vector(6, i, s[7-i]), y), f) + Mod(Pol(vector(6, i, s[13-i]), y), f)*x, x^2+23);
}
\\ Lw element -> [integral 12-vector, denominator]
Lsplit(e) =
{
  my(l = lift(e), d = 1, v);
  v = vector(12, r, polcoef(l, r-1, w));
  d = denominator(content(v));
  [v * d, d];
}
Lvec(v) = Mod(Pol(vector(12, i, v[13-i]), w), g);     \\ integral 12-vector -> Lw element

\\ ---------- series constructors ----------
ls_const(e, v) = my(s = Lsplit(e)); [s[1], s[2], v, LS_INF];     \\ e * z^v, exact
ls_zero(n) = [vector(12), 1, 0, n];                                 \\ 0 + O(z^n)
ls_z() = [vector(12, i, if(i == 1, 1, 0)), 1, 1, LS_INF];           \\ the series z

ls_trunc(S, n) = if(n >= S[4], S, [vector(12, i, if(S[1][i] == 0, 0, S[1][i] % z^n)), S[2], S[3], n]);
ls_shift(S, a) = [S[1], S[2], S[3] + a, S[4]];
ls_neg(S) = [-S[1], S[2], S[3], S[4]];
ls_absprec(S) = if(S[4] >= LS_INF, LS_INF, S[3] + S[4]);

\\ divide out common content of A and d
ls_norm(S) =
{
  my(c = 0);
  for(i = 1, 12, if(S[1][i] != 0, c = gcd(c, content(S[1][i])); if(c == 1, break)));
  if(c == 0, return([S[1], 1, S[3], S[4]]));
  c = gcd(c, S[2]);
  if(c == 1, S, [S[1] / c, S[2] / c, S[3], S[4]]);
}

ls_mul(S1, S2) =
{
  my(n = min(S1[4], S2[4]), A1 = S1[1], A2 = S2[1], C = vector(23), A = vector(12), zn);
  for(j = 1, 12, if(A1[j] != 0, for(k = 1, 12, if(A2[k] != 0, C[j+k-1] += A1[j] * A2[k]))));
  for(i = 1, 12, A[i] = C[i]);
  for(i = 13, 23, if(C[i] != 0, A += C[i] * RED[i-12]));
  if(n < LS_INF, zn = z^n; for(i = 1, 12, if(A[i] != 0, A[i] = A[i] % zn)));
  ls_norm([A, S1[2]*S2[2], S1[3] + S2[3], n]);
}

ls_add(S1, S2) =
{
  my(v = min(S1[3], S2[3]), ab = min(ls_absprec(S1), ls_absprec(S2)), n, d, A1, A2, A, zn);
  n = if(ab >= LS_INF, LS_INF, ab - v);
  d = lcm(S1[2], S2[2]);
  A1 = S1[1] * (d / S1[2]); if(S1[3] > v, A1 = A1 * z^(S1[3] - v));
  A2 = S2[1] * (d / S2[2]); if(S2[3] > v, A2 = A2 * z^(S2[3] - v));
  A = A1 + A2;
  if(n < LS_INF, zn = z^n; for(i = 1, 12, if(A[i] != 0, A[i] = A[i] % zn)));
  ls_norm([A, d, v, n]);
}
ls_sub(S1, S2) = ls_add(S1, ls_neg(S2));
ls_scal(S, e) = ls_mul(S, ls_const(e, 0));          \\ e an Lw element

\\ the valuation (first nonzero coefficient), or [] if zero to precision
ls_val(S) =
{
  my(m = LS_INF);
  for(i = 1, 12, if(S[1][i] != 0, m = min(m, valuation(S[1][i], z))));
  if(m >= LS_INF, [], S[3] + m);
}
\\ coefficient of z^m as an Lw element
ls_coef(S, m) =
{
  if(m < S[3], return(Mod(0, g)));
  if(S[4] < LS_INF && m >= S[3] + S[4], error("ls_coef: beyond precision"));
  my(e = m - S[3]);
  Mod(Pol(vector(12, i, my(a = S[1][13-i]); if(a == 0, 0, polcoef(a, e, z))), w), g) / S[2];
}
\\ strip leading zero coefficients (exact zeros) so that S[3] is the true valuation
ls_strip(S) =
{
  my(v = ls_val(S));
  if(v == [], return(S));
  my(e = v - S[3]); if(e == 0, return(S));
  [vector(12, i, if(S[1][i] == 0, 0, S[1][i] \ z^e)), S[2], v, if(S[4] >= LS_INF, LS_INF, S[4] - e)];
}

\\ inverse of a series whose leading coefficient is a unit of L
ls_inv(S) =
{
  S = ls_strip(S);
  my(n = S[4], a0 = ls_coef(S, S[3]), b0 = 1/a0, X, m = 1, P = [S[1], S[2], 0, S[4]], two = ls_const(Mod(2, g), 0));
  if(n >= LS_INF, n = 2*max(1, vecmax(vector(12, i, if(S[1][i] == 0, 0, poldegree(S[1][i], z))))) + 1; P = ls_trunc(P, n); P[4] = n);
  X = ls_const(b0, 0); X[4] = 1;
  while(m < n,
    m = min(2*m, n);
    X = [X[1], X[2], 0, m];
    X = ls_mul(X, ls_sub(two, ls_mul(ls_trunc(P, m), X))));
  [X[1], X[2], -S[3], m];
}
ls_div(S1, S2) = ls_mul(S1, ls_inv(S2));
ls_pow(S, e) = my(R = ls_const(Mod(1, g), 0)); for(i = 1, e, R = ls_mul(R, S)); R;

\\ ---------- the curve: Q, C in Lw form ----------
ls_curve_init() =
{
  QUADw = matrix(4, 4, i, j, Lw(QUAD[i,j]));
  CUBCw = vector(#CUBC, t, Lw(CUBC[t]));
  PMw = vector(4, i, Lw(PM[i]));
  1;
}
\\ v = vector of 4 series (v[1] = 1): evaluate Q, C, and the partials wrt coordinates 3, 4
\\ returns [q, c, Qy, Qw, Cy, Cw] using a monomial table in (v3, v4)
ls_QC(v) =
{
  my(y1 = v[3], w1 = v[4], zz = ls_z(), one = ls_const(Mod(1, g), 0), tab = matrix(4, 4), q, c, Qy, Qw, Cy, Cw);
  \\ tab[b+1, c+1] = y1^b w1^c for b + c <= 3
  tab[1,1] = one; tab[2,1] = y1; tab[1,2] = w1;
  tab[3,1] = ls_mul(y1, y1); tab[2,2] = ls_mul(y1, w1); tab[1,3] = ls_mul(w1, w1);
  tab[4,1] = ls_mul(tab[3,1], y1); tab[3,2] = ls_mul(tab[3,1], w1); tab[2,3] = ls_mul(tab[1,3], y1); tab[1,4] = ls_mul(tab[1,3], w1);
  my(mono(a, b, c) = ls_shift(tab[b+1, c+1], a));     \\ z^a y^b w^c   (x0 = 1)
  \\ Q = sum_{i<=j} QUAD[i,j] x_i x_j
  q = ls_zero(LS_INF);
  for(i = 1, 4, for(j = i, 4, if(QUADw[i,j] != 0,
    my(m = vector(4)); m[i]++; m[j]++;
    q = ls_add(q, ls_scal(mono(m[2], m[3], m[4]), QUADw[i,j])))));
  Qy = ls_zero(LS_INF); Qw = ls_zero(LS_INF);
  for(i = 1, 4, if(QUADw[i,3] != 0, my(m = vector(4)); m[i]++; Qy = ls_add(Qy, ls_scal(mono(m[2],m[3],m[4]), QUADw[i,3])));
                if(QUADw[i,4] != 0, my(m = vector(4)); m[i]++; Qw = ls_add(Qw, ls_scal(mono(m[2],m[3],m[4]), QUADw[i,4]))));
  Qy = ls_add(Qy, ls_scal(y1, QUADw[3,3])); Qw = ls_add(Qw, ls_scal(w1, QUADw[4,4]));
  c = ls_zero(LS_INF); Cy = ls_zero(LS_INF); Cw = ls_zero(LS_INF);
  for(t = 1, #CUBM, if(CUBCw[t] != 0, my(m = CUBM[t]);
    c = ls_add(c, ls_scal(mono(m[2], m[3], m[4]), CUBCw[t]));
    if(m[3], Cy = ls_add(Cy, ls_scal(mono(m[2], m[3]-1, m[4]), m[3]*CUBCw[t])));
    if(m[4], Cw = ls_add(Cw, ls_scal(mono(m[2], m[3], m[4]-1), m[4]*CUBCw[t])))));
  [q, c, Qy, Qw, Cy, Cw];
}

\\ Newton expansion at P_+ = (1:0:0:0) with z = x1/x0 : returns [1, z, y(z), w(z)] to O(z^prec)
ls_localexp_plus(prec) =
{
  my(v = [ls_const(Mod(1,g),0), ls_z(), ls_zero(1), ls_zero(1)], n = 1);
  while(n < prec,
    n = min(2*n, prec);
    v[3] = [v[3][1], v[3][2], v[3][3], n]; v[4] = [v[4][1], v[4][2], v[4][3], n];
    my(E = ls_QC(v), q = E[1], c = E[2], a = E[3], b = E[4], cc = E[5], d = E[6], det, idet, du, dv);
    det = ls_sub(ls_mul(a, d), ls_mul(b, cc));
    idet = ls_inv(det);
    du = ls_mul(ls_add(ls_mul(ls_neg(q), d), ls_mul(c, b)), idet);
    dv = ls_mul(ls_add(ls_mul(ls_neg(a), c), ls_mul(cc, q)), idet);
    v[3] = ls_add(v[3], du); v[4] = ls_add(v[4], dv);
    v[3] = ls_trunc(v[3], n); v[4] = ls_trunc(v[4], n));
  v;
}
