\\ ===================================================================
\\ Exact verification over L = F(sqrt(-23)) of the descended Petri value U
\\ for the canonical model reconstructed from the 13-adic lift.
\\ Coordinates x_1..x_4 of GP are x_0..x_3 of the paper.
\\ ===================================================================

Qev(v) = my(s = L0); for(i=1,4, for(j=i,4, s += QUAD[i,j]*v[i]*v[j])); s;

Qd(v,k) = my(s = L0); for(i=1,4, s += QUAD[i,k]*v[i]); s + QUAD[k,k]*v[k];

Cev(v) = Pev(CUBM, CUBC, v);

Cd(v,k) =
{
  my(s = L0);
  for(t=1,#CUBM, my(m=CUBM[t]);
    if(m[k] && CUBC[t] != L0,
      my(p = CUBC[t]*m[k]);
      for(i=1,4, my(e = m[i] - (i==k)); if(e, p *= v[i]^e));
      s += p));
  s;
}

Pev(mons, coeffs, v) =
{
  my(s = L0);
  for(t=1,#mons,
    if(coeffs[t] != L0,
      my(m=mons[t], p=coeffs[t]);
      for(i=1,4, if(m[i], p *= v[i]^m[i]));
      s += p));
  s;
}

\\ Newton for the two solved coordinates as series in the free parameter T.
\\ base: the 4 coordinates with the free one set and the solved ones at their constant term.
localexp(base, sol, prec) =
{
  my(v = base, i = sol[1], j = sol[2], q, c, a, b, cc, d, det, du, dv, n = 1);
  while(n < prec,
    n = min(2*n, prec);
    v[i] = truncate(v[i]) + O(T^n); v[j] = truncate(v[j]) + O(T^n);
    q = Qev(v); c = Cev(v);
    a = Qd(v,i); b = Qd(v,j); cc = Cd(v,i); d = Cd(v,j);
    det = a*d - b*cc;
    du = (-q*d + c*b)/det; dv = (-a*c + cc*q)/det;
    v[i] = v[i] + du; v[j] = v[j] + dv;
  );
  v;
}

sermax(u, lo, hi) =
{
  my(s = 0, t = truncate(u));
  for(k = lo, hi, my(c = polcoef(t, k)); if(c != L0 && c != 0, s = max(s, #Str(lift(lift(c))))));
  s;
}

\\ ---------- series helpers ----------
trunc(u, n) = truncate(u + O(T^n));

\\ series of a quintic given by its coefficient vector on the standard monomials
quintser(cf, QS, n) =
{
  my(s = L0 + O(T^n));
  for(k = 1, #cf, if(cf[k] != 0, s += cf[k]*QS[k]));
  s;
}

\\ the 23rd root of u = c*T^(23m)*(1 + a1 T + ...): returns T^m*(1 + a1/23 T + ...)
root23(u, prec) =
{
  my(v = valuation(u,T), c, g, w, n = 1);
  if(v % 23, error("valuation not divisible by 23"));
  c = polcoef(truncate(u + O(T^(v+1))), v);
  g = u/(c*T^v) + O(T^prec);
  w = L1 + O(T^prec);
  while(n < prec,
    n = min(2*n, prec);
    w = trunc(w, n) + O(T^n);
    w = w - w*(w^23 - g)/(23*w^23));
  T^(v/23) * w;
}

\\ the canonical differential trivialisation eta = dT / jac, as a coefficient series
jacob(v) = Qd(v,3)*Cd(v,4) - Qd(v,4)*Cd(v,3);
