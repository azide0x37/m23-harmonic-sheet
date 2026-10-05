\\ ===================================================================
\\ RELQ, exact-model route:  compute the descended Petri coordinate U of the
\\ reconstructed M23 cover over L = F(sqrt(-23)) and compare with Phi.
\\ ===================================================================
T = varhigher("T");
read("model-hat.gp");
read("verify_L.gp");
default(seriesprecision, 200);
NPP = 42;      \\ series order at P_+
NPM = 34;      \\ series order at P_-

print("== 0. the model");
print("  P_+ = (1:0:0:0) on X: ", [Qev([L1,L0,L0,L0]) == L0, Cev([L1,L0,L0,L0]) == L0]);
print("  P_- on X: ", [Qev(PM) == L0, Cev(PM) == L0]);
print("  q00=q01=0, q03=-1, q02=-q11, q12=q22=1: ", [QUAD[1,1] == L0, QUAD[1,2] == L0, QUAD[1,4] == -L1, QUAD[1,3] == -QUAD[2,2], QUAD[2,3] == L1, QUAD[3,3] == L1]);

print("== 1. local expansions");
gettime();
vp = localexp([L1, T + O(T^NPP), L0 + O(T^NPP), L0 + O(T^NPP)], [3,4], NPP);
print("  P_+ to order ", NPP, ": ", gettime(), " ms");
vm = localexp([L1, PM[2] + T + O(T^NPM), PM[3] + O(T^NPM), PM[4] + O(T^NPM)], [3,4], NPM);
print("  P_- to order ", NPM, ": ", gettime(), " ms");
print("  Q,C vanish at P_+: ", [valuation(Qev(vp),T), valuation(Cev(vp),T)]);
print("  Q,C vanish at P_-: ", [valuation(Qev(vm),T), valuation(Cev(vm),T)]);

print("== 2. the model is echelon at P_+ for the parameter z = x1/x0");
print("  ord(x2/x0 - z^2) (need >= 4): ", valuation(vp[3] - T^2 + O(T^NPP), T));
print("  ord(x3/x0 - z^3) (need >= 4): ", valuation(vp[4] - T^3 + O(T^NPP), T));
ETAP = 1/jacob(vp);
ETAP = ETAP/polcoef(trunc(ETAP, NPP), 0);
print("  ord(omega_0/(c dz) - 1) (need >= 3): ", valuation(ETAP - 1 + O(T^NPP), T));

print("== 3. the two quintic systems");
QSP = vector(#STDQ, k, my(m = QUINTM[STDQ[k]+1], p = L1 + O(T^NPP)); for(i=1,4, if(m[i], p *= vp[i]^m[i])); p);
QSM = vector(#STDQ, k, my(m = QUINTM[STDQ[k]+1], p = L1 + O(T^NPM)); for(i=1,4, if(m[i], p *= vm[i]^m[i])); p);
print("  quintic series: ", gettime(), " ms");
print("  rank of the 27 standard quintics (need 27): ", matrank(matrix(31, #STDQ, r, c, polcoef(trunc(QSP[c], NPP), r-1))));
MP = matrix(23, #STDQ, r, c, polcoef(trunc(QSP[c], NPP), r-1));
MM = matrix(23, #STDQ, r, c, polcoef(trunc(QSM[c], NPM), r-1));
VP = matker(MP); VM = matker(MM);
print("  dim H^0(5K-23P_+), dim H^0(5K-23P_-) = ", [#VP, #VM], ": ", gettime(), " ms");
DS = vector(4, i, quintser(VP[,i]~, QSP, NPP));
NS = vector(4, i, quintser(VM[,i]~, QSP, NPP));
print("  ord_{P_+} of the V_+ basis: ", vector(4,i,valuation(DS[i],T)));
print("  ord_{P_+} of the V_- basis: ", vector(4,i,valuation(NS[i],T)));

print("== 3b. base-point-freeness of |5K-23P_+| and |5K-23P_-|");
QUARTM = [[a,b,c,4-a-b-c] | a <- [0..4]; b <- [0..4-a]; c <- [0..4-a-b]];
{
  QP4 = vector(#QUARTM, k, my(m = QUARTM[k], p = L1 + O(T^NPP)); for(i=1,4, if(m[i], p *= vp[i]^m[i])); p);
  QM4 = vector(#QUARTM, k, my(m = QUARTM[k], p = L1 + O(T^NPM)); for(i=1,4, if(m[i], p *= vm[i]^m[i])); p);
}
print("  quartic monomials: ", #QUARTM);
print("  rank of the 23-jet map on quartics at P_+ (need 21 = h^0(4K)): ", matrank(matrix(23, #QUARTM, r, c, polcoef(trunc(QP4[c], NPP), r-1))));
print("  rank of the 23-jet map on quartics at P_- (need 21 = h^0(4K)): ", matrank(matrix(23, #QUARTM, r, c, polcoef(trunc(QM4[c], NPM), r-1))));
print("  ==> h^0(4K-23P) = 0, so |5K-23P| is base-point-free and any g with g H^0(5K-23P_+) = H^0(5K-23P_-) has divisor 23P_- - 23P_+");

print("== 4. the multiplier: the unique g with g H^0(5K-23P_+) = H^0(5K-23P_-)");
{
  PR = matrix(4,4);
  for(i=1,4, for(j=1,4, PR[i,j] = trunc(NS[j]*DS[i], 40)));
  EQ = matrix(96, 16); row = 0;
  for(i=1,3, for(k=i+1,4, for(e=23,38, row++;
    for(j=1,4,
      EQ[row, 4*(i-1)+j] = EQ[row, 4*(i-1)+j] + polcoef(PR[k,j], e);
      EQ[row, 4*(k-1)+j] = EQ[row, 4*(k-1)+j] - polcoef(PR[i,j], e)))));
}
KM = matker(EQ);
print("  solution space (need 1): ", #KM, ": ", gettime(), " ms");
MU = matrix(4,4,i,j, KM[4*(i-1)+j, 1]);
{
  NHAT = vector(4, i, my(v = vector(#STDQ, k, L0));
    for(j=1,4, if(MU[i,j] != 0, for(k=1,#STDQ, v[k] = v[k] + MU[i,j]*VM[k,j]))); v);
  idx = 0;
  for(i=1,4, if(idx == 0 && valuation(DS[i],T) == 23, idx = i));
}
DVEC = VP[,idx]~; NVEC = NHAT[idx];
DP = quintser(DVEC, QSP, NPP); NPL = quintser(NVEC, QSP, NPP);
DM = quintser(DVEC, QSM, NPM); NM = quintser(NVEC, QSM, NPM);
print("  index ", idx, ": ord_{P_+}(D,N) = ", [valuation(DP,T), valuation(NPL,T)]);
print("           ord_{P_-}(D,N) = ", [valuation(DM,T), valuation(NM,T)]);

print("== 5. tau at P_+ and admissibility of z");
INV = DP/NPL + O(T^(NPP-2));
c0 = polcoef(trunc(INV, 27), 23);
print("  ord_{P_+}(1/tau) = ", valuation(INV,T));
print("  a1, a2 in 1/tau = c z^23 (1 + a1 z + a2 z^2 + ...): ", [polcoef(trunc(INV,27),24)/c0, polcoef(trunc(INV,27),25)/c0]);
JPLUS = QUAD[2,4]/QUAD[3,3];
print("  J_+ = q13/q22 matches the LLL value: ", JPLUS == JPLUS_recognised);

print("== 6. the echelon basis at P_- and J_-");
TAUM = NM/DM;
print("  ord_{P_-} tau = ", valuation(TAUM,T));
WM = root23(TAUM, 7);
TW = serreverse(WM + O(T^7));
JMSER = jacob(vm);
DTDW = deriv(TW, T);
OM = vector(4, j, subst(trunc(vm[j]/JMSER, 7), T, TW) * DTDW + O(T^6));
A4 = matrix(4,4,k,i, polcoef(trunc(OM[i],6), k-1));
B = matsolve(A4, matid(4))~;
{
  for(k=1,4, my(o = L0 + O(T^6)); for(i=1,4, o = o + B[k,i]*OM[i]);
     print("  echelon row ", k, ": ord(omega'_k - w'^", k-1, ") = ", valuation(o - T^(k-1) + O(T^6), T)));
}
AS = matrix(4,4,i,j, if(i==j, QUAD[i,i], QUAD[i,j]/2));
Binv = matsolve(B, matid(4));
ASP = Binv~ * AS * Binv;
JMINUS = 2*ASP[2,4]/ASP[3,3];

print("== 7. U and its minimal polynomial");
UU = JPLUS + JMINUS + JPLUS*JMINUS;
print("  U has no sqrt(-23) part: ", polcoef(lift(UU), 1) == 0);
uF = polcoef(lift(UU), 0);
MPU = minpoly(uF);
PHI = x^6 + (206305139624557234*x^5 + 353749018385880145*x^4 - 55027953306246948*x^3 - 759255603013669319*x^2 - 868575878997768238*x - 139453053861313073)/50019668646075143;
print("  minpoly(U) = ", MPU);
print("  ==> g_U = Phi : ", MPU/polcoef(MPU,6) == PHI);
print("  polredabs(minpoly(U)) = ", polredabs(MPU));
print("  J_- = ", JMINUS);
print("  U   = ", uF);
print("  digits in the coefficients of D-hat, N-hat: ", [vecmax(vector(#STDQ,k,#Str(lift(lift(DVEC[k]))))), vecmax(vector(#STDQ,k,#Str(lift(lift(NVEC[k])))))]);
write("relq-receipt.txt", "minpoly(U) = ", MPU);
write("relq-receipt.txt", "polredabs   = ", polredabs(MPU));
write("relq-receipt.txt", "J_plus = ", JPLUS);
write("relq-receipt.txt", "J_minus = ", JMINUS);
write("relq-receipt.txt", "U = ", uF);
