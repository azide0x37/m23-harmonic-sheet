\\ ===================================================================
\\ Eight-point test, stage 1: the exact quintics D, N with N/D = tau
\\ (identical to sections 1-5 of run_verify.gp; results saved with writebin).
\\ ===================================================================
T = varhigher("T");
read("model-hat.gp");
read("verify_L.gp");
default(seriesprecision, 200);
NPP = 42;
NPM = 34;

print("== stage 1: D, N (run_verify.gp sections 1-5)");
print("  model sha256 7423ca42... assumed; P_+ on X: ", [Qev([L1,L0,L0,L0]) == L0, Cev([L1,L0,L0,L0]) == L0]);
print("  P_- on X: ", [Qev(PM) == L0, Cev(PM) == L0]);
gettime();
vp = localexp([L1, T + O(T^NPP), L0 + O(T^NPP), L0 + O(T^NPP)], [3,4], NPP);
vm = localexp([L1, PM[2] + T + O(T^NPM), PM[3] + O(T^NPM), PM[4] + O(T^NPM)], [3,4], NPM);
print("  expansions: ", gettime(), " ms");
QSP = vector(#STDQ, k, my(m = QUINTM[STDQ[k]+1], p = L1 + O(T^NPP)); for(i=1,4, if(m[i], p *= vp[i]^m[i])); p);
QSM = vector(#STDQ, k, my(m = QUINTM[STDQ[k]+1], p = L1 + O(T^NPM)); for(i=1,4, if(m[i], p *= vm[i]^m[i])); p);
MP = matrix(23, #STDQ, r, c, polcoef(trunc(QSP[c], NPP), r-1));
MM = matrix(23, #STDQ, r, c, polcoef(trunc(QSM[c], NPM), r-1));
VP = matker(MP); VM = matker(MM);
print("  dim H^0(5K-23P_+), dim H^0(5K-23P_-) = ", [#VP, #VM], ": ", gettime(), " ms");
DS = vector(4, i, quintser(VP[,i]~, QSP, NPP));
NS = vector(4, i, quintser(VM[,i]~, QSP, NPP));
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
print("  multiplier solution space (need 1): ", #KM, ": ", gettime(), " ms");
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
print("  ord_{P_+}(D,N) = ", [valuation(DP,T), valuation(NPL,T)], "  ord_{P_-}(D,N) = ", [valuation(DM,T), valuation(NM,T)]);
INV = DP/NPL + O(T^(NPP-2));
c0 = polcoef(trunc(INV, 27), 23);
print("  1/tau = c0 z^23 (1 + a1 z + a2 z^2 + ...), [a1,a2] = ", [polcoef(trunc(INV,27),24)/c0, polcoef(trunc(INV,27),25)/c0]);
print("  digits of c0: ", #Str(lift(lift(c0))));
system("rm -f stage1-DN.bin");   \\ writebin appends to an existing file
writebin("stage1-DN.bin", [DVEC, NVEC, c0]);
print("  written stage1-DN.bin  [DVEC, NVEC, c0]");
