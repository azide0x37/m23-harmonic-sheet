import cypari2, time, sys
pari = cypari2.Pari(); pari.allocatemem(4*10**9)
pari('T = varhigher("T"); w = varhigher("w")')
pari('read("model-hat.gp")')
pari('read("lser.gp")')
print("init:", pari('ls_init(f)'), pari('g'), "k=", pari('LS_k')); sys.stdout.flush()
pari('ls_curve_init()')
# round trip test
print("roundtrip:", pari('Lnest(Lw(QUAD[4,4])) == QUAD[4,4]'), pari('Lnest(Lw(PM[2])) == PM[2]'))
# series tests: (1+z)*(1-z) = 1 - z^2 ; inverse
pari('S1 = ls_add(ls_const(Mod(1,g),0), ls_z()); S2 = ls_sub(ls_const(Mod(1,g),0), ls_z())')
print("mul:", pari('ls_mul(S1,S2)'))
pari('S3 = ls_trunc(S1, 8); I3 = ls_inv(S3)'); print("inv(1+z):", pari('I3'))
print("check:", pari('ls_mul(S3, I3)'))
# compare the expansion with verify_L at order 20
t=time.time(); pari('vp2 = ls_localexp_plus(20)'); print("localexp 20:", time.time()-t); sys.stdout.flush()
pari('read("verify_L.gp"); default(seriesprecision, 100)')
pari('vp = localexp([L1, T + O(T^20), L0 + O(T^20), L0 + O(T^20)], [3,4], 20)')
ok = all(pari(f'Lnest(ls_coef(vp2[3], {k})) == polcoef(vp[3], {k})') and pari(f'Lnest(ls_coef(vp2[4], {k})) == polcoef(vp[4], {k})') for k in range(20))
print("agreement with verify_L to order 20:", ok)
for N in (40, 80, 160):
    t=time.time(); pari(f'vp2 = ls_localexp_plus({N})'); print(f"localexp {N}: {time.time()-t:.1f} s; digits of top coeff:", pari('#Str(vp2[3][1][1]) ')); sys.stdout.flush()
