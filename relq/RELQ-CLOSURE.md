> **Update 2026-10-05.** Item 4 of the referee pass at the end of this note is
> withdrawn: the inequality |G_eta| <= |G_s| it relies on compares G_eta with the
> stabilizer of a special component of the Galois closure (a group containing
> inertia), not with the monodromy of the special fibre, and it fails for
> t = x^4 - 2x^2 + sx. The branch count is now an exact certificate, (V9)-(V10),
> in `../eightpoint/` (plane relation H_0(tau_1,U_0) = 0 of bidegree (5,23) and the
> exact factorization H_0(b,U) = q_7 d_8^2); monodromy then follows by Grothendieck's
> tame specialization theorem. See Theorem 12.2, Steps 4-5, and Remark 12.5 of the paper.

# RELQ closed: g_U = Phi

**Date** 2026-09-03. **Claim** The characteristic polynomial of the descended Petri
coordinate U on the six-dimensional factor E of H^0(H, O) equals the recognized sextic Phi;
hence E = Q[x]/(x^6 - x^5 - 3x^4 - 3x^3 + x^2 + 5x + 4).

## Result

    minpoly_Q(U) = 50019668646075143 x^6 + 206305139624557234 x^5 + 353749018385880145 x^4
                 -  55027953306246948 x^3 - 759255603013669319 x^2 - 868575878997768238 x
                 - 139453053861313073
                 = 50019668646075143 * Phi
    polredabs(minpoly_Q(U)) = x^6 - x^5 - 3x^4 - 3x^3 + x^2 + 5x + 4 = f
    U = r(beta), the embedding Q[x]/(Phi) -> F recorded in the receipts appendix.

## Route

Exact model, not height bounds and not more primes.

1. `lift13b.py 400` -- the singularity formulation of the plane model H(tau,U) (143 unknowns:
   48 coefficients, the degree-40 node polynomial N, u0 mod N, q7, d8) has an invertible
   Jacobian mod 13, so the characteristic-13 point lifts uniquely to W_400 = W(F_13^6)/13^400.
2. `canon13.py 400` -- the canonical model in the intrinsic echelon coordinates at P_+:
   the Petri quadric (normalised q12 = q22 = 1), the Petri cubic, the second wild point P_-.
3. `recon_canon.py 400` -- lattice reduction over L = F(sqrt(-23)) recognises those:
   quadric 53 digits (margin 13^137), cubic 122 digits (margin 13^75), P_- 53 digits
   (margin 13^137), against an expected shortest vector of 13^(6*400/13) = 13^184.6.
   The quintics N, D are NOT recognisable and were abandoned: exactly, their coefficients
   have 7823 and 13742 decimal digits.
4. `run_verify.gp` (PARI, 150 s) -- everything else is computed exactly over L from the
   recognised model:
   (V1) P_+, P_- on X; the Petri normalisation; 13-integrality.
   (V2) rank 21 = h^0(4K) for the 23-jet map on quartics at both points, so
        h^0(4K - 23P) = 0 and |5K - 23P| is base-point-free.
   (V3) dim H^0(5K - 23P_+) = dim H^0(5K - 23P_-) = 4; the space of M with
        M(D)D' = M(D')D on Taylor coefficients 23..38 at P_+ is exactly 1-dimensional.
   (V4) ord_{P_+}(D, M(D)) = (23, 0), ord_{P_-}(D, M(D)) = (0, 23).
   (V5) the model is the echelon basis for z = x1/x0: x2/x0 - z^2, x3/x0 - z^3 vanish to
        order 4 and omega_0/dz is constant to order 3.
   (V6) 1/tau = c z^23 (1 + 0*z + 0*z^2 + ...), so z is an admissible parameter modulo
        O(z^4) and J_+ = q13/q22 is read off the quadric.
   (V7) at P_- the adapted parameter is tau^(1/23); echelonising and transforming the
        quadric gives J_-, the sqrt(-23)-conjugate of J_+.
   (V8) U = J_+ + J_- + J_+ J_- lies in F and has minimal polynomial Phi.

The geometry that makes this a point of H -- deg tau = 23 (from (V2)+(V3)), three branch
points, monodromy M23, the inertia labels -- is Steps 2-6 of Theorem 12.2 in the paper:
base-point-freeness over L, then specialisation to the characteristic-13 fibre, then
"M23 contains no transposition" to force the eight simple ramification points over a single
third branch value.

No height bound. No second prime. 13 is used only to certify that the special fibre is
smooth and to supply its monodromy and inertia labels.

## Discipline

Recognition proposed Q-hat, C-hat, P_-. Binding is (V1)-(V8) plus Steps 1-8: exact
computations with the printed data and statements about the object those data define.
How the model was found plays no role in the proof.

## Files

    model-hat.gp             the exact model over L (sha256 7423ca42...)
    run_verify.gp            (V1)-(V8), PARI
    verify_L.gp              series and evaluation helpers
    relq-verification.log    the run
    lift13b.py               13-adic lift (singularity formulation)
    canon13.py               canonical model in echelon coordinates over W_K
    recon_canon.py           lattice reduction over L
    emit_L.py                emits model-hat.gp

## Referee pass (same day)

Four things were found and fixed in Theorem 12.2.

1. **Flatness of the integral model was assumed.** Now argued: the model is a complete
   intersection in a regular scheme, hence Cohen-Macaulay and unmixed; every component
   has dimension 2 while the special fibre has dimension 1, so no associated point lies in
   the special fibre.
2. **tau was treated as a morphism without justification.** N and D share the residual
   divisor E, so [N : D] is not defined there. Now argued from the divisor: after scaling
   the two quintics primitive, div(tau) on the total space is the difference of the two
   DISJOINT sections 23 P_+ and 23 P_-, so by normality every point of the model misses one
   of them and tau is regular or has a pole there; proper with finite fibres gives finite,
   and finite between regular schemes of equal dimension gives flat.
3. **The reduction of tau is the INVERSE of the char-13 map, not a multiple of it**
   (tau has a pole where t has a zero); they differ by an automorphism of the target and
   define the same cover.
4. **The monodromy step was the real gap.** It rested on Zariski connectedness of the
   Galois closure, which needs geometric connectedness the argument did not supply; and the
   obvious repair -- SGA 1 Exp. XIII specialization -- does NOT apply, because every
   statement of that theorem requires the branch divisor to be ETALE over the base, which
   is exactly what cannot be assumed while the eight simple branch points might still
   collide in the special fibre. Replaced by an elementary argument:

   Let W be the Galois closure of the cover over the complement of the branch locus, G its
   group, and W-hat the normalisation of P^1_O in W.  The geometric fibres of W-hat have
   [G : G_eta] and [G : G_s] components; by Stein factorization both counts are counts of
   points of one finite flat O-algebra A, the first equal to the rank of A and the second at
   most the rank.  So |G_eta| <= |G_s| = |M23|.  And G_eta is transitive of degree 23 and
   contains an involution with at most eight transpositions.  Of the seven transitive groups
   of degree 23, A_23 and S_23 are too large, C_23 and C_23:C_11 have odd order, and the
   involutions of D_23 and AGL(1,23) fix one point and so have eleven transpositions.  Hence
   G_eta = M23, and then "every involution of M23 has cycle type 2^8 1^7" forces the eight
   simple ramification points over a single branch value.

References added: Stacks Project tags 0BWA and 0BW7 (the different is an effective Cartier
divisor whose formation commutes with base change), Serre Local Fields Ch. III Sec. 6 (tame
different = e - 1), Mueller Arch. Math. 85 (2005) (Burnside: transitive of prime degree is
solvable or 2-transitive), Guralnick J. Algebra 81 (1983) (the 2-transitive groups of degree
23 are M23, A23, S23), Hulpke J. Symbolic Comput. 39 (2005) (the tables), ATLAS (M23 has one
class of involutions, type 2^8 1^7).
