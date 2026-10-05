# The eight-point certificate: (V9) and (V10) of Theorem 12.2

Exact, over L = F(sqrt(-23)), F = Q[y]/(y^6 - y^5 - 3y^4 - 3y^3 + y^2 + 5y + 4),
for the canonical model `relq/model-hat.gp` (SHA-256 `7423ca42...`).

**Claim certified.** The degree-23 function tau on X (Step 2 of Theorem 12.2) is
branched over exactly three points: 0, infinity, and one further value. Equivalently,
the eight simple ramification points of tau lie over a single value. Riemann–Hurwitz
then leaves no room for anything else, and with the branch divisor étale at 13 the
monodromy of the generic fibre equals that of the special fibre (Step 5).

**Why it is needed.** The 3 September draft derived the branch count from an
inequality |G_eta| <= |G_s| that is false in general (t = x^4 - 2x^2 + sx: D_4 at s = 0,
S_4 for s != 0). Reduction modulo 13 cannot see whether the eight critical values
merely coincide modulo 13. See Remark 12.5 of the paper.

## What runs

| stage | script | does | output | time (2 cores) |
|---|---|---|---|---|
| 1 | `stage1-ND.gp` (gp) | sections 1–5 of `relq/run_verify.gp`: the quintics D, N with N/D = tau and the leading coefficient c0 of 1/tau | `stage1-DN.bin`, `stage1.log` | 5 min |
| 2 | `stage2.py` (+ `lser.gp`) | P+ expansion to order 160; tau_1 = c0 tau; U_0 in L(5P+) with U_0(P-) = 0 via the residual divisor E_3 of the osculating plane; the plane relation H_0(tau_1, U_0) = 0 by back-substitution (72 coefficients, 44 gap checks, remainder zero to z^22) | `eightpoint-H0.txt`, `stage2-H0.bin`, `stage2.log` | 12 min |
| 3 | `stage3.py` (+ `stage3-lib.gp`) | b_3 mod 13 from the discriminant; Hensel lift of (b, q7, d8) to 13^2500; LLL recognition in L; **exact** checks H_0(b,U) = q7 d8^2, d8 squarefree of degree 8, gcd(q7,d8) = 1, gcd(d8, dH_0/dtau(b,U)) = 1, b a 13-unit, q7 squarefree, c0 a 13-unit | `eightpoint-witness.txt`, `stage3-b3.bin`, `stage3.log` | 150 s |

The proof is the exact arithmetic in stages 2 and 3 (the gap checks and the
vanishing of the remainder in stage 2; the eight `[1]`–`[8]` checks in stage 3; both scripts exit 1 on any failure). The
mod-13 reduction, the Hensel lift and the lattice reduction only *find* the witness.

`lser.gp` is a small exact Laurent-series library over L written for this test:
coefficients are kept integral with one denominator per series, over the primitive
element w = sqrt(-23) - beta (g = w^12 + 2w^11 + ... + 194457104), so that a series
product is 144 integer polynomial products plus the integral reduction modulo g.
It is checked against `relq/verify_L.gp` to order 20 inside `stage2.py`'s test
harness (`test_lser.py`).

## Replay

    pip install cypari2            # PARI 2.17.x; gp is needed for stage 1 only
    gp -q --default parisizemax=6000000000 stage1-ND.gp > stage1.log
    NPP2=160  python3 stage2.py
    NPREC=2500 python3 stage3.py   # NPREC=1500 already recognises everything (worst margin 13^17)

Expected last line of `stage3.log`:

    ==> EIGHT-POINT TEST: PASS -- the eight simple ramification points of tau lie over the single value b_3; tau has exactly three branch points

## Files

- `eightpoint-H0.txt` — H_0 as a 6 x 24 coefficient matrix over L (rows j = 0..5, columns
  k = 0..23; 49 nonzero entries, numerators/denominators up to 1588 digits), the quadric
  A_U, and U_0(P-).
- `eightpoint-witness.txt` — b_3 (396-digit coordinates over a 386-digit denominator) and
  the coefficient vectors of q7 (up to 606 digits) and d8 (up to 692 digits).
- `*.bin` — the same objects in PARI binary form (`read()`), used between stages.
- `stage*.log` — the runs on 2026-10-05 (cypari2 2.2.4 / PARI 2.17.2; gp 2.17.4).

Modulo 13 the discriminant of H_0 in U has degree 110 = 22 + 8 + 2·40: the point
b = 0 (multiplicity 22), the third branch value (multiplicity 8), and the 40 nodes of the
plane model (double roots consistent with 40 nodes) in Frobenius orbits 1 + 3 + 3 + 9 + 24.
