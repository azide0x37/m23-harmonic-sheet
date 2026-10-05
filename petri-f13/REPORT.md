# Characteristic-13 Petri connectedness route: bounded closeout

Date: 2026-09-02

## Current outcome after exact finite-field construction

The earlier numerical-recognition boundary below has been superseded by a
direct exact Riemann--Roch construction from the marked canonical model.
The following are now exact over

```text
k = F_13[a]/(a^6-a^5-3a^4-3a^3+a^2+5a+4) = F_(13^6).
```

* `exact_f13_map.py` constructs the two four-dimensional spaces
  `H0(5K-23P0)` and `H0(5K-23Pinf)` by exact jets.  The multiplication
  compatibility matrix has rank 15 in dimension 16, so its kernel is one
  dimensional and yields an invertible multiplier.  The resulting exact
  function has divisor `23P0-23Pinf`.
* `sage_passport_probe.py` verifies globally and exactly the passport
  `(23),(2^8 1^7),(23)`, separability, and Riemann--Hurwitz closure.
* `petri_binding_mod13.py` computes the marked Petri value by adapted formal
  parameters and obtains exactly

  ```text
  3 + 4a^2 + 2a^3 + 9a^4 + 11a^5.
  ```

  Its Frobenius degree is six and its minimal polynomial is the reduction of
  the certified Petri sextic.
* `sage_function_field_polynomial.py` eliminates the exact map to a unique
  monic, irreducible, separable degree-23 polynomial over `k(t)`.

The sole candidate-membership obstruction is now the exact upper monodromy
gate.  The branch data and primitive-group classification give the sharp
dichotomy

```text
G_geom = M23 or A23.
```

No artifact in this directory presently excludes `A23`.  Exhaustive bounded
Frobenius scans and timed Magma computations are corroboration only.

Magma constructs an exact compact `A23`-relative `M23` invariant at
`Worklevel=5`, cost 1012:

```text
J(x1,...,x23) = sum over 253 Witt heptads B of (sum_(i in B) xi)^5.
```

The machine receipt verifies the Steiner `S(4,7,23)` incidence property and
that the design automorphism group has order 10,200,960.  Applying this
separator still requires a certified global root labeling; unseeded search
would range over `[A23:M23]=1,267,136,462,592,000` conjugates.

The selected deterministic fallback is the direct unordered-four-set tower.
`MON13-FOURSET-ALGORITHM.md` constructs its degree-8855 field from the
quartic-divisor scheme and specifies the exact decisive factorization
`3+16` versus `19`.  This is not the ordered degree-212520 root tower.

## Historical outcome before the exact Riemann--Roch construction

The remainder of this report records the earlier bounded recognition stage
and its rejected lanes.  Statements below saying that no exact map existed
describe that earlier checkpoint and are retained for provenance.

The arithmetic endpoint is exact, but the actual-sheet binding is not yet
closed.

* The exact residue field is
  `k = F_13[alpha]/(alpha^6-alpha^5-3alpha^4-3alpha^3+alpha^2+5alpha+4)`;
  the modulus is irreducible, so `k = F_(13^6)`.
* At the prime `(13,s-4)`, `s^2=-23`, the conditional exact Petri quadric
  reduces integrally.  The chosen cubic gauge must first be multiplied by
  13; its four apparent denominator poles cancel only after substituting the
  13-adic lift of `s=4`.  Componentwise numerator/denominator reduction would
  be wrong.
* The reduced canonical candidate `C=V(Q,C3)` and all 10+20 coefficients are
  recorded exactly in `osculation-probe.json`.  The point `[1:0:0:0]` is on
  the curve and has nonzero Jacobian there.
* The exact descended Petri value is

      theta_bar = 3 + 4 alpha^2 + 2 alpha^3 + 9 alpha^4 + 11 alpha^5.

  Its Frobenius orbit has length six, and its minimal polynomial is exactly

      T^6 + 7 T^5 + T^4 + 6 T^3 + 7 T^2 + T + 1  over F_13.

  Thus `theta_bar` generates `F_(13^6)` and is a root of the reduced
  certified Petri sextic.
* No exact degree-23 map was emitted.  The numerical map is exceptionally
  stable, but its nontrivial coefficients do not pass the bounded exact
  recognition gate.  Consequently there is not yet an *actual* finite-field
  Hurwitz point to which `theta_bar` has been bound.

## Numerical map reconstruction

The upstream `Belyi/Code/triangle_phi.m` formula for the
`Delta(23,2,23)` Hauptmodul was ported exactly at the series level.  With

    A=1/4, B=27/92, C=24/23,
    A2=19/92, B2=1/4, C2=22/23,

write

    v = t (2F1(A,B;C;t) / 2F1(A2,B2;C2;t))^23

and revert to obtain `t=u(v)`.  The upstream canonical forms are evaluated
at `kappa*w` while `TrianglePhi` is evaluated at `w`, where

    kappa = 0.8330747178655918374843657780180158729...

is real.  Therefore the solver's original disc coordinate uses
`phi(w)=u((w/kappa)^23)`.

On the normalized canonical curve, 27 independent quintic monomials span
`H^0(5K)`.  Solving

    N(w) = phi(w) D(w),   N,D in H^0(5K),

gives matrix rank 50 and kernel dimension 4, exactly as
`h^0(5K-23P)=4` predicts.  At 100 p200 jets the maximum graph-kernel
residual is `2.20e-188`.

This also explains why the tempting `4K` route is invalid.  Exact reduction
at 13 gives rank 21 for the 23-jet map on `H^0(4K)`; its 14-dimensional
ambient kernel is precisely the degree-four canonical ideal, leaving quotient
kernel dimension zero.  Degree five is the first uniform Riemann--Roch gauge.

## Reproducible gauge and two-precision check

For the p200 label-1/label-2 pair, the representative was selected only from
denominator coordinates:

1. Column-normalize each four-dimensional graph kernel.
2. Among all four-row restrictions of the 27 denominator coordinates,
   maximize the minimum smallest singular value across both labels.
3. Set three selected coordinates to zero and one to one; choose the unit row
   minimizing the largest coefficient.
4. Choose at 100 jets, freeze the rows, and recompute independently at 120
   jets.

The selected section rows are `(13,16,19,20)`, with unit row 13 under the
correct raw-label scalar convention.  The maximum 100/120 drift is
`1.96e-186` for label 1 and `4.71e-186` for label 2.  This establishes a
precision-stable numerical graph representative, not an exact map.

## Semilinear convention audit

Two operations must not be conflated:

* scalar coefficient exactification uses
  `tau(c_l) = conjugate(c_kappa(l))`, with no block swap;
* target compatibility uses
  `tau(Graph(phi)) = Graph(1/phi)`, which swaps the `N,D` blocks and requires
  a separate target/source normalization comparison.

An earlier single-embedding 13-term relation and an earlier all-six
interpolation with a premature block swap were rejected and retained as
negative receipts.  They must not be used as candidates.

## Exact-recognition boundary

With the correct no-swap scalar rule, rows 0 and 1 of the p200 graph pass a
two-threshold `A(alpha),B(alpha)` recognition.  Row 2 fails at every tested
adjacent threshold from 190 down to 120.  Its successive relation heights
fall monotonically from about `9.5e53` to `1.5e34`, with no repeated exact
coefficient vector.  This is the signature of successive approximants, not a
certified relation.

The corrected all-six p100 interpolation is also decisive within its bound:

* a common 27-monomial basis was chosen with worst smallest singular value
  `3.86e-4` across all six embeddings;
* every graph matrix has rank 50 and residual between `1.29e-94` and
  `9.30e-87`;
* all interpolated power-basis coefficients are real to about `1e-110`, so
  the scalar Galois pairing is consistent;
* only seven forced zero/unit rows recognize.  Zero nontrivial rows pass the
  two-tolerance rational gate with coefficient bound `1e50`; 47 rows fail.

Therefore relaxing tolerances would manufacture candidates and is not
permitted.  If the observed p200 relation height near `1e54` is genuine, a
conservative seven-term integer-relation budget is roughly 300--400 reliable
decimal digits, versus about 185 currently available.  The clean next input
is either p350/p400 data (preferably all six sheets, so interpolation reduces
to rational recognition) or an exact algebraic solve.  More searches at the
current precision are not justified.

## Concrete exact construction and verification algorithm

Once an exact coefficient vector is available, the characteristic-13 route
is finite and deterministic:

1. **Integral reduction.** Express every coefficient of `Q,C3,N,D` as
   `A(alpha)+s B(alpha)`.  At the 13-adic prime `s=4`, combine `A+sB` before
   measuring valuation.  Multiply all coefficients in each projective
   equation/pencil by one common power of 13 clearing the minimum valuation;
   then reduce modulo `(13,s-4,g(alpha))`.  Reject if a denominator remains or
   either map half becomes zero.
2. **Genus-four M22 quotient.** Form
   `C=Proj k[x0,x1,x2,x3]/(Q,C3)`.  Verify the homogeneous Jacobian ideal has
   empty projective zero locus.  A smooth `(2,3)` complete intersection has
   genus 4.  This curve is the point-stabilizer quotient of the sought
   M23-Galois cover; the index of M22 in M23 is 23.
3. **Base divisor.** In the coordinate ring of `C`, compute the saturated
   common zero scheme `E=gcd_C(N,D)`.  It must have length 7 because a
   quintic canonical section has degree 30 and the residual map has degree
   23.  Verify no residual common base point remains.
4. **Two totally ramified fibers.** Primary-decompose/saturate the schemes
   `(Q,C3,N):E` and `(Q,C3,D):E`.  Each must have length 23 supported at one
   smooth point, with completed-local-ring order exactly 23.  These are the
   oriented 23A and 23B fibers.
5. **Involution fiber.** Decompose `(Q,C3,N-D):E`.  It must consist of eight
   length-2 local factors and seven reduced points, giving cycle type
   `2^8 1^7` (the degree-23 2A action).
6. **No hidden ramification.** The certified contributions are
   `22+22+8=52`, equal to
   `2g(C)-2+2*23 = 52`.  Since 13 is coprime to 2 and 23, the cover is tame;
   Riemann--Hurwitz then excludes any additional branch point.
7. **M23, not passport alone.** Eliminate to a degree-23 polynomial over
   `k(t)` and certify the geometric monodromy group, or bind the exact map to
   the already fixed Nielsen triple by a rigorous uniqueness/continuation
   certificate.  Cycle types alone do not identify M23.  Check that the
   point stabilizer has order `|M23|/23=443520` and is isomorphic to M22.
8. **Petri binding.** Evaluate the exact descended Petri function on this
   exact Hurwitz point and prove it equals
   `3+4alpha^2+2alpha^3+9alpha^4+11alpha^5`.  The existing exact receipt then
   proves its residue field is all of `F_(13^6)`.  Frobenius supplies one
   orbit of six nonharmonic sheets.

Steps 1--7 construct and certify the actual cover and its genus-four M22
quotient.  Step 8 is the theorem-critical horizontal moduli arrow.  At
present, Step 1 cannot start for `N,D`; only `Q,C3` and the arithmetic target
value are exact.

## Artifacts

* `numerical_map_probe.py`, `numerical-map-probe.json`: p200 TrianglePhi
  reconstruction and 4-dimensional `5K` graph kernel.
* `exact_map_candidate.py`: deterministic p200 gauge and bounded recognition
  program.  It exits nonzero at the first unstable coefficient.
* `p200-correct-scalar-recognition-failure.json`: authoritative p200 failure
  receipt for the corrected scalar convention.
* `p200-recognition-failure.json`: explicitly superseded wrong-swap negative.
* `p100_map_interpolation.py`, `p100-map-interpolation.json`: corrected
  all-six interpolation, common-basis conditioning, and 47-row failure
  receipt.
* `osculation_probe.py`, `osculation-probe.json`: exact reduced `Q,C3`, local
  marked point, and exact rejection of the `4K` shortcut.
* `petri_value_mod13.py`, `petri-value-mod13.json`: exact degree-six Petri
  value and Frobenius certificate.

## Historical fail-closed statement (superseded)

There is currently an exact genus-four *candidate model* over `F_(13^6)` and
an exact primitive Petri target value, plus a precision-stable numerical
degree-23 map.  There is not yet an exact `N/D`, exact passport, M23
monodromy certificate, or equality between an actual Hurwitz point's Petri
coordinate and that target.  Hence the Petri connectedness proposition is
not proved by these artifacts.
