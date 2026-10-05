# MON13: direct unordered four-set tower

## Status

This is a proof-grade construction plan, not a completed monodromy
certificate.  It constructs the degree-`8855` fixed field directly and never
constructs the degree-`212520` ordered four-root tower.

An exact implementation checkpoint is now available at the split prime
`p=9319`.  It proves the quartic-divisor algebra is zero-dimensional of rank
`8855`, but the public Magma service runs out of memory while constructing the
univariate eliminant.  Consequently the checkpoint retains

```text
certified_M23 = false.
```

The `p=9319` calculation is a computationally lean proxy for the same
four-set construction.  Transferring a future positive result to `MON13`
still requires the separate same-global-cover, good-tame-specialization, and
orientation gates; the rank checkpoint alone does not supply them.

Let

```text
K = F_(13^6)(t),
h(Z) in K[Z]
```

be the exact monic, irreducible, separable degree-23 sheet polynomial already
extracted from the exact cover.

## 0. Exact `p=9319` rank checkpoint and service boundary

The executable wrapper

```text
../petri-rr-map/magma_p9319_four_set_tower_probe.py
```

has a lean second-stage mode.  That mode does not reconstruct or retain the
curve, its function field, divisor objects, or passport objects.  Instead it
extracts the already printed exact polynomial `h(T,w)` from the hash-bound
Magma response

```text
petri-model-exactify/split-prime-p9319-a3983-s273-pair3-exact-passport.xml
SHA-256 d6526152c48e2541589cc7c056372ac03c7b035ab397fdab3672fa6bda892eda
```

and also requires the exact pair-3 source hash

```text
petri-model-exactify/split-prime-p9319-a3983-s273-pair3.m
SHA-256 b1113174bf4a1b78473e5534741e47e0e5772d0c89f0744a7269c879453ca455
```

before running.  The serialized polynomial expression itself has SHA-256

```text
38d0069465d628212cd2d0fba547fac079094fdd42d62d09ae1984f937b71826.
```

The lean Magma replay rechecks that `h` is monic, irreducible, separable, and
of degree 23.  It then forms the four coefficient equations of the remainder
of `h` modulo the universal monic quartic.  Exact output is:

```text
FOUR_SET_STAGE polynomial_ready
FOUR_SET_STAGE ideal_defined
Time: 48.880
FOUR_SET_STAGE zero_dimensional
Time: 0.000
FOUR_SET_QUOTIENT_DIMENSION 8855
```

Thus the rank-`8855` gate in Section 2 is now passed exactly for the `p=9319`
polynomial.  The immediately following operation was

```text
UnivariateEliminationIdealGenerator(J,4).
```

The service stopped it with the exact diagnostic

```text
Current total memory usage: 203.4MB, failed memory request: 299.2MB
System Error: User memory limit has been reached
```

Removing the live curve/passport construction did not change those numbers,
so this is the eliminant's service-memory boundary rather than retained
geometric objects.  No eliminant, selected quartic, degree-19 complement, or
`3+16` factorization was produced.

The canonical receipt is

```text
../petri-rr-map/magma-p9319-four-set-tower-lean.json
SHA-256 41539dabde214284919e65723718623031b989784d6b7df6ae937e86459d5563
```

Its exact structured output fields are:

```text
quotient_dimension                                  8855
calculator_reported_zero_dimensional_seconds        48.88
calculator_reported_memory_failure.operation        UnivariateEliminationIdealGenerator(J,4)
calculator_reported_memory_failure.current_total_mb 203.4
calculator_reported_memory_failure.failed_request_mb 299.2
completed_stages                                    [polynomial_ready,
                                                     ideal_defined,
                                                     zero_dimensional]
complement_factor_degrees                           null
certified_M23                                       false
```

The emitted lean Magma input has SHA-256
`38efa929aa5a780ebab137970751d9e99f78930da1ea2b8bc69b51e9b10d1977`;
the raw response has SHA-256
`31244edeb929a3e1ce3461fa03fcb185dfb2ab8b587565d4347c5da1083e9b34`.

Separately, the exact collision probe at `T=4` splits `h(T,w)` over
`F_(9319^6)` with factor degrees `1,2,2,3,3,6,6` and verifies that all
`8855` plain sums of four roots are distinct.  Therefore `Y=u1`, the plain
four-root sum, is already a valid generic separating coordinate at `p=9319`;
no Tschirnhausen ladder is needed there.  This collision result is not a
substitute for the missing generic eliminant and terminal factorization.

## 1. The direct degree-8855 algebra

Introduce four independent symmetric coordinates and the universal monic
quartic

```text
Q_U(Z) = Z^4 - u1 Z^3 + u2 Z^2 - u3 Z + u4.
```

Compute

```text
Rem_Z(h, Q_U) = r3 Z^3 + r2 Z^2 + r1 Z + r0
```

in `K[u1,u2,u3,u4][Z]`, and define

```text
A4 = K[u1,u2,u3,u4] / (r0,r1,r2,r3).
```

Over a splitting field of `h`, a point of `Spec(A4)` is exactly a monic
degree-four divisor of `h`, hence exactly an unordered four-subset of its 23
roots.  Since `h` is separable, this factor-selection scheme is finite etale
of rank

```text
binomial(23,4) = 8855.
```

The ordered construction adjoins four distinct roots and has degree

```text
23*22*21*20 = 212520 = 24*8855.
```

The direct construction above instead records only the elementary symmetric
functions of those roots.  It is the `S4`-invariant algebra of the ordered
tower.

The remainder should be generated by recurrence rather than generic symbolic
long division.  If `v_n` is the four-component coefficient vector of
`Z^n mod Q_U`, then for `n >= 4`,

```text
v_n = u1*v_(n-1) - u2*v_(n-2) + u3*v_(n-3) - u4*v_(n-4).
```

Then `r = sum_n coefficient(h,n)*v_n`.

## 2. Exact construction gates

The computation is accepted only if all of these gates pass.

1. Recheck that `h` is monic of degree 23, irreducible, and separable over
   `K`.
2. Compute a Groebner basis of `I4=(r0,r1,r2,r3)`.  Require a zero-dimensional
   quotient with exactly 8855 standard monomials.
3. Certify reducedness/etaleness.  The paper proof may invoke the standard
   separable-factor-selection lemma.  A machine receipt should additionally
   verify that the determinant of the `4 x 4` Jacobian
   `d(r0,r1,r2,r3)/d(u1,u2,u3,u4)` is a unit modulo `I4`, equivalently that
   `1` lies in `I4 + (det(Jacobian))`.
4. Select a deterministic primitive element.  One reproducible ladder is

   ```text
   Y_c = u1 + c*u2 + c^2*u3 + c^3*u4,
   c = 0,1,a,a+1,...
   ```

   in a fixed serialized order on `F_(13^6)`.  Stop at the first `c` for
   which the exact minimal polynomial `R4` of multiplication by `Y_c` has
   degree 8855 and satisfies `gcd(R4,R4')=1`.
5. The degree-8855 minimal polynomial is the no-collision certificate: all
   8855 geometric four-subsets have distinct `Y_c` values.  Do not accept
   repeated precision stability as a substitute.
6. Prove `R4` irreducible over `K`.  An exact generic function-field
   factorization is sufficient.  Alternatively, the already-established
   dichotomy `G_geom in {M23,A23}` proves irreducibility because both groups
   are transitive on unordered four-sets, once the no-collision gate has
   passed.  A constant-`t` specialization cannot be expected to give an
   irreducible degree-8855 polynomial: it records one Frobenius permutation,
   and neither candidate group contains an element acting as an 8855-cycle.

## 3. Recover the selected quartic

Run FGLM/change-of-primitive-element on

```text
I4 + (Y - Y_c(u1,u2,u3,u4)).
```

It must produce exact identities

```text
ui = Ui(Y) mod R4(Y),  i=1,2,3,4.
```

Define

```text
E = K[Y]/(R4),
Q(Z) = Z^4 - U1(Y) Z^3 + U2(Y) Z^2 - U3(Y) Z + U4(Y).
```

The following independent checks are mandatory:

```text
degree(Q) = 4,
Rem_Z(h,Q) = 0 in E[Z],
h = Q*H19 exactly,
degree(H19) = 19,
gcd(Q,H19) = 1,
Q and H19 are separable.
```

This is the exact recovery of the selected quartic.  The quartic lives over
the degree-8855 field `E`; it is not a record of four sequential root
adjunctions.

## 4. The decisive factorization

Factor `H19` completely over `E` and multiply all returned factors back.

Acceptance for `MON13=M23` requires exactly two multiplicity-one irreducible
factors of degrees

```text
3, 16.
```

The `A23` outcome is an irreducible degree-19 polynomial.  Given the existing
`M23/A23` dichotomy, the set stabilizer in `M23` has remaining-point orbits
`3+16`, whereas the set stabilizer in `A23` (and in a possible arithmetic
`S23` overgroup) is transitive on the remaining 19 points.  Thus an exact
`3+16` factorization rules out `A23` and certifies geometric monodromy `M23`.
The `23A/23B` orientation remains a separate gate.

## 5. Implementation and resource boundary

Do not call `GaloisSubgroup` to create this field: that route presupposes a
completed Galois-group computation and is circular for `MON13`.

A credible implementation is:

1. grevlex F4/F5 for the four-variable ideal;
2. a sparse multiplication operator for `Y_c` on the 8855 standard
   monomials;
3. block Wiedemann (with independently verified minimal polynomial) rather
   than dense linear algebra;
4. sparse FGLM or transposed-Krylov recovery of the four `Ui(Y)`;
5. exact factorization of the degree-19 cofactor over the resulting simple
   algebraic function field.

A single dense `8855 x 8855` matrix has 78,411,025 entries.  Several dense
matrices with CAS object overhead can consume tens of gigabytes, so dense
FGLM is not the default.  The Groebner/change-of-order phase is likely the
dominant cost.  For the measured `p=9319` replay, start the emitted script in
licensed Magma with at least 4 GB available.  For the more expensive direct
`F_(13^6)(t)` calculation, begin on a checkpointed 32-GB worker and keep a
128-GB worker as the escalation tier rather than treating it as a measured
minimum.  Provision hours, with an hours-to-days stop window.  The public
Magma calculator did complete the zero-dimensional/rank gate within its
60-second limit, but its next `299.2MB` request at `203.4MB` current use was
denied, so it cannot finish the eliminant stage.

Modular evaluation/interpolation in `t` is permitted as an accelerator only
if the reconstructed generic objects are finally checked by exact reduction
of all four remainder equations over `K`, exact minimal-polynomial identities,
and exact multiplication-back of the `3+16` factors.  Agreement across a few
specializations is not the certificate.
