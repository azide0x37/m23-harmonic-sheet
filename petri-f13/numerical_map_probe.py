#!/usr/bin/env python3
"""Recover the degree-23 triangle map numerically from the p200 1-form basis.

The Petri pipeline stopped after reconstructing the canonical curve.  The
ambient triangle group nevertheless supplies its Hauptmodul explicitly: the
inverse of a hypergeometric Schwarz map (``TrianglePhi`` in the upstream
Belyi package).  This probe ports just that formula for Delta(23,2,23), then
solves

    N(w) = phi(w) D(w)

with N,D in H^0(X,5K).  Degree 5 is the first uniform choice for which a
denominator can absorb the unique order-23 pole: deg(5K)=30 and
``h^0(5K-23P)=4``.  Hence the correct kernel is expected to have dimension 4.

This file is discovery-only.  The returned complex kernel must still be
recognized over the exact coefficient field and verified after reduction.
"""

from __future__ import annotations

import argparse
import importlib
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from mpmath import fabs, gamma, mp, mpc, mpf, pi, sqrt
from sympy import QQ
from sympy.polys.ring_series import (
    rs_series_inversion,
    rs_series_reversion,
    rs_trunc,
)
from sympy.polys.rings import ring


MONOMIALS_2 = tuple(itertools.combinations_with_replacement(range(4), 2))
MONOMIALS_3 = tuple(itertools.combinations_with_replacement(range(4), 3))
MONOMIALS_5 = tuple(itertools.combinations_with_replacement(range(4), 5))


def hypergeometric_coefficients(
    a: Fraction, b: Fraction, c: Fraction, count: int
) -> list[Fraction]:
    coefficients = [Fraction(1)]
    for n in range(1, count):
        coefficients.append(
            coefficients[-1]
            * (a + n - 1)
            * (b + n - 1)
            / ((c + n - 1) * n)
        )
    return coefficients


def triangle_phi_coefficients(count: int) -> list[Fraction]:
    """Return u(v) with v=(kappa*w)^23 and phi(w)=u(v)."""

    a = Fraction(1, 4)
    b = Fraction(27, 92)
    c = Fraction(24, 23)
    a2 = Fraction(19, 92)
    b2 = Fraction(1, 4)
    c2 = Fraction(22, 23)
    f1 = hypergeometric_coefficients(a, b, c, count)
    f2 = hypergeometric_coefficients(a2, b2, c2, count)
    polynomial_ring, u, v = ring("u, v", QQ)
    f1_polynomial = sum(
        (
            polynomial_ring.domain.convert(value.numerator)
            / polynomial_ring.domain.convert(value.denominator)
        )
        * u**index
        for index, value in enumerate(f1)
    )
    f2_polynomial = sum(
        (
            polynomial_ring.domain.convert(value.numerator)
            / polynomial_ring.domain.convert(value.denominator)
        )
        * u**index
        for index, value in enumerate(f2)
    )
    # Schwarz map w = t^(1/23) F1(t)/F2(t), so
    # v=w^23=t*(F1/F2)^23.  Revert this ordinary series.
    ratio = rs_trunc(
        f1_polynomial * rs_series_inversion(f2_polynomial, u, count), u, count
    )
    forward = rs_trunc(u * ratio**23, u, count)
    inverse = rs_series_reversion(forward, u, count, v)
    coefficients = []
    for exponent in range(count):
        coefficient = inverse.coeff(v**exponent)
        coefficients.append(Fraction(int(coefficient.numerator), int(coefficient.denominator)))
    if coefficients[0] != 0 or coefficients[1] != 1:
        raise AssertionError("unexpected normalized TrianglePhi series")
    return coefficients


def triangle_vertices_and_kappas() -> tuple[mpc, mpc, list[mpc]]:
    a = 23
    b = 2
    c = 23
    ca, sa = mp.cos(pi / a), mp.sin(pi / a)
    cb, sb = mp.cos(pi / b), mp.sin(pi / b)
    cc = mp.cos(pi / c)
    ell = (ca * cb + cc) / (sa * sb)
    tau = ell + sqrt(ell * ell - 1)
    da = mp.matrix([[ca, sa], [-sa, ca]])
    db = mp.matrix([[cb, tau * sb], [-sb / tau, cb]])
    dc = (da * db) ** -1
    c11, c12, c21, c22 = dc[0, 0], dc[0, 1], dc[1, 0], dc[1, 1]
    discriminant = sqrt((c22 - c11) ** 2 + 4 * c21 * c12)
    roots = [
        ((c11 - c22) + discriminant) / (2 * c21),
        ((c11 - c22) - discriminant) / (2 * c21),
    ]
    vc = next(root for root in roots if mp.im(root) > 0)
    vc_mirror = -mp.conj(vc)
    z0 = mpc(0, 1)

    A = mpf(1) / 4
    B = mpf(27) / 92
    C = mpf(24) / 23
    A2 = mpf(19) / 92
    B2 = mpf(1) / 4
    C2 = mpf(22) / 23

    def at_one(aa: mpf, bb: mpf, cc_: mpf) -> mpf:
        return gamma(cc_) * gamma(cc_ - aa - bb) / (gamma(cc_ - aa) * gamma(cc_ - bb))

    gamma_ratio = at_one(A2, B2, C2) / at_one(A, B, C)
    # FundamentalDomain returns [z_a, z_c, z_b, z_c*], and TrianglePhi
    # specifically uses entries 1 and 3 (Magma indexing): z_a=i and
    # z_b=i*tau.  The c-vertex is *not* used in kappa.
    zb = mpc(0, tau)
    kappa = (zb - z0) / (zb - mp.conj(z0)) * gamma_ratio
    kappas = [kappa]
    return vc, vc_mirror, kappas


def phi_series(kappa: mpc, length: int) -> list[mpc]:
    coefficients = triangle_phi_coefficients(length // 23 + 2)
    result = [mpc(0) for _ in range(length)]
    for exponent, coefficient in enumerate(coefficients):
        degree = 23 * exponent
        if degree >= length:
            break
        # Upstream evaluates the canonical forms as f(kappa*w_new) while
        # leaving TrianglePhi in w_new.  Therefore in the solver's original
        # disc coordinate w_old the Hauptmodul is phi(w_old/kappa).
        result[degree] = (
            mpf(coefficient.numerator) / coefficient.denominator / kappa**degree
        )
    return result


def series_mul(left: list[mpc], right: list[mpc], length: int) -> list[mpc]:
    return [
        sum(left[index] * right[degree - index] for index in range(degree + 1))
        for degree in range(length)
    ]


def monomial_series(
    coordinates: list[list[mpc]], monomial: tuple[int, ...], length: int
) -> list[mpc]:
    result = [mpc(0) for _ in range(length)]
    result[0] = 1
    for variable in monomial:
        result = series_mul(result, coordinates[variable], length)
    return result


def numeric_pivots(matrix: np.ndarray, tolerance: float = 1e-10) -> list[int]:
    work = matrix.copy()
    row = 0
    pivots = []
    for column in range(work.shape[1]):
        candidates = np.abs(work[row:, column])
        if not candidates.size:
            break
        offset = int(np.argmax(candidates))
        if candidates[offset] < tolerance:
            continue
        pivot = row + offset
        work[[row, pivot], :] = work[[pivot, row], :]
        work[row, :] /= work[row, column]
        for other in range(work.shape[0]):
            if other != row:
                work[other, :] -= work[other, column] * work[row, :]
        pivots.append(column)
        row += 1
    return pivots


def section_basis(
    basis: list[list[mpc]], canonical_curve: Any, length: int
) -> tuple[list[tuple[int, ...]], list[list[mpc]], mpc]:
    indices, coefficients = canonical_curve.quadric(basis, n_use=length)
    raw_quadric = dict(zip(indices, coefficients, strict=True))
    lam = raw_quadric[(1, 2)] / raw_quadric[(2, 2)]
    normalized = [
        [value / lam**index for value in basis[index]][:length]
        for index in range(4)
    ]
    all_series = [
        monomial_series(normalized, monomial, length) for monomial in MONOMIALS_5
    ]
    # The 56 degree-five monomials have a 29-dimensional ideal of relations,
    # hence rank 27 on a (2,3) canonical complete intersection.
    evaluation = np.array(
        [[complex(series[degree]) for series in all_series] for degree in range(length)],
        dtype=np.complex128,
    )
    pivots = numeric_pivots(evaluation)
    if len(pivots) != 27:
        raise AssertionError(f"degree-five section rank {len(pivots)} != 27")
    return [MONOMIALS_5[index] for index in pivots], [all_series[index] for index in pivots], lam


def kernel_residual(
    columns: list[list[mpc]], vector: list[mpc], length: int
) -> mpf:
    return max(
        fabs(sum(vector[column] * columns[column][degree] for column in range(len(columns))))
        for degree in range(length)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pipeline-dir",
        type=Path,
        default=Path("outputs/hm-degeneration-20260819/phaseB-open"),
    )
    parser.add_argument(
        "--label1",
        type=Path,
        default=Path("/private/tmp/m23-p200-recovered/label1-out-p200.m"),
    )
    parser.add_argument("--length", type=int, default=120)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "outputs/m23-proof-closeout-20260902/agents/petri-f13/"
            "numerical-map-probe.json"
        ),
    )
    args = parser.parse_args()
    mp.dps = 220
    if args.length < 80:
        raise ValueError("length must be at least 80")

    sys.path.insert(0, str(Path("outputs/m23-proof-closeout-20260902").resolve()))
    sys.path.insert(0, str(args.pipeline_dir.resolve()))
    p200 = importlib.import_module("petri_p200_field_recognition")
    canonical_curve = importlib.import_module("canonical_curve")
    p200.mp.dps = 220
    canonical_curve.mp.dps = 220
    basis = p200.load_p200_basis(
        args.pipeline_dir.resolve(), args.label1.resolve(), nuse=args.length
    )
    monomials, sections, lam = section_basis(basis, canonical_curve, args.length)
    vc, vc_mirror, kappas = triangle_vertices_and_kappas()
    trials = []
    successful_kernel = None
    for orientation, kappa in enumerate(kappas):
        phi = phi_series(kappa, args.length)
        phi_sections = [series_mul(phi, section, args.length) for section in sections]
        columns = sections + [
            [-value for value in series] for series in phi_sections
        ]
        kernel, rank = canonical_curve.nullspace(
            columns,
            args.length,
            expect=None,
            tol=mpf("1e-150"),
        )
        residuals = [kernel_residual(columns, vector, args.length) for vector in kernel]
        trials.append(
            {
                "orientation": orientation,
                "kappa": mp.nstr(kappa, 215),
                "kappa_abs": mp.nstr(abs(kappa), 50),
                "rank": rank,
                "kernel_dimension": len(kernel),
                "max_kernel_residual": (
                    mp.nstr(max(residuals), 30) if residuals else None
                ),
            }
        )
        print(
            "orientation",
            orientation,
            "rank",
            rank,
            "kernel",
            len(kernel),
            "residual",
            max(residuals) if residuals else None,
            flush=True,
        )
        if len(kernel) == 4:
            successful_kernel = {
                "orientation": orientation,
                "vectors": [
                    [mp.nstr(value, 215) for value in vector] for vector in kernel
                ],
            }

    payload = {
        "status": "numerical_discovery_probe",
        "length": args.length,
        "degree_five_section_dimension": len(sections),
        "degree_five_monomial_basis": [str(monomial) for monomial in monomials],
        "canonical_coordinate_lambda": mp.nstr(lam, 215),
        "triangle": {
            "signature": [23, 2, 23],
            "upper_c_vertex": mp.nstr(vc, 215),
            "reflected_upper_c_vertex": mp.nstr(vc_mirror, 215),
            "phi_normalization": "0 at a-vertex, 1 at b-vertex, infinity at c-vertex",
        },
        "trials": trials,
        "kernel": successful_kernel,
        "boundary": {
            "proves": "nothing exact; this only reconstructs a numerical candidate map",
            "next_gate": (
                "recognize one kernel pencil over F(sqrt(-23)), reduce it at 13, "
                "and verify the exact divisor/passport identities"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
