#!/usr/bin/env python3
"""Recognize a conditional exact degree-23 map over F(sqrt(-23)).

This is the exactification gate following ``numerical_map_probe.py``.  It
uses the graph equation

    N = phi D,       N,D in H^0(5K),

whose numerical solution space has dimension four.  A representative is
chosen without inspecting any prospective algebraic coefficients:

* among the 27 denominator coordinates, choose the four rows whose 4-by-4
  restrictions have the largest common smallest singular value for labels
  1 and 2;
* set three of those coordinates to zero and one to one, selecting the unit
  row which minimizes the largest numerical coefficient.

The row choice is made at ``--discovery-length`` and frozen for the
independent ``--check-length`` run.  Coefficients in

    L = Q(alpha,s),  g(alpha)=0,  s^2=-23,

are recognized after the scalar semilinear split

    tau(c_1) = conjugate(c_2),
    A=(c_1+tau(c_1))/2, B=(c_1-tau(c_1))/(2s).

There is no graph-half swap in this scalar-field operation.  The separate
law ``tau(Graph(phi))=Graph(1/phi)`` does swap N and D and remains a distinct
target/source normalization gate.  Each F-basis relation must be returned at
two LLL thresholds and agree with both numerical precisions.

The output is deliberately called a *candidate*: integer relations plus
high-precision q-series matching do not prove the q-series are an algebraic
Hurwitz cover.  Exact curve/map/passport verification is a separate gate.
"""

from __future__ import annotations

import argparse
import importlib
import itertools
import json
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from cypari2 import Pari
from mpmath import fabs, lu_solve, matrix, mp, mpc, mpf

import numerical_map_probe as numerical


COMPACT_ASCENDING = (4, 5, 1, -3, -3, -1, 1)
RELATION_DIGITS = (190, 180, 170, 160, 150, 140, 130, 120)
ZERO_TOL = mpf("1e-155")
CHECK_TOL = mpf("1e-145")
MAX_RELATION_INTEGER = 10**60


def fixed_sections(
    basis: list[list[mpc]],
    canonical_curve: Any,
    monomials: list[tuple[int, ...]],
    length: int,
) -> tuple[list[list[mpc]], mpc]:
    """Evaluate a previously selected monomial basis of H^0(5K)."""

    indices, coefficients = canonical_curve.quadric(basis, n_use=length)
    raw_quadric = dict(zip(indices, coefficients, strict=True))
    lam = raw_quadric[(1, 2)] / raw_quadric[(2, 2)]
    normalized = [
        [value / lam**index for value in basis[index]][:length]
        for index in range(4)
    ]
    return (
        [
            numerical.monomial_series(normalized, monomial, length)
            for monomial in monomials
        ],
        lam,
    )


def graph_kernel(
    p200: Any,
    canonical_curve: Any,
    pipeline_dir: Path,
    path: Path,
    length: int,
    fixed_monomials_: list[tuple[int, ...]] | None = None,
) -> dict[str, Any]:
    """Return the four-dimensional kernel of [sections|-phi*sections]."""

    basis = p200.load_p200_basis(pipeline_dir, path, nuse=length)
    if fixed_monomials_ is None:
        monomials, sections, lam = numerical.section_basis(
            basis, canonical_curve, length
        )
    else:
        monomials = fixed_monomials_
        sections, lam = fixed_sections(
            basis, canonical_curve, monomials, length
        )
    _, _, kappas = numerical.triangle_vertices_and_kappas()
    if len(kappas) != 1:
        raise ValueError(f"expected one TrianglePhi scale, found {len(kappas)}")
    phi = numerical.phi_series(kappas[0], length)
    phi_sections = [
        numerical.series_mul(phi, section, length) for section in sections
    ]
    columns = sections + [[-value for value in series] for series in phi_sections]
    kernel, rank = canonical_curve.nullspace(
        columns, length, expect=4, tol=mpf("1e-150")
    )
    residuals = [
        numerical.kernel_residual(columns, vector, length) for vector in kernel
    ]
    if rank != 50:
        raise ValueError(f"graph matrix rank {rank}, expected 50")
    return {
        "kernel": kernel,
        "rank": rank,
        "residual": max(residuals),
        "monomials": list(monomials),
        "lambda": lam,
        "kappa": kappas[0],
    }


def swap_graph_halves(kernel: list[list[mpc]]) -> list[list[mpc]]:
    """Align label 2 with label 1 when phi^tau=1/phi."""

    result = []
    for vector in kernel:
        if len(vector) != 54:
            raise ValueError(f"unexpected graph-vector length {len(vector)}")
        result.append(vector[27:] + vector[:27])
    return result


def normalized_kernel_array(kernel: list[list[mpc]]) -> np.ndarray:
    array = np.array(
        [
            [complex(kernel[column][row]) for column in range(4)]
            for row in range(54)
        ],
        dtype=np.complex128,
    )
    norms = np.linalg.norm(array, axis=0)
    if np.any(norms == 0) or not np.all(np.isfinite(norms)):
        raise ValueError(f"invalid kernel column norms {norms}")
    return array / norms


def solve_denominator_gauge(
    kernel: list[list[mpc]],
    selected_denominator_rows: tuple[int, ...],
    normalization_denominator_row: int,
) -> list[mpc]:
    """Impose four rational coordinate conditions on the D half."""

    absolute_rows = tuple(27 + row for row in selected_denominator_rows)
    system = matrix(
        [[kernel[column][row] for column in range(4)] for row in absolute_rows]
    )
    target = matrix(
        [
            1 if row == normalization_denominator_row else 0
            for row in selected_denominator_rows
        ]
    )
    weights = lu_solve(system, target)
    return [
        sum(kernel[column][row] * weights[column] for column in range(4))
        for row in range(54)
    ]


def choose_denominator_gauge(
    kernels: dict[int, list[list[mpc]]],
    monomials: list[tuple[int, ...]],
) -> tuple[tuple[int, ...], int, dict[int, list[mpc]], dict[str, Any]]:
    """Choose a deterministic, simultaneous, denominator-only gauge."""

    arrays = {label: normalized_kernel_array(kernel) for label, kernel in kernels.items()}
    best_rows: tuple[int, ...] | None = None
    best_singular = -1.0
    for rows in itertools.combinations(range(27), 4):
        absolute = [27 + row for row in rows]
        smallest = min(
            float(np.linalg.svd(array[absolute, :], compute_uv=False)[-1])
            for array in arrays.values()
        )
        if smallest > best_singular:
            best_singular = smallest
            best_rows = rows
    if best_rows is None:
        raise ValueError("no denominator gauge rows found")

    best_norm: int | None = None
    best_vectors: dict[int, list[mpc]] | None = None
    best_max = float("inf")
    trials = []
    for norm_row in best_rows:
        vectors = {
            label: solve_denominator_gauge(kernel, best_rows, norm_row)
            for label, kernel in kernels.items()
        }
        largest = max(
            float(max(fabs(value) for value in vector))
            for vector in vectors.values()
        )
        trials.append(
            {
                "normalization_denominator_row": norm_row,
                "normalization_monomial": str(monomials[norm_row]),
                "largest_coefficient": largest,
            }
        )
        if largest < best_max:
            best_max = largest
            best_norm = norm_row
            best_vectors = vectors
    if best_norm is None or best_vectors is None:
        raise ValueError("no denominator normalization found")
    return (
        best_rows,
        best_norm,
        best_vectors,
        {
            "selection_rule": (
                "maximize the minimum, across label1 and label2, of "
                "the smallest singular value of a four-row restriction of "
                "the column-normalized graph kernels; then minimize the "
                "largest coefficient among the four possible unit-row gauges"
            ),
            "selected_denominator_rows": list(best_rows),
            "selected_denominator_monomials": [
                str(monomials[row]) for row in best_rows
            ],
            "minimum_selected_submatrix_singular_value": best_singular,
            "normalization_denominator_row": best_norm,
            "normalization_monomial": str(monomials[best_norm]),
            "zero_denominator_rows": [row for row in best_rows if row != best_norm],
            "trials": trials,
        },
    )


def evaluate_component(coefficients: list[Fraction], alpha: mpc) -> mpc:
    result = mpc(0)
    for coefficient in reversed(coefficients):
        result = result * alpha + mpf(coefficient.numerator) / coefficient.denominator
    return result


def evaluate_l_element(
    a_coefficients: list[Fraction],
    b_coefficients: list[Fraction],
    alpha: mpc,
    s: mpc,
) -> mpc:
    return evaluate_component(a_coefficients, alpha) + s * evaluate_component(
        b_coefficients, alpha
    )


def relation_to_coefficients(relation: Any) -> tuple[list[Fraction], list[Fraction]]:
    integers = [int(entry) for entry in relation]
    if len(integers) != 13 or integers[12] == 0:
        raise ValueError(f"invalid L relation {integers}")
    if max(abs(value) for value in integers) > MAX_RELATION_INTEGER:
        raise ValueError(
            "relation exceeds the fail-closed height bound: "
            f"{max(abs(value) for value in integers)}"
        )
    a = [Fraction(-integers[index], integers[12]) for index in range(6)]
    b = [Fraction(-integers[6 + index], integers[12]) for index in range(6)]
    return a, b


def recognize_f_element(
    pari: Pari,
    p200: Any,
    value: mpc,
    alpha: mpc,
) -> tuple[list[Fraction], dict[str, Any]]:
    """Require the same direct compact-sextic relation at two thresholds."""

    alpha_pari = p200.pari_complex(pari, alpha)
    z_pari = p200.pari_complex(pari, value)
    basis = [alpha_pari**index for index in range(6)]
    attempts = []
    recognized: list[list[Fraction]] = []
    for digits in RELATION_DIGITS:
        relation = pari.lindep(pari(basis + [z_pari]), digits)
        integers = [int(entry) for entry in relation]
        if len(integers) != 7 or integers[6] == 0:
            raise ValueError(f"invalid F relation {integers}")
        if max(abs(entry) for entry in integers) > MAX_RELATION_INTEGER:
            raise ValueError(
                "relation exceeds the fail-closed height bound: "
                f"{max(abs(entry) for entry in integers)}"
            )
        coefficients = [
            Fraction(-integers[index], integers[6]) for index in range(6)
        ]
        error = fabs(evaluate_component(coefficients, alpha) - value)
        attempts.append(
            {
                "digits": digits,
                "relation_height": max(abs(entry) for entry in integers),
                "embedding_error": str(error),
            }
        )
        recognized.append(coefficients)
        if len(recognized) >= 2 and recognized[-1] == recognized[-2]:
            stable_attempts = attempts[-2:]
            if max(mpf(item["embedding_error"]) for item in stable_attempts) >= CHECK_TOL:
                raise ValueError(
                    f"stable F relation embedding error too large: {stable_attempts}"
                )
            return recognized[-1], {
                "method": "direct_PARI_lindep_in_compact_sextic_basis",
                "threshold_stability": (
                    "identical exact coefficient vectors at two adjacent thresholds"
                ),
                "stable_thresholds": [item["digits"] for item in stable_attempts],
                "attempts": attempts,
            }
    raise ValueError(f"F relation changed at every tested threshold: {attempts}")


def recognize_semilinear_pair(
    pari: Pari,
    p200: Any,
    value: mpc,
    tau_value: mpc,
    alpha: mpc,
    s: mpc,
) -> tuple[list[Fraction], list[Fraction], dict[str, Any]]:
    """Split c and tau(c) first, then recognize A(alpha), B(alpha).

    A single complex embedding of the 12-dimensional field L cannot by
    itself select a trustworthy small relation.  The second orientation,
    without a graph-half swap, supplies tau(c), making the two F-valued components
    A=(c+tau(c))/2 and B=(c-tau(c))/(2s) separately recognizable.
    """

    a_value = (value + tau_value) / 2
    b_value = (value - tau_value) / (2 * s)
    a, a_receipt = recognize_f_element(pari, p200, a_value, alpha)
    b, b_receipt = recognize_f_element(pari, p200, b_value, alpha)
    return a, b, {
        "method": "semilinear_split_then_two_compact_sextic_lindep_relations",
        "A": a_receipt,
        "B": b_receipt,
    }


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def serialize_element(
    a: list[Fraction], b: list[Fraction]
) -> dict[str, list[str]]:
    return {
        "A_coefficients_ascending": [fraction_text(value) for value in a],
        "B_coefficients_ascending": [fraction_text(value) for value in b],
    }


def exact_rational(value: int) -> tuple[list[Fraction], list[Fraction]]:
    return (
        [Fraction(value), *([Fraction(0)] * 5)],
        [Fraction(0)] * 6,
    )


def parse_complex_text(text: str) -> mpc:
    """Parse the long ``mp.nstr`` complex form without a float round-trip."""

    match = re.fullmatch(
        r"\(\s*([+-]?[0-9.]+(?:[Ee][+-]?[0-9]+)?)\s*"
        r"([+-])\s*([0-9.]+(?:[Ee][+-]?[0-9]+)?)j\s*\)",
        text,
    )
    if match is None:
        raise ValueError(f"cannot parse complex text {text!r}")
    imaginary = mpf(match.group(3))
    if match.group(2) == "-":
        imaginary = -imaginary
    return mpc(mpf(match.group(1)), imaginary)


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
    parser.add_argument(
        "--label2",
        type=Path,
        default=Path("/private/tmp/m23-p200-recovered/label2-out-p200.m"),
    )
    parser.add_argument("--discovery-length", type=int, default=100)
    parser.add_argument("--check-length", type=int, default=120)
    parser.add_argument(
        "--field-receipt",
        type=Path,
        default=Path(
            "outputs/m23-proof-closeout-20260902/artifacts/"
            "petri-p200-field-recognition.json"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "outputs/m23-proof-closeout-20260902/agents/petri-f13/"
            "exact-map-candidate.json"
        ),
    )
    args = parser.parse_args()
    if args.discovery_length >= args.check_length:
        raise ValueError("discovery length must be smaller than check length")
    if args.discovery_length < 80:
        raise ValueError("discovery length must be at least 80")

    mp.dps = 220
    pipeline_dir = args.pipeline_dir.resolve()
    sys.path.insert(0, str(Path("outputs/m23-proof-closeout-20260902").resolve()))
    sys.path.insert(0, str(pipeline_dir))
    p200 = importlib.import_module("petri_p200_field_recognition")
    canonical_curve = importlib.import_module("canonical_curve")
    p200.mp.dps = 220
    canonical_curve.mp.dps = 220

    # Discovery run.  Only label 1 selects the degree-five monomial basis.
    discovery_1 = graph_kernel(
        p200,
        canonical_curve,
        pipeline_dir,
        args.label1.resolve(),
        args.discovery_length,
    )
    monomials = discovery_1["monomials"]
    discovery_2_raw = graph_kernel(
        p200,
        canonical_curve,
        pipeline_dir,
        args.label2.resolve(),
        args.discovery_length,
        monomials,
    )
    discovery_kernels = {1: discovery_1["kernel"], 2: discovery_2_raw["kernel"]}
    selected_rows, norm_row, discovery_vectors, gauge = choose_denominator_gauge(
        discovery_kernels, monomials
    )
    print(
        "discovery gauge",
        selected_rows,
        "unit",
        norm_row,
        flush=True,
    )

    # Independent longer-jet run: reuse both the monomial basis and gauge.
    check_1 = graph_kernel(
        p200,
        canonical_curve,
        pipeline_dir,
        args.label1.resolve(),
        args.check_length,
        monomials,
    )
    check_2_raw = graph_kernel(
        p200,
        canonical_curve,
        pipeline_dir,
        args.label2.resolve(),
        args.check_length,
        monomials,
    )
    check_kernels = {1: check_1["kernel"], 2: check_2_raw["kernel"]}
    check_vectors = {
        label: solve_denominator_gauge(kernel, selected_rows, norm_row)
        for label, kernel in check_kernels.items()
    }
    precision_errors = {
        str(label): str(
            max(
                fabs(discovery_vectors[label][row] - check_vectors[label][row])
                for row in range(54)
            )
        )
        for label in (1, 2)
    }
    if max(mpf(value) for value in precision_errors.values()) >= CHECK_TOL:
        raise ValueError(f"gauge is not precision-stable: {precision_errors}")
    print("precision errors", precision_errors, flush=True)

    field_receipt = json.loads(args.field_receipt.read_text())
    alpha = parse_complex_text(field_receipt["field"]["label1_compact_root"])
    compact_residual = fabs(
        sum(mpf(coefficient) * alpha**index for index, coefficient in enumerate(COMPACT_ASCENDING))
    )
    if compact_residual >= mpf("1e-180"):
        raise ValueError(f"compact-root residual {compact_residual}")
    s = mpc(0, mp.sqrt(23))
    pari = Pari()
    pari("\\p 250")

    coefficients = []
    recognition_errors = []
    exact_gauge_rows = {27 + row for row in selected_rows}
    for row in range(54):
        high_1 = check_vectors[1][row]
        high_2_tau_target = check_vectors[2][row].conjugate()
        low_1 = discovery_vectors[1][row]
        low_2_tau_target = discovery_vectors[2][row].conjugate()
        if row in exact_gauge_rows:
            expected = 1 if row == 27 + norm_row else 0
            a, b = exact_rational(expected)
            recognition = {
                "method": "exact_denominator_gauge",
                "value": str(expected),
            }
        elif max(fabs(value) for value in (high_1, high_2_tau_target, low_1, low_2_tau_target)) < ZERO_TOL:
            a, b = exact_rational(0)
            recognition = {
                "method": "precision_stable_numerical_zero",
                "threshold": str(ZERO_TOL),
            }
        else:
            a, b, recognition = recognize_semilinear_pair(
                pari,
                p200,
                high_1,
                high_2_tau_target,
                alpha,
                s,
            )

        value_1 = evaluate_l_element(a, b, alpha, s)
        value_tau = evaluate_l_element(a, b, alpha, -s)
        errors = {
            "label1_discovery": fabs(value_1 - low_1),
            "label1_check": fabs(value_1 - high_1),
            "label2_scalar_tau_discovery": fabs(value_tau - low_2_tau_target),
            "label2_scalar_tau_check": fabs(value_tau - high_2_tau_target),
        }
        if max(errors.values()) >= CHECK_TOL:
            raise ValueError(
                f"row {row} fails precision/semilinear embedding checks: {errors}"
            )
        recognition_errors.append(
            {
                "graph_row": row,
                **{name: str(error) for name, error in errors.items()},
            }
        )
        coefficients.append(
            {
                "graph_row": row,
                "half": "N" if row < 27 else "D",
                "section_row": row if row < 27 else row - 27,
                "monomial": str(monomials[row if row < 27 else row - 27]),
                "exact": serialize_element(a, b),
                "recognition": recognition,
            }
        )
        print("recognized row", row, flush=True)

    max_recognition_error = max(
        mpf(item[name])
        for item in recognition_errors
        for name in (
            "label1_discovery",
            "label1_check",
            "label2_scalar_tau_discovery",
            "label2_scalar_tau_check",
        )
    )
    payload = {
        "status": "conditional_exact_map_candidate",
        "field": {
            "compact_polynomial_coefficients_ascending": list(COMPACT_ASCENDING),
            "quadratic_generator_relation": "s^2+23=0",
            "label1_alpha_embedding": mp.nstr(alpha, 215),
            "compact_root_residual": str(compact_residual),
        },
        "semilinear_convention": {
            "scalar_action": "tau(c_1)=conjugate(c_2), with no graph-half swap",
            "separate_target_action": (
                "tau(Graph(phi))=Graph(1/phi) swaps (N,D), and is a distinct "
                "target/source normalization gate not used in scalar recognition"
            ),
            "coefficient_pairing": (
                "label1 c=A(alpha)+sB(alpha); conjugate(label2) "
                "must equal A(alpha)-sB(alpha)"
            ),
        },
        "numerical_reconstruction": {
            "discovery_length": args.discovery_length,
            "check_length": args.check_length,
            "degree_five_monomial_basis": [str(value) for value in monomials],
            "discovery_graph_ranks": {
                "label1": discovery_1["rank"],
                "label2": discovery_2_raw["rank"],
            },
            "check_graph_ranks": {
                "label1": check_1["rank"],
                "label2": check_2_raw["rank"],
            },
            "graph_kernel_dimensions": {"label1": 4, "label2": 4},
            "graph_residuals": {
                "discovery_label1": str(discovery_1["residual"]),
                "discovery_label2": str(discovery_2_raw["residual"]),
                "check_label1": str(check_1["residual"]),
                "check_label2": str(check_2_raw["residual"]),
            },
            "coordinate_lambdas": {
                "discovery_label1": mp.nstr(discovery_1["lambda"], 80),
                "discovery_label2": mp.nstr(discovery_2_raw["lambda"], 80),
                "check_label1": mp.nstr(check_1["lambda"], 80),
                "check_label2": mp.nstr(check_2_raw["lambda"], 80),
            },
            "triangle_kappa": mp.nstr(discovery_1["kappa"], 215),
            "precision_errors_after_frozen_gauge": precision_errors,
        },
        "denominator_gauge": gauge,
        "recognition": {
            "basis": "1,alpha,...,alpha^5,s,s*alpha,...,s*alpha^5",
            "lindep_thresholds": list(RELATION_DIGITS),
            "height_bound": str(MAX_RELATION_INTEGER),
            "embedding_tolerance": str(CHECK_TOL),
            "max_four_way_embedding_error": str(max_recognition_error),
            "per_coefficient_errors": recognition_errors,
        },
        "map": {
            "interpretation": "phi=N/D in the normalized canonical coordinates",
            "coefficients": coefficients,
        },
        "boundary": {
            "proves": (
                "a reproducibly gauged exact L-valued candidate (N,D), with "
                "identical integer relations at two LLL thresholds, agreement "
                "at two jet lengths, and the label2 scalar semilinear check"
            ),
            "does_not_prove": (
                "that the analytic q-series are exact algebraic Hurwitz covers, "
                "that Q=C=N/D defines a smooth degree-23 cover, that its exact "
                "passport is (2A,23A,23B), that the separate N/D block-swap "
                "target normalization is compatible, or that reduction is good at 13"
            ),
            "next_gate": (
                "reduce the exact Q,C,N,D at the prime (13,s-4), then verify "
                "smoothness, degree, the complete three-fiber divisor passport, "
                "and the Petri value in exact finite-field arithmetic"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)
    print("max recognition error", max_recognition_error)


if __name__ == "__main__":
    main()
