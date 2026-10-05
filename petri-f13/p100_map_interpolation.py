#!/usr/bin/env python3
"""Attempt six-embedding exactification of the reconstructed Belyi map.

The two p200 orientations are enough to make a numerical graph kernel very
stable, but not enough to distinguish a genuine sextic-field coefficient
from successive high-quality LLL approximants.  This script instead uses all
six nonharmonic p100 solutions.

For every label ``l`` it selects ``c_l`` by four rational conditions on the
denominator D.  Scalar-field conjugation is then

    tau(c_l) = conjugate(c_kappa(l))

with *no* graph-half swap.  Thus A=(c+tau(c))/2 and
B=(c-tau(c))/(2s) can be interpolated over the six
compact-sextic roots.  Only rational power-basis coefficients that are
stable at two PSLQ tolerances and reproduce every embedding are accepted.
If this succeeds, label 1 and 2 are recomputed from the independent p200
files as an out-of-sample precision check.

The block swap is a separate target compatibility law
``tau(Graph(phi))=Graph(1/phi)``.  It is not used in scalar recognition and
is deliberately left as a distinct fail-closed gate.

The result remains conditional until Q,C,N,D receive exact algebraic
smoothness, divisor/passport, and reduction checks.
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from mpmath import fabs, lu_solve, matrix, mp, mpc, mpf, pslq

import exact_map_candidate as p200_map
import numerical_map_probe as numerical


NONHARMONIC = (1, 2, 3, 4, 5, 7)
KAPPA = {1: 2, 2: 1, 3: 3, 4: 4, 5: 7, 6: 6, 7: 5}
GAUGE_MONOMIALS = (
    (0, 0, 1, 2, 3),
    (0, 1, 1, 1, 2),
    (0, 1, 1, 2, 3),
    (0, 1, 2, 2, 2),
)
GAUGE_UNIT_MONOMIAL = (0, 1, 2, 2, 2)
PSLQ_TOLERANCES = (mpf("1e-70"), mpf("1e-64"))
PSLQ_BOUND = 10**50
EMBEDDING_TOL = mpf("1e-57")


@dataclass(frozen=True)
class RationalResult:
    value: Fraction
    attempts: list[dict[str, str]]


def graph_kernel_from_basis(
    basis: list[list[mpc]],
    canonical_curve: Any,
    monomials: list[tuple[int, ...]],
    length: int,
) -> dict[str, Any]:
    sections, lam = p200_map.fixed_sections(
        basis, canonical_curve, monomials, length
    )
    _, _, kappas = numerical.triangle_vertices_and_kappas()
    phi = numerical.phi_series(kappas[0], length)
    phi_sections = [
        numerical.series_mul(phi, section, length) for section in sections
    ]
    columns = sections + [[-value for value in series] for series in phi_sections]
    kernel, rank = canonical_curve.nullspace(
        columns, length, expect=4, tol=mpf("1e-82")
    )
    residual = max(
        numerical.kernel_residual(columns, vector, length) for vector in kernel
    )
    if rank != 50:
        raise ValueError(f"graph rank {rank}, expected 50")
    indices, q_values = canonical_curve.quadric(basis, n_use=length)
    invariants = {
        name: mpc(value)
        for name, value in canonical_curve.invariants(indices, q_values)
    }
    return {
        "kernel": kernel,
        "rank": rank,
        "residual": residual,
        "lambda": lam,
        "invariants": invariants,
    }


def solve_half_gauge(
    kernel: list[list[mpc]],
    half_offset: int,
    monomials: list[tuple[int, ...]],
) -> list[mpc]:
    gauge_rows = tuple(monomials.index(monomial) for monomial in GAUGE_MONOMIALS)
    gauge_unit_row = monomials.index(GAUGE_UNIT_MONOMIAL)
    absolute_rows = [half_offset + row for row in gauge_rows]
    system = matrix(
        [[kernel[column][row] for column in range(4)] for row in absolute_rows]
    )
    target = matrix([1 if row == gauge_unit_row else 0 for row in gauge_rows])
    weights = lu_solve(system, target)
    return [
        sum(kernel[column][row] * weights[column] for column in range(4))
        for row in range(54)
    ]


def common_monomial_basis(
    bases: dict[int, list[list[mpc]]],
    canonical_curve: Any,
    length: int,
) -> tuple[list[tuple[int, ...]], dict[str, Any]]:
    """Greedily maximize worst conditioning across all six embeddings."""

    all_monomials = list(numerical.MONOMIALS_5)
    evaluations: dict[int, np.ndarray] = {}
    for label, basis in bases.items():
        sections, _ = p200_map.fixed_sections(
            basis, canonical_curve, all_monomials, length
        )
        array = np.array(
            [
                [complex(sections[column][degree]) for column in range(56)]
                for degree in range(length)
            ],
            dtype=np.complex128,
        )
        norms = np.linalg.norm(array, axis=0)
        if np.any(norms == 0) or not np.all(np.isfinite(norms)):
            raise ValueError(f"label {label} has invalid monomial column norms")
        evaluations[label] = array / norms

    selected = [all_monomials.index(monomial) for monomial in GAUGE_MONOMIALS]
    trace = []
    while len(selected) < 27:
        best_index: int | None = None
        best_score = -1.0
        for candidate in range(56):
            if candidate in selected:
                continue
            trial = selected + [candidate]
            score = min(
                float(
                    np.linalg.svd(
                        evaluation[:, trial], compute_uv=False
                    )[-1]
                )
                for evaluation in evaluations.values()
            )
            if score > best_score:
                best_score = score
                best_index = candidate
        if best_index is None:
            raise ValueError("common monomial basis selection stalled")
        selected.append(best_index)
        trace.append(
            {
                "dimension": len(selected),
                "added_monomial": str(all_monomials[best_index]),
                "worst_smallest_singular_value": best_score,
            }
        )
    final_singular = {
        str(label): float(
            np.linalg.svd(evaluation[:, selected], compute_uv=False)[-1]
        )
        for label, evaluation in evaluations.items()
    }
    if min(final_singular.values()) <= 1e-10:
        raise ValueError(f"common monomial basis remains ill-conditioned: {final_singular}")
    monomials = [all_monomials[index] for index in selected]
    return monomials, {
        "selection_rule": (
            "start with the four p200 denominator-gauge monomials, then "
            "greedily add the monomial maximizing the worst smallest "
            "singular value across the six column-normalized p100 evaluation matrices"
        ),
        "trace": trace,
        "final_smallest_singular_values": final_singular,
    }


def evaluate(coefficients: list[Fraction], alpha: mpc) -> mpc:
    result = mpc(0)
    for coefficient in reversed(coefficients):
        result = result * alpha + mpf(coefficient.numerator) / coefficient.denominator
    return result


def recognize_rational(value: mpc) -> RationalResult:
    attempts: list[dict[str, str]] = []
    candidates: list[Fraction] = []
    for tolerance in PSLQ_TOLERANCES:
        if fabs(value.imag) >= tolerance:
            raise ValueError(
                f"interpolated coefficient has imaginary part {value.imag} "
                f"at tolerance {tolerance}"
            )
        nearest = int(mp.nint(value.real))
        if fabs(value.real - nearest) < tolerance:
            candidate = Fraction(nearest)
        else:
            relation = pslq(
                [mpf(1), value.real],
                maxcoeff=PSLQ_BOUND,
                maxsteps=2_000_000,
                tol=tolerance,
            )
            if relation is None or relation[1] == 0:
                raise ValueError(
                    f"no rational relation at tolerance {tolerance}: {value}"
                )
            candidate = Fraction(-int(relation[0]), int(relation[1]))
        residual = fabs(
            value - mpf(candidate.numerator) / candidate.denominator
        )
        attempts.append(
            {
                "tolerance": str(tolerance),
                "candidate": str(candidate),
                "residual": str(residual),
            }
        )
        candidates.append(candidate)
    if candidates[0] != candidates[1]:
        raise ValueError(f"rational relation changed with tolerance: {attempts}")
    return RationalResult(candidates[0], attempts)


def serialize_fractions(values: list[Fraction]) -> list[str]:
    return [str(value) for value in values]


def parse_alpha(path: Path) -> mpc:
    payload = json.loads(path.read_text())
    return p200_map.parse_complex_text(payload["field"]["label1_compact_root"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pipeline-dir",
        type=Path,
        default=Path("outputs/hm-degeneration-20260819/phaseB-open"),
    )
    parser.add_argument("--length", type=int, default=100)
    parser.add_argument(
        "--p200-label1",
        type=Path,
        default=Path("/private/tmp/m23-p200-recovered/label1-out-p200.m"),
    )
    parser.add_argument(
        "--p200-label2",
        type=Path,
        default=Path("/private/tmp/m23-p200-recovered/label2-out-p200.m"),
    )
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
            "p100-map-interpolation.json"
        ),
    )
    args = parser.parse_args()
    if args.length < 80 or args.length > 120:
        raise ValueError("p100 length must lie between 80 and 120")

    mp.dps = 110
    pipeline_dir = args.pipeline_dir.resolve()
    sys.path.insert(0, str(Path("outputs/m23-proof-closeout-20260902").resolve()))
    sys.path.insert(0, str(pipeline_dir))
    parse_series = importlib.import_module("parse_series")
    canonical_curve = importlib.import_module("canonical_curve")
    petri_probe = importlib.import_module("petri_exactify_probe")
    petri_algebra = importlib.import_module("petri_quadric_exactify")
    parse_series.mp.dps = 110
    canonical_curve.mp.dps = 110

    bases = {
        label: parse_series.get_basis(
            label,
            nuse=args.length,
            path=str(pipeline_dir / f"label{label}-out-p100.m"),
        )
        for label in NONHARMONIC
    }
    monomials, monomial_basis_receipt = common_monomial_basis(
        bases, canonical_curve, args.length
    )
    print(
        "common monomial basis worst singular",
        min(monomial_basis_receipt["final_smallest_singular_values"].values()),
        flush=True,
    )

    graphs: dict[int, dict[str, Any]] = {}
    primary: dict[int, list[mpc]] = {}
    invariants: dict[int, dict[str, mpc]] = {}
    for label in NONHARMONIC:
        graph = graph_kernel_from_basis(
            bases[label], canonical_curve, monomials, args.length
        )
        graphs[label] = graph
        invariants[label] = graph["invariants"]
        primary[label] = solve_half_gauge(graph["kernel"], 27, monomials)
        print(
            "p100 label",
            label,
            "rank",
            graph["rank"],
            "residual",
            graph["residual"],
            flush=True,
        )

    # Compact-root matching uses the already certified descended Petri value.
    j1_values = {label: invariants[label]["J1"] for label in NONHARMONIC}
    tau_j1 = {
        label: j1_values[KAPPA[label]].conjugate() for label in NONHARMONIC
    }
    theta = {
        label: j1_values[label]
        + tau_j1[label]
        + j1_values[label] * tau_j1[label]
        for label in NONHARMONIC
    }
    alphas, alpha_errors = petri_algebra.match_compact_roots(theta)
    xs = [alphas[label] for label in NONHARMONIC]
    s = mpc(0, mp.sqrt(23))

    tau_targets: dict[int, list[mpc]] = {}
    for label in NONHARMONIC:
        partner = KAPPA[label]
        tau_targets[label] = [value.conjugate() for value in primary[partner]]

    exact_coefficients: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    embedding_errors: list[dict[str, Any]] = []
    for row in range(54):
        a_values = [
            (primary[label][row] + tau_targets[label][row]) / 2
            for label in NONHARMONIC
        ]
        b_values = [
            (primary[label][row] - tau_targets[label][row]) / (2 * s)
            for label in NONHARMONIC
        ]
        a_numeric = petri_probe.interpolate(xs, a_values)
        b_numeric = petri_probe.interpolate(xs, b_values)
        try:
            a_results = [recognize_rational(value) for value in a_numeric]
            b_results = [recognize_rational(value) for value in b_numeric]
            a = [result.value for result in a_results]
            b = [result.value for result in b_results]
            row_errors = []
            for label in NONHARMONIC:
                value = evaluate(a, alphas[label]) + s * evaluate(b, alphas[label])
                tau_value = evaluate(a, alphas[label]) - s * evaluate(
                    b, alphas[label]
                )
                error = fabs(value - primary[label][row])
                tau_error = fabs(tau_value - tau_targets[label][row])
                if max(error, tau_error) >= EMBEDDING_TOL:
                    raise ValueError(
                        f"embedding error label {label}: c={error}, tau={tau_error}"
                    )
                row_errors.append(
                    {
                        "label": label,
                        "c": str(error),
                        "tau_c": str(tau_error),
                    }
                )
            embedding_errors.extend(
                {"graph_row": row, **entry} for entry in row_errors
            )
            exact_coefficients.append(
                {
                    "graph_row": row,
                    "half": "N" if row < 27 else "D",
                    "section_row": row if row < 27 else row - 27,
                    "monomial": str(monomials[row if row < 27 else row - 27]),
                    "A_coefficients_ascending": serialize_fractions(a),
                    "B_coefficients_ascending": serialize_fractions(b),
                    "recognition": {
                        "A": [result.attempts for result in a_results],
                        "B": [result.attempts for result in b_results],
                    },
                }
            )
            print("interpolated row", row, flush=True)
        except Exception as error:  # preserve a bounded negative, then continue
            failures.append(
                {
                    "graph_row": row,
                    "half": "N" if row < 27 else "D",
                    "monomial": str(monomials[row if row < 27 else row - 27]),
                    "error": str(error),
                    "A_numeric_coefficients": [mp.nstr(value, 105) for value in a_numeric],
                    "B_numeric_coefficients": [mp.nstr(value, 105) for value in b_numeric],
                }
            )
            print("FAILED row", row, error, flush=True)

    p200_checks: dict[str, Any] | None = None
    if not failures:
        # Out-of-sample validation from the independent p200 files.
        mp.dps = 220
        p200 = importlib.import_module("petri_p200_field_recognition")
        p200.mp.dps = 220
        canonical_curve.mp.dps = 220
        check_1 = p200_map.graph_kernel(
            p200,
            canonical_curve,
            pipeline_dir,
            args.p200_label1.resolve(),
            120,
            monomials,
        )
        check_2 = p200_map.graph_kernel(
            p200,
            canonical_curve,
            pipeline_dir,
            args.p200_label2.resolve(),
            120,
            monomials,
        )
        p200_primary = solve_half_gauge(check_1["kernel"], 27, monomials)
        p200_partner = solve_half_gauge(check_2["kernel"], 27, monomials)
        p200_tau = [value.conjugate() for value in p200_partner]
        alpha1 = parse_alpha(args.field_receipt)
        errors = []
        for entry in exact_coefficients:
            row = entry["graph_row"]
            a = [Fraction(value) for value in entry["A_coefficients_ascending"]]
            b = [Fraction(value) for value in entry["B_coefficients_ascending"]]
            value = evaluate(a, alpha1) + s * evaluate(b, alpha1)
            tau_value = evaluate(a, alpha1) - s * evaluate(b, alpha1)
            errors.append(
                {
                    "graph_row": row,
                    "label1": str(fabs(value - p200_primary[row])),
                    "label2_scalar_tau": str(fabs(tau_value - p200_tau[row])),
                }
            )
        max_error = max(
            mpf(entry[key])
            for entry in errors
            for key in ("label1", "label2_scalar_tau")
        )
        if max_error >= mpf("1e-145"):
            raise ValueError(f"p200 out-of-sample error {max_error}")
        p200_checks = {
            "jet_length": 120,
            "max_error": str(max_error),
            "per_coefficient": errors,
        }

    max_embedding_error = (
        max(
            mpf(entry[key])
            for entry in embedding_errors
            for key in ("c", "tau_c")
        )
        if embedding_errors
        else None
    )
    payload = {
        "status": (
            "conditional_exact_map_candidate" if not failures else "failed_closed"
        ),
        "construction": {
            "p100_labels": list(NONHARMONIC),
            "series_length": args.length,
            "semilinear_pairing": (
                "c_l uses D-gauge; scalar tau(c_l) is the complex conjugate "
                "of the D-gauged c_kappa(l), with no graph-half swap"
            ),
            "separate_target_action": (
                "tau(Graph(phi))=Graph(1/phi) requires a separate N/D block "
                "swap and target/source normalization check; it is not used here"
            ),
            "degree_five_monomial_basis": [str(value) for value in monomials],
            "monomial_basis_selection": monomial_basis_receipt,
            "gauge_monomials": [str(value) for value in GAUGE_MONOMIALS],
            "gauge_unit_monomial": str(GAUGE_UNIT_MONOMIAL),
        },
        "compact_root_matching_errors": {
            str(label): str(error)
            for label, error in zip(NONHARMONIC, alpha_errors, strict=True)
        },
        "graph_checks": {
            str(label): {
                "rank": graphs[label]["rank"],
                "kernel_dimension": 4,
                "residual": str(graphs[label]["residual"]),
                "lambda": mp.nstr(graphs[label]["lambda"], 105),
            }
            for label in NONHARMONIC
        },
        "rational_recognition": {
            "tolerances": [str(value) for value in PSLQ_TOLERANCES],
            "coefficient_bound": str(PSLQ_BOUND),
            "embedding_tolerance": str(EMBEDDING_TOL),
            "accepted_rows": len(exact_coefficients),
            "failed_rows": len(failures),
            "max_p100_embedding_error": (
                str(max_embedding_error) if max_embedding_error is not None else None
            ),
        },
        "map": {"coefficients": exact_coefficients},
        "failures": failures,
        "p200_out_of_sample_check": p200_checks,
        "boundary": {
            "proves": (
                "if status is conditional_exact_map_candidate: all 54 graph "
                "coefficients interpolate to tolerance-stable rational A,B "
                "power-basis coefficients across all six p100 embeddings and "
                "also match the independent p200 label1/label2 pair"
            ),
            "does_not_prove": (
                "that the analytic solutions are exact covers, or the exact "
                "smoothness, degree, three-fiber passport, M23 monodromy, and "
                "good reduction at 13; it also does not close the separate "
                "N/D block-swap target-normalization compatibility gate"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)
    print("status", payload["status"], "failed rows", len(failures))


if __name__ == "__main__":
    main()
