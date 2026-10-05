#!/usr/bin/env python3
"""Reduce the conditional Petri model mod 13 and probe 23-fold osculation.

This is deliberately a bounded discovery probe.  It reduces the exact
quadric/cubic certificates at the prime ``(13, alpha^6-alpha^5-...+4,
s-4)`` and works in ``k = F_13[alpha]``.  At the marked point
``P+ = [1:0:0:0]`` it constructs the canonical local parameter ``t=x1/x0``
and solves the curve equations for ``x2/x0`` and ``x3/x0`` through order 26.

The decisive diagnostic is the kernel of the 0..22 jet map on H^0(4K).
For an actual degree-23 cover normalized by

    beta = (phi-s)/(phi+s),  div(beta) = 23(P+ - P-),

there must be a nonzero quartic section with divisor ``23 P+ + R``.  Since
the canonical curve is a (2,3) complete intersection, H^0(4K) has dimension
21.  The probe computes this kernel exactly over k, after quotienting out
the 14-dimensional degree-four part of the canonical ideal.

No numerical approximation is used after loading the rational certificates.
Existence of the expected one-dimensional osculating kernel is necessary and
strong evidence for the candidate curve, but is not by itself a cover
certificate: the opposite ramification point and the map still need to be
constructed and the complete passport verified.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from math import lcm
from pathlib import Path
from typing import Any

import galois
import numpy as np


P = 13
G_COEFFICIENTS_DESCENDING = [1, -1, -3, -3, 1, 5, 4]
QUADRATIC_MONOMIALS = tuple(itertools.combinations_with_replacement(range(4), 2))
CUBIC_MONOMIALS = tuple(itertools.combinations_with_replacement(range(4), 3))
QUARTIC_MONOMIALS = tuple(itertools.combinations_with_replacement(range(4), 4))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def make_field() -> tuple[type, Any]:
    prime = galois.GF(P)
    polynomial = galois.Poly(
        [coefficient % P for coefficient in G_COEFFICIENTS_DESCENDING],
        field=prime,
    )
    if not polynomial.is_irreducible():
        raise AssertionError("the compact sextic is reducible modulo 13")
    field = galois.GF(P**6, irreducible_poly=polynomial)
    alpha = field(P)  # coefficient vector [0,0,0,0,1,0], i.e. x mod g
    value = field(0)
    for coefficient in G_COEFFICIENTS_DESCENDING:
        value = value * alpha + field(coefficient % P)
    if value != 0:
        raise AssertionError("chosen field element is not the compact root")
    return field, alpha


def valuation(integer: int) -> int:
    if integer == 0:
        return 10**9
    result = 0
    while integer % P == 0:
        integer //= P
        result += 1
    return result


def hensel_sqrt_minus_23(residue: int, digits: int = 100) -> int:
    root = residue
    modulus = P
    for _ in range(1, digits):
        next_modulus = modulus * P
        lifts = [
            root + digit * modulus
            for digit in range(P)
            if ((root + digit * modulus) ** 2 + 23) % next_modulus == 0
        ]
        if len(lifts) != 1:
            raise AssertionError("sqrt(-23) Hensel lift is not unique")
        root = lifts[0]
        modulus = next_modulus
    return root


def quadratic_residue_mod_13(
    a_text: str, b_text: str, root: int, scale: int, field: type
) -> Any:
    """Return the residue of 13^scale * (a + b*s) at s -> root.

    The rational A/B components often have denominators divisible by 13 even
    when their sum is integral at the chosen split prime.  Reduction must
    therefore happen *after* combining them in Q_13, not componentwise.
    """

    a = Fraction(a_text)
    b = Fraction(b_text)
    denominator = lcm(a.denominator, b.denominator)
    a_integer = a.numerator * (denominator // a.denominator)
    b_integer = b.numerator * (denominator // b.denominator)
    numerator = a_integer + b_integer * root
    numerator_valuation = valuation(numerator)
    denominator_valuation = valuation(denominator)
    total_valuation = scale + numerator_valuation - denominator_valuation
    if total_valuation < 0:
        raise ZeroDivisionError(
            "coefficient remains nonintegral after equation scaling: "
            f"v13={total_valuation}, a={a_text}, b={b_text}"
        )
    if total_valuation > 0 or numerator == 0:
        return field(0)
    numerator_unit = numerator // (P**numerator_valuation)
    denominator_unit = denominator // (P**denominator_valuation)
    return field(numerator_unit % P) / field(denominator_unit % P)


def coefficient_mod(
    payload: dict[str, Any], field: type, alpha: Any, root: int, scale: int
) -> Any:
    result = field(0)
    power = field(1)
    for a_text, b_text in zip(
        payload["A_coefficients_ascending"],
        payload["B_coefficients_ascending"],
        strict=True,
    ):
        result += quadratic_residue_mod_13(
            a_text, b_text, root, scale, field
        ) * power
        power *= alpha
    return result


def parse_monomial(text: str) -> tuple[int, ...]:
    return tuple(int(part.strip()) for part in text.strip("()").split(","))


def load_model(
    quadric_path: Path, cubic_path: Path, field: type, alpha: Any, s: int
) -> tuple[dict[tuple[int, ...], Any], dict[tuple[int, ...], Any]]:
    quadric_payload = json.loads(quadric_path.read_text())
    cubic_payload = json.loads(cubic_path.read_text())
    root = hensel_sqrt_minus_23(s)
    quadric = {
        parse_monomial(key): coefficient_mod(value, field, alpha, root, 0)
        for key, value in quadric_payload["normalized_quadric"]["coefficients"].items()
    }
    cubic = {
        # This chosen cubic gauge has minimum v_13 = -1 at s -> 4.
        # Multiplication by 13 is a harmless projective equation scaling and
        # gives its primitive integral reduction.
        parse_monomial(key): coefficient_mod(value, field, alpha, root, 1)
        for key, value in cubic_payload["cubic_coefficients"].items()
    }
    if set(quadric) != set(QUADRATIC_MONOMIALS):
        raise AssertionError("quadric monomial support is incomplete")
    if set(cubic) != set(CUBIC_MONOMIALS):
        raise AssertionError("cubic monomial support is incomplete")
    return quadric, cubic


def series_mul(left: Any, right: Any) -> Any:
    field = type(left)
    length = left.size
    result = field.Zeros(length)
    for degree in range(length):
        result[degree] = sum(
            (left[index] * right[degree - index] for index in range(degree + 1)),
            start=field(0),
        )
    return result


def series_monomial(coordinates: list[Any], monomial: tuple[int, ...]) -> Any:
    field = type(coordinates[0])
    result = field.Zeros(coordinates[0].size)
    result[0] = 1
    for variable in monomial:
        result = series_mul(result, coordinates[variable])
    return result


def polynomial_series(
    polynomial: dict[tuple[int, ...], Any], coordinates: list[Any]
) -> Any:
    field = type(coordinates[0])
    result = field.Zeros(coordinates[0].size)
    for monomial, coefficient in polynomial.items():
        result += coefficient * series_monomial(coordinates, monomial)
    return result


def local_coordinates(
    quadric: dict[tuple[int, ...], Any],
    cubic: dict[tuple[int, ...], Any],
    field: type,
    order: int,
) -> list[Any]:
    one = field.Zeros(order)
    one[0] = 1
    parameter = field.Zeros(order)
    parameter[1] = 1
    z = field.Zeros(order)
    w = field.Zeros(order)
    jacobian = field(
        [
            [quadric[(0, 2)], quadric[(0, 3)]],
            [cubic[(0, 0, 2)], cubic[(0, 0, 3)]],
        ]
    )
    if np.linalg.det(jacobian) == 0:
        raise AssertionError("P+ is singular or x1 is not a local parameter")
    inverse = np.linalg.inv(jacobian)
    coordinates = [one, parameter, z, w]
    for degree in range(1, order):
        q_residual = polynomial_series(quadric, coordinates)[degree]
        c_residual = polynomial_series(cubic, coordinates)[degree]
        correction = -(inverse @ field([q_residual, c_residual]))
        z[degree] = correction[0]
        w[degree] = correction[1]
    if np.any(polynomial_series(quadric, coordinates) != 0):
        raise AssertionError("quadric local expansion did not converge")
    if np.any(polynomial_series(cubic, coordinates) != 0):
        raise AssertionError("cubic local expansion did not converge")
    return coordinates


def vector_for_product(
    left: dict[tuple[int, ...], Any],
    right_monomial: tuple[int, ...],
    index: dict[tuple[int, ...], int],
    field: type,
) -> Any:
    result = field.Zeros(len(index))
    for monomial, coefficient in left.items():
        product = tuple(sorted((*monomial, *right_monomial)))
        result[index[product]] += coefficient
    return result


def rank(matrix: Any) -> int:
    return int(matrix.row_space().shape[0])


def degree_four_ideal(
    quadric: dict[tuple[int, ...], Any],
    cubic: dict[tuple[int, ...], Any],
    field: type,
) -> Any:
    index = {monomial: position for position, monomial in enumerate(QUARTIC_MONOMIALS)}
    rows = [
        vector_for_product(quadric, monomial, index, field)
        for monomial in QUADRATIC_MONOMIALS
    ]
    rows.extend(
        vector_for_product(cubic, (variable,), index, field) for variable in range(4)
    )
    ideal = field(rows)
    if rank(ideal) != 14:
        raise AssertionError(f"degree-four ideal rank is {rank(ideal)}, expected 14")
    return ideal.row_space()


def serialize_field_element(value: Any) -> list[int]:
    return [int(coefficient) for coefficient in value.vector()]


def serialize_polynomial(vector: Any) -> dict[str, list[int]]:
    return {
        str(monomial): serialize_field_element(vector[index])
        for index, monomial in enumerate(QUARTIC_MONOMIALS)
        if vector[index] != 0
    }


def first_nonzero(series: Any) -> int | None:
    indices = np.flatnonzero(series != 0)
    return int(indices[0]) if indices.size else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quadric",
        type=Path,
        default=Path(
            "outputs/m23-proof-closeout-20260902/artifacts/"
            "petri-p200-quadric-certificate.json"
        ),
    )
    parser.add_argument(
        "--cubic",
        type=Path,
        default=Path(
            "outputs/m23-proof-closeout-20260902/artifacts/"
            "petri-p200-cubic-certificate.json"
        ),
    )
    parser.add_argument("--s", type=int, choices=(4, 9), default=4)
    parser.add_argument("--order", type=int, default=27)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "outputs/m23-proof-closeout-20260902/agents/petri-f13/"
            "osculation-probe.json"
        ),
    )
    args = parser.parse_args()
    if args.order < 25:
        raise ValueError("order must be at least 25")

    field, alpha = make_field()
    quadric, cubic = load_model(args.quadric, args.cubic, field, alpha, args.s)
    coordinates = local_coordinates(quadric, cubic, field, args.order)
    monomial_series = [
        series_monomial(coordinates, monomial) for monomial in QUARTIC_MONOMIALS
    ]
    jet_matrix = field(
        [
            [series[degree] for series in monomial_series]
            for degree in range(23)
        ]
    )
    kernel = jet_matrix.null_space()
    ideal = degree_four_ideal(quadric, cubic, field)
    combined_rank = rank(field(np.vstack([ideal, kernel])))
    quotient_kernel_dimension = combined_rank - rank(ideal)
    if rank(kernel) != kernel.shape[0]:
        raise AssertionError("null-space basis is rank deficient")
    if rank(jet_matrix) + kernel.shape[0] != len(QUARTIC_MONOMIALS):
        raise AssertionError("rank-nullity failed")
    if rank(field(np.vstack([ideal, kernel]))) > 15:
        raise AssertionError("unexpectedly large osculating quotient kernel")

    osculating = None
    span = ideal
    for vector in kernel:
        trial = field(np.vstack([span, vector]))
        if rank(trial) > rank(span):
            osculating = vector
            break
    if quotient_kernel_dimension and osculating is None:
        raise AssertionError("failed to extract osculating section")

    section_series = None
    section_order = None
    if osculating is not None:
        section_series = sum(
            (
                osculating[index] * monomial_series[index]
                for index in range(len(QUARTIC_MONOMIALS))
            ),
            start=field.Zeros(args.order),
        )
        section_order = first_nonzero(section_series)

    payload = {
        "status": "exact_finite_field_probe",
        "inputs": {
            "quadric": str(args.quadric.resolve()),
            "quadric_sha256": sha256(args.quadric),
            "cubic": str(args.cubic.resolve()),
            "cubic_sha256": sha256(args.cubic),
        },
        "field": {
            "characteristic": P,
            "degree": 6,
            "order": P**6,
            "compact_polynomial_coefficients_descending_mod_13": [
                coefficient % P for coefficient in G_COEFFICIENTS_DESCENDING
            ],
            "compact_polynomial_irreducible": True,
            "alpha_vector_high_to_low": serialize_field_element(alpha),
            "sqrt_minus_23": args.s,
            "sqrt_check": (args.s * args.s + 23) % P == 0,
        },
        "reduced_model": {
            "ambient": "P^3 over F_(13^6)",
            "curve_definition": "Q=0 and C=0",
            "quadric_coefficients": {
                str(monomial): serialize_field_element(value)
                for monomial, value in quadric.items()
            },
            "cubic_coefficients": {
                str(monomial): serialize_field_element(value)
                for monomial, value in cubic.items()
            },
            "integral_scaling": {
                "quadric": "no power of 13",
                "cubic": "multiply the characteristic-zero chosen gauge by 13 before reduction",
                "prime": "(13,s-4)",
            },
        },
        "marked_point": {
            "projective_coordinates": [[1], [0], [0], [0]],
            "local_parameter": "t=x1/x0",
            "jacobian_nonzero": True,
            "curve_equations_vanish_through_order": args.order - 1,
            "x2_series_coefficients_ascending": [
                serialize_field_element(value) for value in coordinates[2]
            ],
            "x3_series_coefficients_ascending": [
                serialize_field_element(value) for value in coordinates[3]
            ],
        },
        "degree_four_osculation": {
            "ambient_quartic_monomials": len(QUARTIC_MONOMIALS),
            "degree_four_ideal_rank": rank(ideal),
            "h0_4K_dimension": len(QUARTIC_MONOMIALS) - rank(ideal),
            "jet_orders_imposed": [0, 22],
            "jet_matrix_rank": rank(jet_matrix),
            "ambient_kernel_dimension": kernel.shape[0],
            "quotient_kernel_dimension": quotient_kernel_dimension,
            "osculating_section": (
                serialize_polynomial(osculating) if osculating is not None else None
            ),
            "osculating_section_vanishing_order": section_order,
            "section_series_coefficients_23_onward": (
                [
                    serialize_field_element(section_series[degree])
                    for degree in range(23, args.order)
                ]
                if section_series is not None
                else None
            ),
        },
        "boundary": {
            "proves": (
                "the exact reduced candidate curve has the reported 4K jet "
                "kernel at the marked point over F_(13^6)"
            ),
            "does_not_prove": (
                "that the candidate is an actual Hurwitz cover; the opposite "
                "23-ramification point, degree-23 map, full tame passport, and "
                "M23 monodromy still require exact construction"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)
    print(
        "ranks:",
        f"ideal={rank(ideal)}",
        f"jet={rank(jet_matrix)}",
        f"kernel={kernel.shape[0]}",
        f"quotient-kernel={quotient_kernel_dimension}",
        f"vanishing-order={section_order}",
    )


if __name__ == "__main__":
    main()
