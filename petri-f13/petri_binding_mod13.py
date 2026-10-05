#!/usr/bin/env python3
"""Bind the exact characteristic-13 cover to the descended Petri value.

For each totally ramified point this verifier constructs the distinguished
local parameter

    t_0^23 = f,       t_inf^23 = 1/f.

Because gcd(23, 13^6-1)=1, the leading 23rd root is unique in the residue
field.  The four canonical coordinate series are composed with this adapted
parameter, echelonized through order three, and their unique quadratic Petri
relation is recovered from exact 0..12 jets.  The marked invariant is

    J = q_13 / q_22.

The unmarked descended value U1=J0+Jinf+J0*Jinf is then compared exactly with
the previously certified primitive sextic generator modulo 13.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import galois
import numpy as np

import exact_f13_map as core


HERE = Path(__file__).resolve().parent
JET_COUNT = 13


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def deserialize_element(vector: list[int], field: type, alpha: Any) -> Any:
    value = field(0)
    for coefficient in vector:
        value = value * alpha + field(coefficient % core.PRIME)
    return value


def deserialize_form(
    payload: dict[str, list[int]], degree: int, field: type, alpha: Any
) -> Any:
    index = {
        monomial: position for position, monomial in enumerate(core.MONOMIALS[degree])
    }
    result = field.Zeros(len(index))
    for monomial_text, vector in payload.items():
        monomial = core.parse_monomial(monomial_text)
        result[index[monomial]] = deserialize_element(vector, field, alpha)
    return result


def series_divide(numerator: Any, denominator: Any) -> Any:
    field = type(numerator)
    if denominator[0] == 0:
        raise ZeroDivisionError("series denominator is not a unit")
    result = field.Zeros(numerator.size)
    for degree in range(numerator.size):
        accumulated = field(0)
        for index in range(degree):
            accumulated = accumulated + result[index] * denominator[degree - index]
        result[degree] = (numerator[degree] - accumulated) / denominator[0]
    if np.any(core.series_mul(result, denominator) != numerator):
        raise AssertionError("formal series division failed")
    return result


def series_power(series: Any, exponent: int) -> Any:
    field = type(series)
    result = field.Zeros(series.size)
    result[0] = 1
    base = series.copy()
    power = exponent
    while power:
        if power & 1:
            result = core.series_mul(result, base)
        power >>= 1
        if power:
            base = core.series_mul(base, base)
    return result


def unit_nth_root(series: Any, exponent: int) -> Any:
    field = type(series)
    if series[0] == 0:
        raise ValueError("unit root requested for nonunit series")
    group_order = field.order - 1
    if np.gcd(exponent, group_order) != 1:
        raise ValueError("leading nth root is not unique in the coefficient field")
    root_exponent = pow(exponent, -1, group_order)
    result = field.Zeros(series.size)
    result[0] = series[0] ** root_exponent
    if result[0] ** exponent != series[0]:
        raise AssertionError("finite-field leading root failed")
    derivative = field(exponent % field.characteristic) * result[0] ** (exponent - 1)
    if derivative == 0:
        raise ValueError("nth-root Hensel derivative vanishes")
    for degree in range(1, series.size):
        known = series_power(result, exponent)[degree]
        result[degree] = (series[degree] - known) / derivative
    if np.any(series_power(result, exponent) != series):
        raise AssertionError("formal unit nth root failed")
    return result


def compose(outer: Any, inner: Any) -> Any:
    field = type(outer)
    result = field.Zeros(outer.size)
    power = field.Zeros(outer.size)
    power[0] = 1
    for coefficient in outer:
        result += coefficient * power
        power = core.series_mul(power, inner)
    return result


def compositional_inverse(series: Any) -> Any:
    field = type(series)
    if series[0] != 0 or series[1] == 0:
        raise ValueError("series is not composition-invertible at the origin")
    result = field.Zeros(series.size)
    result[1] = field(1) / series[1]
    for degree in range(2, series.size):
        known = compose(series, result)[degree]
        result[degree] = -known / series[1]
    identity = compose(series, result)
    expected = field.Zeros(series.size)
    expected[1] = 1
    if np.any(identity != expected):
        raise AssertionError("formal compositional inverse failed")
    return result


def adapted_parameter(map_series: Any, field: type, jet_count: int) -> dict[str, Any]:
    order = map_series.size
    first = core.first_nonzero(map_series)
    if first != 23:
        raise ValueError(f"map series has order {first}, expected 23")
    unit = field.Zeros(jet_count)
    available = min(jet_count, order - 23)
    unit[:available] = map_series[23 : 23 + available]
    root = unit_nth_root(unit, 23)
    t_of_u = field.Zeros(jet_count + 1)
    t_of_u[1:] = root
    u_of_t = compositional_inverse(t_of_u)
    if np.any(series_power(t_of_u, 23)[23:] != map_series[23 : 23 + jet_count - 22]):
        # This direct check covers the coefficients available in the shorter
        # adapted series; the root identity above is the primary certificate.
        raise AssertionError("adapted parameter does not reproduce the map jet")
    return {
        "unit": unit,
        "unit_root": root,
        "t_of_u": t_of_u,
        "u_of_t": u_of_t,
    }


def derivative_polynomial(
    polynomial: dict[tuple[int, ...], Any], variable: int, field: type
) -> dict[tuple[int, ...], Any]:
    result: dict[tuple[int, ...], Any] = {}
    for monomial, coefficient in polynomial.items():
        multiplicity = monomial.count(variable)
        if multiplicity == 0:
            continue
        reduced = list(monomial)
        reduced.remove(variable)
        key = tuple(reduced)
        result[key] = result.get(key, field(0)) + field(multiplicity) * coefficient
    return result


def marked_invariant(
    expansion: core.LocalExpansion,
    map_series: Any,
    quadric: dict[tuple[int, ...], Any],
    cubic: dict[tuple[int, ...], Any],
    field: type,
) -> tuple[Any, dict[str, Any]]:
    adapted = adapted_parameter(map_series, field, JET_COUNT)
    u_of_t = adapted["u_of_t"]
    projective_coordinate_series = []
    for coordinate in expansion.coordinates:
        short = field.Zeros(JET_COUNT + 1)
        short[:] = coordinate[: JET_COUNT + 1]
        projective_coordinate_series.append(compose(short, u_of_t))

    # Adjunction for a (2,3) complete intersection.  In the affine chart and
    # with u=x_free, the canonical differential represented by x_i has
    # coefficient x_i/J with respect to du, where J is the implicit-function
    # Jacobian determinant.  Relative to dt this becomes x_i*(du/dt)/J.
    left, right = expansion.solved_variables
    q_left = derivative_polynomial(quadric, left, field)
    q_right = derivative_polynomial(quadric, right, field)
    c_left = derivative_polynomial(cubic, left, field)
    c_right = derivative_polynomial(cubic, right, field)
    jacobian_u = (
        core.series_mul(
            core.polynomial_series(q_left, expansion.coordinates),
            core.polynomial_series(c_right, expansion.coordinates),
        )
        - core.series_mul(
            core.polynomial_series(q_right, expansion.coordinates),
            core.polynomial_series(c_left, expansion.coordinates),
        )
    )
    jacobian_short = field.Zeros(JET_COUNT + 1)
    jacobian_short[:] = jacobian_u[: JET_COUNT + 1]
    jacobian_t = compose(jacobian_short, u_of_t)
    du_dt = field.Zeros(JET_COUNT + 1)
    for degree in range(JET_COUNT):
        du_dt[degree] = field((degree + 1) % field.characteristic) * u_of_t[degree + 1]
    common_differential_factor = series_divide(du_dt, jacobian_t)
    coordinate_series = [
        core.series_mul(common_differential_factor, coordinate)
        for coordinate in projective_coordinate_series
    ]

    leading = field(
        [
            [coordinate_series[row][degree] for degree in range(4)]
            for row in range(4)
        ]
    )
    if np.linalg.det(leading) == 0:
        raise ValueError("canonical differential leading-jet matrix is singular")
    transform = np.linalg.inv(leading)
    raw = field(
        [
            [coordinate_series[row][degree] for degree in range(JET_COUNT)]
            for row in range(4)
        ]
    )
    echelon = transform @ raw
    if np.any(echelon[:, :4] != field.Identity(4)):
        raise AssertionError("canonical differential echelonization failed")

    indices = tuple(itertools.combinations_with_replacement(range(4), 2))
    products = [core.series_mul(echelon[left], echelon[right]) for left, right in indices]
    jet_matrix = field(
        [[product[degree] for product in products] for degree in range(JET_COUNT)]
    )
    kernel = jet_matrix.null_space()
    if kernel.shape[0] != 1 or core.rank(jet_matrix) != 9:
        raise ValueError(
            f"quadratic relation has rank {core.rank(jet_matrix)} and "
            f"kernel dimension {kernel.shape[0]}, expected 9 and 1"
        )
    q = {index: kernel[0, position] for position, index in enumerate(indices)}
    if q[(2, 2)] == 0:
        raise ZeroDivisionError("q22 vanishes in adapted Petri relation")
    invariant = q[(1, 3)] / q[(2, 2)]
    return invariant, {
        "local_chart": expansion.chart,
        "original_local_parameter_variable": expansion.free_variable,
        "adapted_parameter_definition": "t^23=f_local",
        "adapted_parameter_unique": True,
        "uniqueness_reason": f"gcd(23,{field.order}-1)=1",
        "map_unit": [core.serialize_element(value) for value in adapted["unit"]],
        "map_unit_23rd_root": [
            core.serialize_element(value) for value in adapted["unit_root"]
        ],
        "t_as_series_in_original_parameter": [
            core.serialize_element(value) for value in adapted["t_of_u"]
        ],
        "original_parameter_as_series_in_t": [
            core.serialize_element(value) for value in adapted["u_of_t"]
        ],
        "canonical_leading_matrix": core.serialize_matrix(leading),
        "canonical_echelon_transform": core.serialize_matrix(transform),
        "adjunction_jacobian_in_adapted_parameter": [
            core.serialize_element(value) for value in jacobian_t
        ],
        "du_dt": [core.serialize_element(value) for value in du_dt],
        "common_differential_factor_du_dt_over_jacobian": [
            core.serialize_element(value) for value in common_differential_factor
        ],
        "echelon_series_orders_0_through_12": [
            [core.serialize_element(value) for value in row] for row in echelon
        ],
        "quadratic_jet_matrix_shape": list(jet_matrix.shape),
        "quadratic_jet_rank": core.rank(jet_matrix),
        "quadratic_kernel_dimension": int(kernel.shape[0]),
        "quadratic_coefficients": {
            str(index): core.serialize_element(value) for index, value in q.items()
        },
        "q13": core.serialize_element(q[(1, 3)]),
        "q22": core.serialize_element(q[(2, 2)]),
        "J": core.serialize_element(invariant),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model",
        type=Path,
        default=Path(
            "outputs/m23-proof-closeout-20260902/agents/"
            "petri-model-exactify/marked-point-cubic-completion.json"
        ),
    )
    parser.add_argument("--map", type=Path, default=HERE / "exact-f13-map.json")
    parser.add_argument("--output", type=Path, default=HERE / "petri-binding-mod13.json")
    args = parser.parse_args()

    model_payload = json.loads(args.model.read_text())
    map_payload = json.loads(args.map.read_text())
    if map_payload["status"] != "exact_f13_degree23_candidate_constructed":
        raise ValueError("exact map receipt is not accepted")

    field, alpha = core.make_field()
    quadric = core.reduce_polynomial(
        model_payload["quadric_coefficients"], field, alpha
    )
    cubic = core.reduce_polynomial(
        model_payload["cubic_coefficients"], field, alpha
    )
    p0 = field([1, 0, 0, 0])
    pinfinity = field(
        [
            core.reduce_field_element(coordinate, field, alpha)
            for coordinate in model_payload["marked_point"]["coordinates"]
        ]
    )
    order = 64
    expansion0 = core.local_expansion(quadric, cubic, p0, 0, order, field)
    expansion_infinity = core.local_expansion(
        quadric, cubic, pinfinity, 3, order, field
    )
    numerator = deserialize_form(map_payload["map"]["numerator"], 5, field, alpha)
    denominator = deserialize_form(map_payload["map"]["denominator"], 5, field, alpha)
    monomials0 = [
        core.series_monomial(expansion0.coordinates, monomial)
        for monomial in core.MONOMIALS[5]
    ]
    monomials_infinity = [
        core.series_monomial(expansion_infinity.coordinates, monomial)
        for monomial in core.MONOMIALS[5]
    ]
    n0 = core.evaluate_form_series(numerator, monomials0, field)
    d0 = core.evaluate_form_series(denominator, monomials0, field)
    ni = core.evaluate_form_series(numerator, monomials_infinity, field)
    di = core.evaluate_form_series(denominator, monomials_infinity, field)
    f0 = series_divide(n0, d0)
    reciprocal_infinity = series_divide(di, ni)

    j0, receipt0 = marked_invariant(expansion0, f0, quadric, cubic, field)
    jinfinity, receipt_infinity = marked_invariant(
        expansion_infinity, reciprocal_infinity, quadric, cubic, field
    )
    u1 = j0 + jinfinity + j0 * jinfinity
    theta = field(3) + field(4) * alpha**2 + field(2) * alpha**3
    theta += field(9) * alpha**4 + field(11) * alpha**5
    matches = bool(u1 == theta)
    orbit_size = next(
        degree for degree in range(1, 7) if u1 ** (core.PRIME**degree) == u1
    )
    minimal = u1.minimal_poly()

    result = {
        "status": (
            "exact_petri_binding_succeeded"
            if matches and orbit_size == 6
            else "exact_petri_binding_failed"
        ),
        "inputs": {
            "model": str(args.model.resolve()),
            "model_sha256": sha256(args.model),
            "map": str(args.map.resolve()),
            "map_sha256": sha256(args.map),
        },
        "field": {
            "characteristic": core.PRIME,
            "degree": 6,
            "order": field.order,
            "gcd_23_q_minus_1": int(np.gcd(23, field.order - 1)),
        },
        "P0": receipt0,
        "Pinf": receipt_infinity,
        "descent": {
            "formula": "U1=J0+Jinf+J0*Jinf",
            "J0": core.serialize_element(j0),
            "Jinf": core.serialize_element(jinfinity),
            "U1": core.serialize_element(u1),
            "expected_theta": core.serialize_element(theta),
            "equals_expected_theta": matches,
            "frobenius_orbit_size": orbit_size,
            "minimal_polynomial_coefficients_descending_mod_13": [
                int(coefficient) for coefficient in minimal.coeffs
            ],
        },
        "boundary": {
            "proves": (
                "the exact reduced cover's adapted marked Petri coordinate "
                "equals the certified primitive degree-six residue generator"
            ),
            "does_not_prove": (
                "by itself, characteristic-zero descent, M23 geometric "
                "monodromy, or existence of a chosen integral Hurwitz model"
            ),
        },
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(args.output)
    print("status", result["status"])
    print("J0", result["descent"]["J0"])
    print("Jinf", result["descent"]["Jinf"])
    print("U1", result["descent"]["U1"])
    print("theta", result["descent"]["expected_theta"])
    print("orbit", orbit_size, "minimal", result["descent"]["minimal_polynomial_coefficients_descending_mod_13"])


if __name__ == "__main__":
    main()
