#!/usr/bin/env python3
"""Construct the characteristic-13 degree-23 map from exact Petri data.

This verifier uses only the exact marked-point quadric/cubic certificate.  It
reduces at the prime ``(13, s - 4)`` into

    k = F_13[a] / (a^6-a^5-3a^4-3a^3+a^2+5a+4),

constructs exact formal parameters at the two marked points, and computes

    V_0   = H^0(5K - 23 P_0),
    V_inf = H^0(5K - 23 P_inf)

as kernels of the first 23 Taylor coefficients in the 27-dimensional
degree-five canonical quotient.  It then solves the degree-ten identities

    A_i B_j - A_j B_i = 0

for a multiplier matrix from V_0 to V_inf.  A one-dimensional multiplier
kernel and the exact local orders produce a candidate function f=N/D.

This script deliberately stops before claiming a Belyi passport.  The output
records only algebraic equalities actually checked here; common-base-divisor,
third-fibre, monodromy, and moduli-lifting checks are separate gates.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

import galois
import numpy as np


PRIME = 13
S_RESIDUE = 4
COMPACT_DESCENDING = (1, -1, -3, -3, 1, 5, 4)


def monomials(degree: int) -> tuple[tuple[int, ...], ...]:
    return tuple(itertools.combinations_with_replacement(range(4), degree))


MONOMIALS = {degree: monomials(degree) for degree in range(11)}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_monomial(text: str) -> tuple[int, ...]:
    return tuple(int(value.strip()) for value in text.strip("()").split(","))


def make_field() -> tuple[type, Any]:
    prime_field = galois.GF(PRIME)
    compact = galois.Poly(
        [coefficient % PRIME for coefficient in COMPACT_DESCENDING],
        field=prime_field,
    )
    if not compact.is_irreducible():
        raise AssertionError("compact sextic is not irreducible modulo 13")
    field = galois.GF(PRIME**6, irreducible_poly=compact)
    # In galois' integer representation, PRIME represents the residue class x.
    alpha = field(PRIME)
    value = field(0)
    for coefficient in COMPACT_DESCENDING:
        value = value * alpha + field(coefficient % PRIME)
    if value != 0:
        raise AssertionError("selected alpha does not satisfy the compact sextic")
    return field, alpha


def rational_mod(text: str) -> int:
    value = Fraction(text)
    denominator = value.denominator % PRIME
    if denominator == 0:
        raise ZeroDivisionError(f"nonintegral coefficient at 13: {text}")
    return value.numerator % PRIME * pow(denominator, -1, PRIME) % PRIME


def reduce_field_element(payload: dict[str, Any], field: type, alpha: Any) -> Any:
    result = field(0)
    power = field(1)
    for a_text, b_text in zip(
        payload["A_coefficients_ascending"],
        payload["B_coefficients_ascending"],
        strict=True,
    ):
        residue = (rational_mod(a_text) + S_RESIDUE * rational_mod(b_text)) % PRIME
        result += field(residue) * power
        power *= alpha
    return result


def reduce_polynomial(
    payload: dict[str, Any], field: type, alpha: Any
) -> dict[tuple[int, ...], Any]:
    return {
        parse_monomial(monomial): reduce_field_element(coefficient, field, alpha)
        for monomial, coefficient in payload.items()
    }


def evaluate_polynomial(
    polynomial: dict[tuple[int, ...], Any], point: Any, field: type
) -> Any:
    result = field(0)
    for monomial, coefficient in polynomial.items():
        # ``galois`` extension-field scalars are mutable 0-D arrays.  Never
        # apply an in-place operation to an alias of a stored coefficient.
        term = field(coefficient)
        for variable in monomial:
            term = term * point[variable]
        result = result + term
    return result


def derivative_value(
    polynomial: dict[tuple[int, ...], Any], variable: int, point: Any, field: type
) -> Any:
    result = field(0)
    for monomial, coefficient in polynomial.items():
        multiplicity = monomial.count(variable)
        if multiplicity == 0:
            continue
        removed = list(monomial)
        removed.remove(variable)
        term = field(multiplicity) * coefficient
        for index in removed:
            term *= point[index]
        result += term
    return result


def series_mul(left: Any, right: Any) -> Any:
    # galois implements numpy convolution in compiled finite-field kernels.
    # This is exactly the same truncated Cauchy product as the scalar loop,
    # but is orders of magnitude faster for extension-field series.
    return np.convolve(left, right)[: left.size]


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


@dataclass
class LocalExpansion:
    chart: int
    free_variable: int
    solved_variables: tuple[int, int]
    normalized_point: Any
    coordinates: list[Any]
    jacobian_determinant: Any


def local_expansion(
    quadric: dict[tuple[int, ...], Any],
    cubic: dict[tuple[int, ...], Any],
    point: Any,
    chart: int,
    order: int,
    field: type,
) -> LocalExpansion:
    if point[chart] == 0:
        raise ValueError(f"point is not in chart x{chart} != 0")
    normalized = point / point[chart]
    affine_variables = tuple(index for index in range(4) if index != chart)

    choice = None
    for free_variable in affine_variables:
        solved = tuple(index for index in affine_variables if index != free_variable)
        q0 = derivative_value(quadric, solved[0], normalized, field)
        q1 = derivative_value(quadric, solved[1], normalized, field)
        c0 = derivative_value(cubic, solved[0], normalized, field)
        c1 = derivative_value(cubic, solved[1], normalized, field)
        determinant = q0 * c1 - q1 * c0
        if determinant != 0:
            choice = (free_variable, solved, (q0, q1, c0, c1), determinant)
            break
    if choice is None:
        raise ValueError("no implicit-function Jacobian minor is invertible")

    free_variable, solved, jacobian, determinant = choice
    q0, q1, c0, c1 = jacobian
    coordinates = [field.Zeros(order) for _ in range(4)]
    for variable in range(4):
        coordinates[variable][0] = normalized[variable]
    coordinates[free_variable][1] = 1

    for degree in range(1, order):
        q_residual = polynomial_series(quadric, coordinates)[degree]
        c_residual = polynomial_series(cubic, coordinates)[degree]
        # Inverse of [[q0,q1],[c0,c1]], applied to -residual.
        coordinates[solved[0]][degree] = (
            -c1 * q_residual + q1 * c_residual
        ) / determinant
        coordinates[solved[1]][degree] = (
            c0 * q_residual - q0 * c_residual
        ) / determinant

    q_check = polynomial_series(quadric, coordinates)
    c_check = polynomial_series(cubic, coordinates)
    if np.any(q_check != 0) or np.any(c_check != 0):
        raise AssertionError("formal implicit-function expansion failed")
    return LocalExpansion(
        chart=chart,
        free_variable=free_variable,
        solved_variables=solved,
        normalized_point=normalized,
        coordinates=coordinates,
        jacobian_determinant=determinant,
    )


def rank(matrix: Any) -> int:
    return int(matrix.row_space().shape[0])


def polynomial_vector(
    polynomial: dict[tuple[int, ...], Any], degree: int, field: type
) -> Any:
    index = {monomial: position for position, monomial in enumerate(MONOMIALS[degree])}
    vector = field.Zeros(len(index))
    for monomial, coefficient in polynomial.items():
        vector[index[monomial]] += coefficient
    return vector


def multiply_by_monomial(
    polynomial: dict[tuple[int, ...], Any],
    factor: tuple[int, ...],
    degree: int,
    field: type,
) -> Any:
    index = {monomial: position for position, monomial in enumerate(MONOMIALS[degree])}
    result = field.Zeros(len(index))
    for monomial, coefficient in polynomial.items():
        result[index[tuple(sorted((*monomial, *factor)))]] += coefficient
    return result


@dataclass
class QuotientReducer:
    degree: int
    ideal_rank: int
    rows: Any
    pivots: tuple[int, ...]
    free_columns: tuple[int, ...]

    def reduce(self, vector: Any) -> Any:
        result = vector.copy()
        for row, pivot in zip(self.rows, self.pivots, strict=True):
            result -= result[pivot] * row
        if np.any(result[list(self.pivots)] != 0):
            raise AssertionError("RREF quotient reduction retained pivot entries")
        return result[list(self.free_columns)]


def canonical_quotient(
    quadric: dict[tuple[int, ...], Any],
    cubic: dict[tuple[int, ...], Any],
    degree: int,
    field: type,
) -> QuotientReducer:
    if degree < 3:
        raise ValueError("canonical ideal computation requires degree at least 3")
    rows = [
        multiply_by_monomial(quadric, factor, degree, field)
        for factor in MONOMIALS[degree - 2]
    ]
    rows.extend(
        multiply_by_monomial(cubic, factor, degree, field)
        for factor in MONOMIALS[degree - 3]
    )
    matrix = field(rows)
    rref = matrix.row_reduce()
    nonzero_rows = []
    pivots = []
    for row in rref:
        nonzero = np.flatnonzero(row != 0)
        if nonzero.size:
            nonzero_rows.append(row)
            pivots.append(int(nonzero[0]))
    pivot_set = set(pivots)
    free_columns = tuple(
        index for index in range(len(MONOMIALS[degree])) if index not in pivot_set
    )
    return QuotientReducer(
        degree=degree,
        ideal_rank=len(nonzero_rows),
        rows=field(nonzero_rows),
        pivots=tuple(pivots),
        free_columns=free_columns,
    )


def degree_five_space(
    expansion: LocalExpansion,
    reducer: QuotientReducer,
    field: type,
) -> tuple[Any, dict[str, Any], list[Any]]:
    if reducer.degree != 5:
        raise ValueError("wrong quotient reducer degree")
    all_series = [
        series_monomial(expansion.coordinates, monomial) for monomial in MONOMIALS[5]
    ]
    quotient_series = [all_series[index] for index in reducer.free_columns]
    jets = field(
        [
            [series[degree] for series in quotient_series]
            for degree in range(23)
        ]
    )
    kernel = jets.null_space()
    full_basis = field.Zeros((kernel.shape[0], len(MONOMIALS[5])))
    for row in range(kernel.shape[0]):
        full_basis[row, list(reducer.free_columns)] = kernel[row]
    checks = {
        "jet_matrix_shape": list(jets.shape),
        "jet_rank": rank(jets),
        "kernel_dimension": int(kernel.shape[0]),
        "kernel_rank": rank(kernel),
    }
    return full_basis, checks, all_series


def osculation_kernel(
    expansion: LocalExpansion,
    reducer: QuotientReducer,
    jet_count: int,
    field: type,
) -> tuple[Any, dict[str, Any]]:
    all_series = [
        series_monomial(expansion.coordinates, monomial)
        for monomial in MONOMIALS[reducer.degree]
    ]
    quotient_series = [all_series[index] for index in reducer.free_columns]
    jets = field(
        [
            [series[degree] for series in quotient_series]
            for degree in range(jet_count)
        ]
    )
    kernel = jets.null_space()
    return kernel, {
        "degree": reducer.degree,
        "jet_matrix_shape": list(jets.shape),
        "jet_rank": rank(jets),
        "kernel_dimension": int(kernel.shape[0]),
        "kernel_rank": rank(kernel),
    }


def multiply_forms(left: Any, right: Any, field: type) -> Any:
    index = {monomial: position for position, monomial in enumerate(MONOMIALS[10])}
    result = field.Zeros(len(MONOMIALS[10]))
    left_nonzero = np.flatnonzero(left != 0)
    right_nonzero = np.flatnonzero(right != 0)
    for left_index in left_nonzero:
        for right_index in right_nonzero:
            monomial = tuple(
                sorted((*MONOMIALS[5][int(left_index)], *MONOMIALS[5][int(right_index)]))
            )
            result[index[monomial]] += left[left_index] * right[right_index]
    return result


def recover_multiplier(
    zero_space: Any,
    infinity_space: Any,
    reducer10: QuotientReducer,
    field: type,
) -> tuple[Any, Any, dict[str, Any]]:
    products: dict[tuple[int, int], Any] = {}
    for plus in range(4):
        for minus in range(4):
            products[(plus, minus)] = reducer10.reduce(
                multiply_forms(zero_space[plus], infinity_space[minus], field)
            )

    equations = []
    quotient_dimension = len(reducer10.free_columns)
    for left, right in itertools.combinations(range(4), 2):
        for coordinate in range(quotient_dimension):
            row = field.Zeros(16)
            for plus in range(4):
                row[4 * left + plus] += products[(plus, right)][coordinate]
                row[4 * right + plus] -= products[(plus, left)][coordinate]
            equations.append(row)
    system = field(equations)
    kernel = system.null_space()
    checks = {
        "linear_system_shape": list(system.shape),
        "linear_system_rank": rank(system),
        "kernel_dimension": int(kernel.shape[0]),
        "kernel_rank": rank(kernel),
    }
    if kernel.shape[0] != 1:
        return kernel, system, checks
    vector = kernel[0]
    pivot = int(np.flatnonzero(vector != 0)[0])
    vector /= vector[pivot]
    multiplier = vector.reshape((4, 4))

    maximum_residual_nonzero = False
    for left, right in itertools.combinations(range(4), 2):
        residual = field.Zeros(quotient_dimension)
        for plus in range(4):
            residual += multiplier[left, plus] * products[(plus, right)]
            residual -= multiplier[right, plus] * products[(plus, left)]
        maximum_residual_nonzero |= bool(np.any(residual != 0))
    checks["all_cross_products_zero_in_H0_10K"] = not maximum_residual_nonzero
    checks["normalization_flat_index"] = pivot
    return multiplier, system, checks


def recover_multiplier_by_jets(
    zero_space: Any,
    infinity_space: Any,
    monomial_series_at_zero: list[Any],
    field: type,
) -> tuple[Any, Any, dict[str, Any]]:
    """Recover the multiplier and prove its cross-products globally.

    Every cross-product is a section of 10K, whose degree is 60 on a
    genus-four curve.  Thus a cross-product vanishing through Taylor degree
    60 at one smooth point is identically zero.  This replaces a large
    degree-ten Macaulay reduction by an equivalent exact finite-jet test.
    """

    zero_series = [
        evaluate_form_series(form, monomial_series_at_zero, field)
        for form in zero_space
    ]
    infinity_series = [
        evaluate_form_series(form, monomial_series_at_zero, field)
        for form in infinity_space
    ]
    products = {
        (plus, minus): series_mul(zero_series[plus], infinity_series[minus])
        for plus in range(4)
        for minus in range(4)
    }
    # Coefficients 0,...,60 prove order at least 61.  The first 23 are
    # automatically zero because every V0 section vanishes at P0, but keeping
    # them makes the certificate self-contained.
    equations = []
    for left, right in itertools.combinations(range(4), 2):
        for degree in range(61):
            row = field.Zeros(16)
            for plus in range(4):
                row[4 * left + plus] += products[(plus, right)][degree]
                row[4 * right + plus] -= products[(plus, left)][degree]
            equations.append(row)
    system = field(equations)
    kernel = system.null_space()
    checks = {
        "method": "exact Taylor identities at P0",
        "linear_system_shape": list(system.shape),
        "linear_system_rank": rank(system),
        "kernel_dimension": int(kernel.shape[0]),
        "kernel_rank": rank(kernel),
        "cross_product_bundle": "10K",
        "degree_10K": 60,
        "coefficients_imposed": [0, 60],
        "vanishing_order_lower_bound": 61,
    }
    if kernel.shape[0] != 1:
        return kernel, system, checks
    vector = kernel[0]
    pivot = int(np.flatnonzero(vector != 0)[0])
    vector /= vector[pivot]
    multiplier = vector.reshape((4, 4))

    all_zero = True
    residual_through_available_order = True
    available_order = monomial_series_at_zero[0].size
    for left, right in itertools.combinations(range(4), 2):
        residual = field.Zeros(available_order)
        for plus in range(4):
            residual += multiplier[left, plus] * products[(plus, right)]
            residual -= multiplier[right, plus] * products[(plus, left)]
        all_zero &= bool(np.all(residual[:61] == 0))
        residual_through_available_order &= bool(np.all(residual == 0))
    checks.update(
        {
            "normalization_flat_index": pivot,
            "all_cross_products_zero_through_degree_60": all_zero,
            "all_cross_products_zero_through_available_series_order": (
                residual_through_available_order
            ),
            "available_series_coefficients": available_order,
            "global_cross_product_identity_proved": (
                all_zero and 61 > 60
            ),
            "global_identity_reason": (
                "a nonzero section of 10K has a zero divisor of degree 60, "
                "so vanishing order at least 61 at P0 forces the section to be zero"
            ),
        }
    )
    return multiplier, system, checks


def evaluate_form_series(form: Any, monomial_series: list[Any], field: type) -> Any:
    result = field.Zeros(monomial_series[0].size)
    for index in np.flatnonzero(form != 0):
        result += form[index] * monomial_series[int(index)]
    return result


def first_nonzero(series: Any) -> int | None:
    indices = np.flatnonzero(series != 0)
    return int(indices[0]) if indices.size else None


def serialize_element(value: Any) -> list[int]:
    return [int(coefficient) for coefficient in value.vector()]


def serialize_point(point: Any) -> list[list[int]]:
    return [serialize_element(value) for value in point]


def serialize_form(form: Any, degree: int) -> dict[str, list[int]]:
    return {
        str(MONOMIALS[degree][int(index)]): serialize_element(form[index])
        for index in np.flatnonzero(form != 0)
    }


def serialize_matrix(matrix: Any) -> list[list[list[int]]]:
    return [
        [serialize_element(matrix[row, column]) for column in range(matrix.shape[1])]
        for row in range(matrix.shape[0])
    ]


def local_receipt(expansion: LocalExpansion, order: int) -> dict[str, Any]:
    return {
        "chart": f"x{expansion.chart}=1",
        "local_parameter": (
            f"t=x{expansion.free_variable}/x{expansion.chart}"
            f"-P[x{expansion.free_variable}/x{expansion.chart}]"
        ),
        "free_variable": expansion.free_variable,
        "solved_variables": list(expansion.solved_variables),
        "normalized_point": serialize_point(expansion.normalized_point),
        "jacobian_determinant": serialize_element(expansion.jacobian_determinant),
        "equations_vanish_through_order": order - 1,
        "coordinate_series": {
            f"x{variable}/x{expansion.chart}": [
                serialize_element(value) for value in expansion.coordinates[variable]
            ]
            for variable in range(4)
        },
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
    parser.add_argument("--series-order", type=int, default=64)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("exact-f13-map.json"),
    )
    args = parser.parse_args()
    if args.series_order < 62:
        raise ValueError(
            "series order must be at least 62 to certify degree-ten identities"
        )

    payload = json.loads(args.model.read_text())
    field, alpha = make_field()
    quadric = reduce_polynomial(payload["quadric_coefficients"], field, alpha)
    cubic = reduce_polynomial(payload["cubic_coefficients"], field, alpha)
    p0 = field([1, 0, 0, 0])
    pinfinity = field(
        [
            reduce_field_element(coordinate, field, alpha)
            for coordinate in payload["marked_point"]["coordinates"]
        ]
    )
    input_checks = {
        "Q_at_P0_zero": bool(evaluate_polynomial(quadric, p0, field) == 0),
        "C_at_P0_zero": bool(evaluate_polynomial(cubic, p0, field) == 0),
        "Q_at_Pinf_zero": bool(evaluate_polynomial(quadric, pinfinity, field) == 0),
        "C_at_Pinf_zero": bool(evaluate_polynomial(cubic, pinfinity, field) == 0),
        "points_distinct": bool(np.any(p0 * pinfinity[3] != pinfinity * p0[3])),
    }
    if not all(input_checks.values()):
        raise ValueError(f"reduced marked-point input failed: {input_checks}")

    expansion0 = local_expansion(
        quadric, cubic, p0, chart=0, order=args.series_order, field=field
    )
    expansion_infinity = local_expansion(
        quadric, cubic, pinfinity, chart=3, order=args.series_order, field=field
    )

    reducer4 = canonical_quotient(quadric, cubic, 4, field)
    if reducer4.ideal_rank != 14 or len(reducer4.free_columns) != 21:
        raise ValueError(
            f"degree-four quotient has ideal rank {reducer4.ideal_rank} "
            f"and dimension {len(reducer4.free_columns)}, expected 14 and 21"
        )
    four_k_zero, four_k_zero_checks = osculation_kernel(
        expansion0, reducer4, 23, field
    )
    four_k_infinity, four_k_infinity_checks = osculation_kernel(
        expansion_infinity, reducer4, 23, field
    )
    if four_k_zero.shape[0] != 0 or four_k_infinity.shape[0] != 0:
        raise ValueError(
            "unexpected effective divisor in |4K-23P|: "
            f"dimensions {four_k_zero.shape[0]}, {four_k_infinity.shape[0]}"
        )

    reducer5 = canonical_quotient(quadric, cubic, 5, field)
    if reducer5.ideal_rank != 29 or len(reducer5.free_columns) != 27:
        raise ValueError(
            f"degree-five quotient has ideal rank {reducer5.ideal_rank} "
            f"and dimension {len(reducer5.free_columns)}, expected 29 and 27"
        )
    zero_space, zero_checks, monomial_series0 = degree_five_space(
        expansion0, reducer5, field
    )
    infinity_space, infinity_checks, monomial_series_infinity = degree_five_space(
        expansion_infinity, reducer5, field
    )
    if zero_space.shape[0] != 4 or infinity_space.shape[0] != 4:
        raise ValueError(
            "23-jet spaces did not both have dimension four: "
            f"{zero_space.shape[0]}, {infinity_space.shape[0]}"
        )

    multiplier, _, multiplier_checks = recover_multiplier_by_jets(
        zero_space, infinity_space, monomial_series0, field
    )

    status = "exact_f13_multiplier_gate_failed"
    map_payload = None
    if multiplier_checks["kernel_dimension"] == 1:
        multiplier_determinant = np.linalg.det(multiplier)
        multiplier_checks["determinant"] = serialize_element(multiplier_determinant)
        multiplier_checks["invertible"] = bool(multiplier_determinant != 0)
        if multiplier_determinant == 0:
            raise ValueError("the unique multiplier matrix is singular")
        row_diagnostics = []
        chosen = None
        for row in range(4):
            numerator = field.Zeros(len(MONOMIALS[5]))
            for plus in range(4):
                numerator += multiplier[row, plus] * zero_space[plus]
            denominator = infinity_space[row].copy()
            n0 = evaluate_form_series(numerator, monomial_series0, field)
            d0 = evaluate_form_series(denominator, monomial_series0, field)
            ni = evaluate_form_series(numerator, monomial_series_infinity, field)
            di = evaluate_form_series(denominator, monomial_series_infinity, field)
            diagnostic = {
                "row": row,
                "ord_P0_N": first_nonzero(n0),
                "ord_P0_D": first_nonzero(d0),
                "ord_Pinf_N": first_nonzero(ni),
                "ord_Pinf_D": first_nonzero(di),
            }
            diagnostic["accepted"] = (
                diagnostic["ord_P0_N"] == 23
                and diagnostic["ord_P0_D"] == 0
                and diagnostic["ord_Pinf_N"] == 0
                and diagnostic["ord_Pinf_D"] == 23
            )
            row_diagnostics.append(diagnostic)
            if diagnostic["accepted"] and chosen is None:
                chosen = (row, numerator, denominator, n0, d0, ni, di)

        if chosen is None:
            status = "exact_f13_multiplier_found_but_no_nondegenerate_row"
            map_payload = {"row_diagnostics": row_diagnostics}
        else:
            row, numerator, denominator, n0, d0, ni, di = chosen
            pivot = int(np.flatnonzero(denominator != 0)[0])
            scale = denominator[pivot]
            numerator /= scale
            denominator /= scale
            n0 /= scale
            d0 /= scale
            ni /= scale
            di /= scale
            status = "exact_f13_degree23_candidate_constructed"
            map_payload = {
                "interpretation": "f=N/D, with zero at P0 and pole at Pinf",
                "function_divisor": "23*P0 - 23*Pinf",
                "function_divisor_proved": True,
                "degree": 23,
                "separable": True,
                "separability_reason": (
                    "an inseparable degree in characteristic 13 is divisible by 13, "
                    "whereas the proved map degree is 23"
                ),
                "divisor_reason": (
                    "the invertible multiplier identifies the full spaces "
                    "H0(5K-23Pinf) and H0(5K-23P0); their residual degree-seven "
                    "linear systems are basepoint-free by the two vanishing "
                    "H0(4K-23P) screens, so comparing fixed divisors gives "
                    "div(f)=23P0-23Pinf"
                ),
                "selected_multiplier_row": row,
                "row_diagnostics": row_diagnostics,
                "normalization_denominator_monomial": str(MONOMIALS[5][pivot]),
                "numerator": serialize_form(numerator, 5),
                "denominator": serialize_form(denominator, 5),
                "local_orders": {
                    "ord_P0_N": first_nonzero(n0),
                    "ord_P0_D": first_nonzero(d0),
                    "ord_Pinf_N": first_nonzero(ni),
                    "ord_Pinf_D": first_nonzero(di),
                },
                "local_coefficients_at_orders_0_through_24": {
                    "N_at_P0": [serialize_element(value) for value in n0[:25]],
                    "D_at_P0": [serialize_element(value) for value in d0[:25]],
                    "N_at_Pinf": [serialize_element(value) for value in ni[:25]],
                    "D_at_Pinf": [serialize_element(value) for value in di[:25]],
                },
            }

    result = {
        "status": status,
        "input": {
            "path": str(args.model.resolve()),
            "sha256": sha256(args.model),
            "source_status": payload["status"],
        },
        "field": {
            "characteristic": PRIME,
            "degree": 6,
            "order": PRIME**6,
            "compact_polynomial_coefficients_descending_mod_13": [
                value % PRIME for value in COMPACT_DESCENDING
            ],
            "compact_polynomial_irreducible": True,
            "alpha_vector_high_to_low": serialize_element(alpha),
            "s_residue": S_RESIDUE,
            "s_relation_check": (S_RESIDUE**2 + 23) % PRIME == 0,
        },
        "reduced_input_checks": input_checks,
        "reduced_model": {
            "quadric": {
                str(monomial): serialize_element(value)
                for monomial, value in quadric.items()
            },
            "cubic": {
                str(monomial): serialize_element(value)
                for monomial, value in cubic.items()
            },
            "P0": serialize_point(p0),
            "Pinf": serialize_point(pinfinity),
        },
        "local_expansions": {
            "P0": local_receipt(expansion0, args.series_order),
            "Pinf": local_receipt(expansion_infinity, args.series_order),
        },
        "degree_four_basepoint_screen": {
            "ambient_monomials": len(MONOMIALS[4]),
            "canonical_ideal_rank": reducer4.ideal_rank,
            "quotient_dimension": len(reducer4.free_columns),
            "P0": four_k_zero_checks,
            "Pinf": four_k_infinity_checks,
            "both_H0_4K_minus_23P_zero": True,
            "consequence": (
                "the residual degree-seven systems |5K-23P0| and "
                "|5K-23Pinf| are basepoint-free"
            ),
            "reason": (
                "a base point R of |5K-23P| would force "
                "5K-23P = K+R and hence a section of 4K-23P"
            ),
        },
        "degree_five_quotient": {
            "ambient_monomials": len(MONOMIALS[5]),
            "canonical_ideal_rank": reducer5.ideal_rank,
            "quotient_dimension": len(reducer5.free_columns),
            "pivot_columns": list(reducer5.pivots),
            "free_columns": list(reducer5.free_columns),
            "P0_space": zero_checks,
            "Pinf_space": infinity_checks,
            "P0_basis": [serialize_form(row, 5) for row in zero_space],
            "Pinf_basis": [serialize_form(row, 5) for row in infinity_space],
        },
        "degree_ten_quotient": {
            "ambient_monomials": len(MONOMIALS[10]),
            "quotient_dimension": 57,
            "dimension_reason": "h0(10K)=deg(10K)-g+1=60-4+1=57",
            "identity_test": "61 exact Taylor coefficients at P0",
        },
        "multiplier": {
            **multiplier_checks,
            "matrix": (
                serialize_matrix(multiplier)
                if multiplier_checks["kernel_dimension"] == 1
                else None
            ),
        },
        "map": map_payload,
        "boundary": {
            "proves": (
                "the reported exact finite-field jet spaces, quotient-ring "
                "syzygies, and local orders, if status says the candidate was constructed"
            ),
            "does_not_prove": (
                "the three-point passport; the residual branch fibre, monodromy, "
                "the M22 quotient interpretation, and tame lifting remain unverified"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(args.output)
    print("status", status)
    print("degree-5", reducer5.ideal_rank, len(reducer5.free_columns), zero_checks, infinity_checks)
    print("degree-4", four_k_zero_checks, four_k_infinity_checks)
    print("degree-10 identity test", multiplier_checks)
    print("multiplier", multiplier_checks)
    if map_payload is not None:
        print("local", map_payload.get("local_orders"), map_payload.get("row_diagnostics"))


if __name__ == "__main__":
    main()
