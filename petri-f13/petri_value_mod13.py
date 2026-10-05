#!/usr/bin/env python3
"""Verify the descended Petri generator has degree six modulo 13.

This is independent of the still-missing cover binding.  It reduces the
certified compact root map theta_U=h(alpha), constructs
F_13[alpha]/g(alpha), and checks the exact Frobenius orbit and minimal
polynomial of h(alpha).
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import galois


P = 13
G_ASCENDING = (4, 5, 1, -3, -3, -1, 1)
ROOT_MAP_DENOMINATOR = 50019668646075143
ROOT_MAP_NUMERATORS = (
    -61970328454788839,
    14793393724122520,
    25097895504098152,
    22903915474596992,
    -29458505838561912,
    7568631896021296,
)
PETRI_ASCENDING = (
    Fraction(-139453053861313073, ROOT_MAP_DENOMINATOR),
    Fraction(-868575878997768238, ROOT_MAP_DENOMINATOR),
    Fraction(-759255603013669319, ROOT_MAP_DENOMINATOR),
    Fraction(-55027953306246948, ROOT_MAP_DENOMINATOR),
    Fraction(353749018385880145, ROOT_MAP_DENOMINATOR),
    Fraction(206305139624557234, ROOT_MAP_DENOMINATOR),
    Fraction(1),
)


def rational_mod_p(value: Fraction) -> int:
    denominator = value.denominator % P
    if denominator == 0:
        raise ValueError(f"denominator is not 13-integral: {value}")
    return (value.numerator % P) * pow(denominator, -1, P) % P


def evaluate(coefficients: list[object], value: object) -> object:
    result = value * 0
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def polynomial_coefficients_ascending(poly: object) -> list[int]:
    return [int(value) for value in reversed(poly.coeffs)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "outputs/m23-proof-closeout-20260902/agents/petri-f13/"
            "petri-value-mod13.json"
        ),
    )
    args = parser.parse_args()

    prime_field = galois.GF(P)
    g_mod = [coefficient % P for coefficient in G_ASCENDING]
    g_poly = galois.Poly(list(reversed(g_mod)), field=prime_field)
    if not g_poly.is_irreducible():
        raise ValueError("compact sextic is reducible modulo 13")

    extension = galois.GF(P**6, irreducible_poly=g_poly)
    # Integer P has base-P polynomial digits [0,1], hence represents x.
    alpha = extension(P)
    g_at_alpha = evaluate([extension(value) for value in g_mod], alpha)
    if g_at_alpha != 0:
        raise ValueError(f"chosen alpha is not the residue class of x: {g_at_alpha}")

    denominator_inverse = pow(ROOT_MAP_DENOMINATOR % P, -1, P)
    h_mod = [
        (numerator % P) * denominator_inverse % P
        for numerator in ROOT_MAP_NUMERATORS
    ]
    theta = evaluate([extension(value) for value in h_mod], alpha)
    petri_mod = [rational_mod_p(value) for value in PETRI_ASCENDING]
    petri_at_theta = evaluate([extension(value) for value in petri_mod], theta)
    if petri_at_theta != 0:
        raise ValueError(f"recognized Petri polynomial does not vanish: {petri_at_theta}")

    orbit = []
    current = theta
    for _ in range(6):
        orbit.append(current)
        current = current**P
    if current != theta or len(set(int(value) for value in orbit)) != 6:
        raise ValueError("Petri value does not have Frobenius orbit length six")

    minimal = theta.minimal_poly()
    if minimal.degree != 6 or not minimal.is_irreducible():
        raise ValueError(f"unexpected minimal polynomial {minimal}")
    minimal_ascending = polynomial_coefficients_ascending(minimal)
    if minimal_ascending != petri_mod:
        raise ValueError(
            f"minimal polynomial {minimal_ascending} != Petri reduction {petri_mod}"
        )

    payload = {
        "status": "exact",
        "field": {
            "definition": "F_13[alpha]/g(alpha)",
            "order": P**6,
            "compact_polynomial_coefficients_ascending_mod13": g_mod,
            "compact_polynomial_irreducible": True,
            "g_at_alpha": int(g_at_alpha),
        },
        "petri_value": {
            "definition": "theta_U=h(alpha)",
            "h_coefficients_ascending_mod13": h_mod,
            "theta_integer_encoding": int(theta),
            "theta_polynomial_vector_ascending": [int(value) for value in theta.vector()[::-1]],
            "petri_polynomial_coefficients_ascending_mod13": petri_mod,
            "minimal_polynomial_coefficients_ascending_mod13": minimal_ascending,
            "minimal_polynomial_degree": int(minimal.degree),
            "petri_polynomial_value": int(petri_at_theta),
        },
        "frobenius": {
            "orbit_integer_encodings": [int(value) for value in orbit],
            "orbit_size": len(set(int(value) for value in orbit)),
            "theta_to_13_power_6_equals_theta": bool(current == theta),
            "proper_subfield_tests": {
                str(degree): bool(theta ** (P**degree) == theta)
                for degree in (1, 2, 3)
            },
        },
        "conclusion": (
            "The exact descended Petri value h(alpha) generates F_(13^6) "
            "over F_13 and is a root of the reduced certified Petri sextic."
        ),
        "boundary": (
            "This certifies the arithmetic value once an actual Hurwitz sheet "
            "is shown to have this Petri coordinate; it does not construct or "
            "bind such a sheet."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)
    print(json.dumps(payload["petri_value"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
