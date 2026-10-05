#!/usr/bin/env python3
"""Scan deterministic linear primes for the degree-23 function-field factor.

The output is exact factorisation data over F_(13^6).  It is used only to
choose a low-cost unramified ``Prime := t-c`` for Magma's exact Stauduhar
calculation; it is not itself an upper-bound certificate for monodromy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx


HERE = Path(__file__).resolve().parent


class Field:
    def __init__(self) -> None:
        pc = fmpz_mod_poly_ctx(13)
        modulus = pc([4, 5, 1, -3, -3, -1, 1])
        if not modulus.is_irreducible():
            raise RuntimeError("field modulus is reducible")
        self.ctx = fq_default_ctx(modulus=modulus, var="a", fq_type="FQ_NMOD")
        self.a = self.ctx.gen()

    def zero(self):
        return self.ctx(0)

    def one(self):
        return self.ctx(1)

    def from_high(self, values: list[int]):
        out = self.zero()
        for value in values:
            out = out * self.a + self.ctx(value)
        return out

    def from_index(self, index: int):
        out = self.zero()
        power = self.one()
        for _ in range(6):
            out += self.ctx(index % 13) * power
            index //= 13
            power *= self.a
        return out

    @staticmethod
    def vector(value: Any) -> list[int]:
        result = [int(item) for item in value.to_list()]
        return result + [0] * (6-len(result))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def eval_t_poly(values: list[list[int]], t_value: Any, field: Field):
    out = field.zero()
    for encoded in reversed(values):
        out = out * t_value + field.from_high(encoded)
    return out


def eval_coefficient(value: dict[str, list[list[int]]], t_value: Any, field: Field):
    numerator = eval_t_poly(value["numerator"], t_value, field)
    denominator = eval_t_poly(value["denominator"], t_value, field)
    if denominator == 0:
        raise ZeroDivisionError
    return numerator / denominator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=Path, default=HERE / "function-field-polynomial.json"
    )
    parser.add_argument("--count", type=int, default=4096)
    parser.add_argument(
        "--output", type=Path, default=HERE / "galois-prime-scan.json"
    )
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    factors = [item for item in payload["factors"] if item["degree"] == 23]
    if payload.get("accepted") is not True or len(factors) != 1:
        raise ValueError("expected one accepted degree-23 factor")
    coefficients = factors[0]["coefficients"]

    field = Field()
    polynomial_context = fq_default_poly_ctx(field.ctx)
    scans = []
    histogram: dict[str, int] = {}
    for index in range(args.count):
        t_value = field.from_index(index)
        try:
            specialised = polynomial_context(
                [eval_coefficient(item, t_value, field) for item in coefficients]
            )
        except ZeroDivisionError:
            scans.append({
                "index": index,
                "t_vector_ascending": field.vector(t_value),
                "status": "coefficient_pole",
            })
            continue
        if specialised.degree() != 23 or specialised.derivative().gcd(specialised).degree() > 0:
            scans.append({
                "index": index,
                "t_vector_ascending": field.vector(t_value),
                "status": "bad_or_ramified",
                "degree": specialised.degree(),
            })
            continue
        _content, factorisation = specialised.factor()
        degrees = sorted(
            [factor.degree() for factor, multiplicity in factorisation for _ in range(multiplicity)]
        )
        order = 1
        from math import gcd
        for degree in degrees:
            order = order * degree // gcd(order, degree)
        key = "+".join(str(value) for value in degrees)
        histogram[key] = histogram.get(key, 0) + 1
        scans.append({
            "index": index,
            "t_vector_ascending": field.vector(t_value),
            "status": "unramified",
            "factor_degrees": degrees,
            "frobenius_order": order,
        })

    accepted = [item for item in scans if item["status"] == "unramified"]
    best = sorted(
        accepted,
        key=lambda item: (
            item["frobenius_order"],
            -len(item["factor_degrees"]),
            item["index"],
        ),
    )[:30]
    receipt = {
        "status": "exact_linear_prime_factorisation_scan",
        "scope": "prime selection only; Frobenius lower-bound data cannot exclude A23",
        "input": str(args.input.resolve()),
        "input_sha256": sha256(args.input),
        "requested_count": args.count,
        "unramified_count": len(accepted),
        "bad_or_pole_count": len(scans)-len(accepted),
        "factor_pattern_histogram": dict(sorted(histogram.items())),
        "lowest_order_candidates": best,
        "all_scans": scans,
    }
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: receipt[key] for key in (
        "status", "requested_count", "unramified_count", "bad_or_pole_count",
        "factor_pattern_histogram", "lowest_order_candidates",
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
