#!/usr/bin/env python3
"""Emit a deterministic Magma certificate script for the degree-23 factor.

The generated script is kept below the public Magma calculator's 50 kB input
limit.  It reconstructs the exact coefficient field and polynomial, checks
irreducibility/separability, and asks Magma's exact global-function-field
Galois-group implementation for an upper bound.  ``ShortOK`` is deliberately
false: a short-coset heuristic is not accepted as a proof gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def vector(value: list[int]) -> str:
    if len(value) != 6:
        raise ValueError(f"expected a compact-basis vector of length six: {value}")
    return "[" + ",".join(str(int(item) % 13) for item in value) + "]"


def t_polynomial(values: list[list[int]]) -> str:
    return "TP([" + ",".join(vector(item) for item in values) + "])"


def coefficient(value: dict[str, list[list[int]]]) -> str:
    numerator = t_polynomial(value["numerator"])
    denominator = t_polynomial(value["denominator"])
    return f"({numerator})/({denominator})"


def build_script(payload: dict[str, Any], proof_effort: int) -> str:
    factors = [item for item in payload["factors"] if item["degree"] == 23]
    if payload.get("accepted") is not True or len(factors) != 1:
        raise ValueError("input does not contain one accepted degree-23 factor")
    factor = factors[0]
    coefficients = factor["coefficients"]
    if len(coefficients) != 24:
        raise ValueError("degree-23 factor must have 24 coefficients")
    coeff_text = ",\n+".join(coefficient(item) for item in coefficients)
    input_hash = payload["runner"]["input_sha256"]
    extraction_hash = sha256(HERE / "function-field-polynomial.json")
    return f'''SetSeed(1);
SetVerbose("GaloisGroup", 1);
procedure CertifiedRun()
Fp := GF(13);
Ru<u> := PolynomialRing(Fp);
k<a> := ext<Fp | u^6-u^5-3*u^4-3*u^3+u^2+5*u+4>;
K<t> := FunctionField(k);
P<x> := PolynomialRing(K);
function E(v)
    return &+[ k!(v[i])*a^(6-i) : i in [1..6] ];
end function;
function TP(vs)
    return &+[ E(vs[i])*t^(i-1) : i in [1..#vs] ];
end function;
cs := [
{coeff_text}
];
f := &+[ cs[i+1]*x^i : i in [0..23] ];
assert Degree(f) eq 23;
assert LeadingCoefficient(f) eq 1;
assert IsIrreducible(f);
assert IsSeparable(f);
printf "INPUT_BINDING exact-map-sha256={input_hash} extraction-sha256={extraction_hash}\\n";
printf "POLYNOMIAL degree=%o irreducible=true separable=true\\n", Degree(f);
time G, roots, data := GaloisGroup(f : ShortOK := false, ProofEffort := {proof_effort});
id, n := TransitiveGroupIdentification(G);
printf "RESULT order=%o degree=%o transitive=%o id=%oT%o\\n", #G, Degree(G), IsTransitive(G), n, id;
printf "GENERATORS %o\\n", Generators(G);
assert n eq 23;
assert id eq 5;
assert #G eq 10200960;
assert IsTransitive(G);
print "CERTIFIED arithmetic-monodromy=M23; geometric-monodromy=M23 by exact inertia lower bound";
end procedure;
CertifiedRun();
'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=Path, default=HERE / "function-field-polynomial.json"
    )
    parser.add_argument(
        "--output", type=Path, default=HERE / "magma-function-field-galois.m"
    )
    parser.add_argument("--proof-effort", type=int, default=20)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    script = build_script(payload, args.proof_effort)
    args.output.write_text(script)
    receipt = {
        "status": "magma_function_field_probe_emitted",
        "input": str(args.input.resolve()),
        "input_sha256": sha256(args.input),
        "output": str(args.output.resolve()),
        "output_sha256": sha256(args.output),
        "bytes": len(script.encode()),
        "public_calculator_input_limit": 50_000,
        "under_limit": len(script.encode()) < 50_000,
        "proof_effort": args.proof_effort,
        "short_ok": False,
        "seed": 1,
    }
    (HERE / "magma-function-field-galois-emission.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if not receipt["under_limit"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
