#!/usr/bin/env python3
"""Emit an exact Magma probe after translating an order-23 place to zero.

Magma 2.29-9 hits an internal assertion for nonzero linear ``Prime=t-c`` on
this input.  The exact specialization scan shows that t=5 is unramified and
has factorization [23].  Replacing t by s+5 makes the same place ``Prime=s``;
this both avoids the translated-linear-prime code path and gives the subgroup
sieve a 23-cycle immediately.  The underlying function-field extension is
unchanged.  ``ShortOK`` remains false.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from magma_function_field_probe import build_script


HERE = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=Path, default=HERE / "function-field-polynomial.json"
    )
    parser.add_argument(
        "--output", type=Path, default=HERE / "magma-function-field-shift5.m"
    )
    parser.add_argument("--proof-effort", type=int, default=1)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    script = build_script(payload, args.proof_effort)
    script = script.replace(
        'SetVerbose("GaloisGroup", 1);',
        'SetVerbose("GaloisGroup", 5);\nSetVerbose("Invariant", 3);',
        1,
    )
    script = script.replace(
        "K<t> := FunctionField(k);",
        "K<s> := FunctionField(k);\nt := s + 5;",
        1,
    )
    old_call = (
        "GaloisGroup(f : ShortOK := false, "
        f"ProofEffort := {args.proof_effort})"
    )
    new_call = (
        "GaloisGroup(f : Prime := s, ShortOK := false, "
        f"ProofEffort := {args.proof_effort})"
    )
    if old_call not in script:
        raise ValueError("generator layout changed: GaloisGroup call absent")
    script = script.replace(old_call, new_call, 1)
    input_marker = 'printf "INPUT_BINDING'
    script = script.replace(
        input_marker,
        'printf "BASE_CHANGE original_t=s+5 initial_prime=s expected_pattern=[23]\\n";\n'
        + input_marker,
        1,
    )
    args.output.write_text(script)

    receipt = {
        "status": "magma_shifted_prime_probe_emitted",
        "input": str(args.input.resolve()),
        "input_sha256": sha256(args.input),
        "output": str(args.output.resolve()),
        "output_sha256": sha256(args.output),
        "bytes": len(script.encode()),
        "base_change": "original_t=s+5",
        "prime": "s",
        "specialization_factor_degrees": [23],
        "specialization_source": str((HERE / "galois-prime-scan.json").resolve()),
        "proof_effort": args.proof_effort,
        "short_ok": False,
        "seed": 1,
        "scope": (
            "Only a completed CERTIFIED line proves M23; timeout/error leaves "
            "the exact interval M23 <= G_geom <= A23 unchanged."
        ),
    }
    receipt_path = HERE / "magma-function-field-shift5-emission.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
