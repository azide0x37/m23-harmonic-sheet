#!/usr/bin/env python3
"""Emit a bounded multi-prime Magma Galois-group probe.

The exact quadratic places below were independently checked irreducible and
unramified in ``magma-quadratic-prime-scan.xml``.  Their specializations cover
every Frobenius order seen in the M23 degree-23 action.  Feeding them through
Magma's ``NextPrime`` hook can collapse the transitive-subgroup sieve before
the expensive A23-to-M23 Stauduhar descent.  ``ShortOK`` remains false.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from magma_function_field_probe import build_script


HERE = Path(__file__).resolve().parent

# (b, c, exact factorization order) for q=T^2+(a+b)T+(a+c).
PRIMES = [
    # Start with a 23-cycle: this certifies transitivity and primitivity to
    # the subgroup sieve before any two-set resolvent is attempted.
    (4, 4, 23),
    (4, 11, 3),
    (0, 6, 4),
    (1, 5, 5),
    (1, 3, 6),
    (0, 1, 7),
    (1, 2, 8),
    (0, 0, 11),
    (0, 4, 14),
    (0, 3, 15),
]


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
        "--output", type=Path, default=HERE / "magma-function-field-multiseed.m"
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
    prime_text = ",\n    ".join(
        f"t^2 + (a+{b})*t + (a+{c})" for b, c, _ in PRIMES
    )
    hook = f'''qs := [
    {prime_text}
];
function DeterministicNextPrime(old)
    position := Index(qs, old);
    if position eq 0 or position eq #qs then
        return qs[1];
    end if;
    return qs[position+1];
end function;
printf "PRIME_SUITE orders=%o count=%o\\n", {[order for _, _, order in PRIMES]}, #qs;
'''
    marker = 'printf "INPUT_BINDING'
    if marker not in script:
        raise ValueError("generator layout changed: input marker absent")
    script = script.replace(marker, hook + marker, 1)
    old_call = (
        "GaloisGroup(f : ShortOK := false, "
        f"ProofEffort := {args.proof_effort})"
    )
    new_call = (
        "GaloisGroup(f : Prime := qs[1], NextPrime := DeterministicNextPrime, "
        f"ShortOK := false, ProofEffort := {args.proof_effort})"
    )
    if old_call not in script:
        raise ValueError("generator layout changed: GaloisGroup call absent")
    script = script.replace(old_call, new_call, 1)
    args.output.write_text(script)

    receipt = {
        "status": "magma_multiseed_probe_emitted",
        "input": str(args.input.resolve()),
        "input_sha256": sha256(args.input),
        "output": str(args.output.resolve()),
        "output_sha256": sha256(args.output),
        "bytes": len(script.encode()),
        "proof_effort": args.proof_effort,
        "short_ok": False,
        "seed": 1,
        "prime_suite": [
            {"b": b, "c": c, "frobenius_order": order}
            for b, c, order in PRIMES
        ],
        "prime_suite_source": str(
            (HERE / "magma-quadratic-prime-scan.xml").resolve()
        ),
        "scope": (
            "A successful CERTIFIED line would prove M23; a timeout or any "
            "Magma error proves no group-theoretic upper bound."
        ),
    }
    receipt_path = HERE / "magma-function-field-multiseed-emission.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
