#!/usr/bin/env python3
"""Emit a bounded Magma scan for favourable irreducible quadratic primes."""

from __future__ import annotations

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "magma-function-field-galois.m"
OUTPUT = HERE / "magma-quadratic-prime-scan.m"


def main() -> None:
    source = SOURCE.read_text()
    marker = "time G, roots, data := GaloisGroup"
    if marker not in source:
        raise ValueError("GaloisGroup marker is missing")
    prefix = source.split(marker, 1)[0]
    program = prefix + r'''R0<T> := PolynomialRing(k);
tested := 0;
unramified := 0;
for bi, ci in [0..12] do
    q := T^2 + (a+bi)*T + (a+ci);
    if IsIrreducible(q) then
        tested +:= 1;
        E<b> := ext<k | q>;
        if &and[ Evaluate(Denominator(c),b) ne 0 : c in cs ] then
            PE<z> := PolynomialRing(E);
            fs := PE![ Evaluate(Numerator(c),b)/Evaluate(Denominator(c),b) : c in cs ];
            if Degree(fs) eq 23 and IsSeparable(fs) then
                unramified +:= 1;
                ff := Factorization(fs);
                degs := Sort(&cat[[Degree(pair[1]) : j in [1..pair[2]]] : pair in ff]);
                ord := 1;
                for d in degs do ord := Lcm(ord,d); end for;
                printf "PRIME b=%o c=%o order=%o degrees=%o\n", bi, ci, ord, degs;
            end if;
        end if;
    end if;
end for;
printf "SCAN_DONE tested=%o unramified=%o\n", tested, unramified;
end procedure;
CertifiedRun();
'''
    OUTPUT.write_text(program)
    receipt = {
        "status": "magma_quadratic_prime_scan_emitted",
        "source": str(SOURCE.resolve()),
        "output": str(OUTPUT.resolve()),
        "bytes": len(program.encode()),
        "candidate_family": "T^2+(a+b)T+(a+c), b,c in F13",
        "candidate_count": 169,
    }
    (HERE / "magma-quadratic-prime-scan-emission.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
