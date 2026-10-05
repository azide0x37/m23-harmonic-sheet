#!/usr/bin/env python3
"""Emit bounded Magma probes for the exact p=9319 degree-23 factor.

Input is the raw official-calculator XML receipt produced by the independent
split-prime exactification agent.  Output scripts are purely mechanical
transcriptions plus exact assertions; no mathematical claim is inferred from
calculator timeouts.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = (
    ROOT
    / "outputs/m23-proof-closeout-20260902/agents/petri-model-exactify"
    / "split-prime-p9319-a3983-s273-pair3-resultant-v2.xml"
)


def extract_polynomial() -> str:
    root = ET.parse(SOURCE).getroot()
    lines = [node.text or "" for node in root.findall(".//results/line")]
    begin = lines.index("MARKER DEFINING_POLYNOMIAL_BEGIN") + 1
    end = lines.index("MARKER DEFINING_POLYNOMIAL_END")
    polynomial = "\n".join(lines[begin:end]).strip()
    if not polynomial.startswith("w^23 +"):
        raise RuntimeError("unexpected polynomial prefix")
    return polynomial


def preamble(polynomial: str) -> str:
    return f'''p := 9319;
k := GF(p);
Kt<T> := FunctionField(k);
P<w> := PolynomialRing(Kt);
h := {polynomial};
assert Degree(h) eq 23;
assert IsMonic(h);
assert IsIrreducible(h);
assert IsSeparable(h);
print "EXACT_POLYNOMIAL_OK";
'''


def main() -> None:
    polynomial = extract_polynomial()
    base = preamble(polynomial)
    scan = base + r'''
patterns := AssociativeArray();
unramified := 0;
for ci in [0..p-1] do
    c := k!ci;
    if &and[Evaluate(Denominator(Coefficient(h,i)),c) ne 0 : i in [0..23]] then
        hc := Polynomial([Evaluate(Numerator(Coefficient(h,i)),c) /
                          Evaluate(Denominator(Coefficient(h,i)),c)
                          : i in [0..23]]);
        if Discriminant(hc) ne 0 then
            unramified +:= 1;
            ds := Sort([Degree(v[1]) : v in Factorization(hc)]);
            key := Sprint(ds);
            if not IsDefined(patterns,key) then
                patterns[key] := ci;
                printf "PATTERN c=%o degrees=%o\n", ci, ds;
            end if;
        end if;
    end if;
end for;
printf "SCAN_COMPLETE unramified=%o patterns=%o\n", unramified, #Keys(patterns);
'''
    (HERE / "magma-p9319-linear-prime-scan.m").write_text(scan)

    quadratic_scan = base + r'''
R0<U> := PolynomialRing(k);
patterns := AssociativeArray();
tested := 0;
unramified := 0;
for ci in [1..256] do
    q := U^2 + k!ci;
    if IsIrreducible(q) then
        tested +:= 1;
        E<b> := ext<k | q>;
        if &and[Evaluate(Denominator(Coefficient(h,i)),b) ne 0 : i in [0..23]] then
            PE<z> := PolynomialRing(E);
            hs := PE![Evaluate(Numerator(Coefficient(h,i)),b) /
                       Evaluate(Denominator(Coefficient(h,i)),b)
                       : i in [0..23]];
            if Degree(hs) eq 23 and IsSeparable(hs) then
                unramified +:= 1;
                ds := Sort([Degree(v[1]) : v in Factorization(hs)]);
                key := Sprint(ds);
                if not IsDefined(patterns,key) then
                    patterns[key] := ci;
                    ord := 1;
                    for d in ds do ord := Lcm(ord,d); end for;
                    printf "QPATTERN c=%o order=%o degrees=%o\n", ci, ord, ds;
                end if;
            end if;
        end if;
    end if;
end for;
printf "QSCAN_COMPLETE tested=%o unramified=%o patterns=%o\n",
       tested, unramified, #Keys(patterns);
'''
    (HERE / "magma-p9319-quadratic-prime-scan.m").write_text(quadratic_scan)

    template = base + r'''
SetSeed(1);
SetVerbose("GaloisGroup", 5);
SetVerbose("Invariant", 3);
prime := T - k!__C__;
printf "PRIME %o\n", prime;
assert &and[Evaluate(Denominator(Coefficient(h,i)),k!__C__) ne 0 : i in [0..23]];
hc := Polynomial([Evaluate(Numerator(Coefficient(h,i)),k!__C__) /
                  Evaluate(Denominator(Coefficient(h,i)),k!__C__)
                  : i in [0..23]]);
assert Discriminant(hc) ne 0;
print "PRIME_PATTERN", Sort([Degree(v[1]) : v in Factorization(hc)]);
procedure CertifiedRun()
    time G, roots, data := GaloisGroup(h : Prime := prime, ShortOK := false);
    print "GROUP_ORDER", #G;
    print "GROUP_DEGREE", Degree(G);
    print "GROUP_TRANSITIVE", IsTransitive(G);
    print "GROUP_GENERATORS", Generators(G);
    print "CERTIFIED_GALOISGROUP_RETURNED";
end procedure;
CertifiedRun();
'''
    (HERE / "magma-p9319-galois-template.m").write_text(template)
    seeds = {
        5: "1^5_3^6",
        740: "1^7_2^8",
        2: "23",
        1: "1_2_4_8^2",
    }
    seeded_outputs = []
    for c, label in seeds.items():
        filename = f"magma-p9319-galois-c{c}-{label}.m"
        (HERE / filename).write_text(template.replace("__C__", str(c)))
        seeded_outputs.append(filename)

    quadratic_template = base + r'''
SetSeed(1);
SetVerbose("GaloisGroup", 5);
SetVerbose("Invariant", 3);
R0<U> := PolynomialRing(k);
prime := U^2 + k!__C__;
assert IsIrreducible(prime);
E<b> := ext<k | prime>;
assert &and[Evaluate(Denominator(Coefficient(h,i)),b) ne 0 : i in [0..23]];
PE<z> := PolynomialRing(E);
hs := PE![Evaluate(Numerator(Coefficient(h,i)),b) /
          Evaluate(Denominator(Coefficient(h,i)),b) : i in [0..23]];
assert IsSeparable(hs);
print "QPRIME", prime;
print "QPRIME_PATTERN", Sort([Degree(v[1]) : v in Factorization(hs)]);
procedure CertifiedRun()
    time G, roots, data := GaloisGroup(h : Prime := prime, ShortOK := false);
    print "GROUP_ORDER", #G;
    print "GROUP_DEGREE", Degree(G);
    print "GROUP_TRANSITIVE", IsTransitive(G);
    print "GROUP_GENERATORS", Generators(G);
    print "CERTIFIED_GALOISGROUP_RETURNED";
end procedure;
CertifiedRun();
'''
    (HERE / "magma-p9319-galois-quadratic-template.m").write_text(
        quadratic_template
    )
    quadratic_seeds = {212: "1_5-3_6", 32: "23", 57: "1_3-2_2-4_4"}
    for c, label in quadratic_seeds.items():
        filename = f"magma-p9319-galois-q{c}-{label}.m"
        (HERE / filename).write_text(quadratic_template.replace("__C__", str(c)))
        seeded_outputs.append(filename)

    receipt = {
        "schema": "p9319-magma-probe-emission-v1",
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "polynomial_sha256": hashlib.sha256(polynomial.encode()).hexdigest(),
        "outputs": [
            "magma-p9319-linear-prime-scan.m",
            "magma-p9319-quadratic-prime-scan.m",
            "magma-p9319-galois-template.m",
            "magma-p9319-galois-quadratic-template.m",
            *seeded_outputs,
        ],
        "boundary": "emission only; calculator receipts are separate",
    }
    (HERE / "magma-p9319-probe-emission.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
