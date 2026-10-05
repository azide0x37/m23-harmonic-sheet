#!/usr/bin/env python3
"""Run an exact Sage/Singular passport probe in the pinned remote image.

The local host has no Sage executable, but the user-provided ``pc`` Docker
context already contains SageMath 10.9 (image id ``2401ffa8e9fc``).  This
driver deterministically translates the finite-field JSON certificate into a
self-contained Sage/Python program, sends it over stdin, and preserves the
program, stdout, stderr, and parsed JSON result in this agent directory.

No repository path is mounted into the remote container and the container is
removed after the bounded run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def coefficient(vector: list[int]) -> str:
    """Translate galois' high-to-low six-vector to a Sage field expression."""

    terms = []
    for index, value in enumerate(vector):
        exponent = 5 - index
        value %= 13
        if value == 0:
            continue
        if exponent == 0:
            terms.append(str(value))
        elif exponent == 1:
            terms.append(f"{value}*a")
        else:
            terms.append(f"{value}*a**{exponent}")
    return "(" + (" + ".join(terms) if terms else "0") + ")"


def polynomial(payload: dict[str, list[int]]) -> str:
    terms = []
    for monomial_text, vector in payload.items():
        monomial = tuple(
            int(value.strip()) for value in monomial_text.strip("()").split(",")
        )
        factors = [coefficient(vector), *(f"x{variable}" for variable in monomial)]
        terms.append("*".join(factors))
    return "(" + (" + ".join(terms) if terms else "0") + ")"


def sage_program(receipt: dict[str, Any]) -> str:
    q = polynomial(receipt["reduced_model"]["quadric"])
    c = polynomial(receipt["reduced_model"]["cubic"])
    n = polynomial(receipt["map"]["numerator"])
    d = polynomial(receipt["map"]["denominator"])
    return f'''from sage.all import *
import sage.version
import json
import time

started = time.monotonic()
Fp = GF(13)
Rz = PolynomialRing(Fp, "z")
z = Rz.gen()
g = z**6-z**5-3*z**4-3*z**3+z**2+5*z+4
k = GF(13**6, name="a", modulus=g)
a = k.gen()

def field_vector(value):
    coeffs = list(value.polynomial())
    coeffs += [Fp(0)] * (6-len(coeffs))
    return [int(coeffs[5-i]) for i in range(6)]

def hp_text(ideal):
    return str(ideal.hilbert_polynomial())

R4 = PolynomialRing(k, names=("x0","x1","x2","x3"), order="degrevlex")
x0,x1,x2,x3 = R4.gens()
Q = R4({q})
C = R4({c})
N = R4({n})
D = R4({d})
I = R4.ideal(Q,C)
irrelevant = R4.ideal(x0,x1,x2,x3)
print("curve ideal ready", flush=True)

Braw = I + R4.ideal(N,D)
B, Bexp = Braw.saturation(irrelevant)
print("base saturation", Bexp, hp_text(B), flush=True)

JN, _ = (I + R4.ideal(N)).saturation(irrelevant)
JD, _ = (I + R4.ideal(D)).saturation(irrelevant)
ZN, ZNexp = JN.saturation(B)
ZD, ZDexp = JD.saturation(B)
ZN, _ = ZN.saturation(irrelevant)
ZD, _ = ZD.saturation(irrelevant)
print("zero/pole residual", hp_text(ZN), hp_text(ZD), flush=True)

# Ramification at finite nonzero target values.  The 3x3 minors say that the
# pencil section N-tD is tangent to the smooth complete intersection.
R5 = PolynomialRing(k, names=("x0","x1","x2","x3","t"), order="degrevlex")
x0,x1,x2,x3,t = R5.gens()
Q5,C5,N5,D5 = map(R5,(Q,C,N,D))
F5 = N5-t*D5
variables = (x0,x1,x2,x3)
rows = ([Q5.derivative(v) for v in variables],
        [C5.derivative(v) for v in variables],
        [F5.derivative(v) for v in variables])
minors = []
for omitted in range(4):
    cols = [j for j in range(4) if j != omitted]
    minors.append(matrix(R5, [[rows[i][j] for j in cols] for i in range(3)]).det())
J = R5.ideal(Q5,C5,F5,*minors)
m5 = R5.ideal(x0,x1,x2,x3)
B5 = R5.ideal(Q5,C5,N5,D5)
J, Jirr_exp = J.saturation(m5)
print("ramification irrelevant saturation", Jirr_exp, flush=True)
J, Jbase_exp = J.saturation(B5)
print("ramification base saturation", Jbase_exp, flush=True)
J, Jzero_exp = J.saturation(R5.ideal(t))
print("ramification zero-value saturation", Jzero_exp, flush=True)
E = J.elimination_ideal([x0,x1,x2,x3])
print("elimination", E, flush=True)
egens = [p for p in E.gens() if p != 0]

result = {{
  "status": "sage_exact_passport_probe",
  "sage_version": sage.version.version,
  "field": {{"characteristic": 13, "degree": 6}},
  "curve_ideal_dimension": int(I.dimension()),
  "base": {{
    "saturation_exponent": int(Bexp),
    "dimension": int(B.dimension()),
    "hilbert_polynomial": hp_text(B),
    "groebner_basis_size": len(B.groebner_basis()),
  }},
  "zero_fibre_after_base_removal": {{
    "saturation_exponent": int(ZNexp),
    "dimension": int(ZN.dimension()),
    "hilbert_polynomial": hp_text(ZN),
    "radical_hilbert_polynomial": hp_text(ZN.radical()),
  }},
  "pole_fibre_after_base_removal": {{
    "saturation_exponent": int(ZDexp),
    "dimension": int(ZD.dimension()),
    "hilbert_polynomial": hp_text(ZD),
    "radical_hilbert_polynomial": hp_text(ZD.radical()),
  }},
  "finite_nonzero_ramification": {{
    "irrelevant_saturation_exponent": int(Jirr_exp),
    "base_saturation_exponent": int(Jbase_exp),
    "zero_value_saturation_exponent": int(Jzero_exp),
    "ideal_dimension": int(J.dimension()),
    "elimination_generators": [str(p) for p in egens],
  }},
  "elapsed_seconds": time.monotonic()-started,
}}
if len(egens) == 1 and egens[0].degree(t) == 1:
    poly = egens[0]
    lam = k(-poly.subs({{t:0}}) / poly.monomial_coefficient(t))
    result["finite_nonzero_ramification"]["branch_value"] = field_vector(lam)
    result["finite_nonzero_ramification"]["branch_value_text"] = str(lam)
    Slambda = N-lam*D
    Flambda, _ = (I + R4.ideal(Slambda)).saturation(irrelevant)
    Fibre, fibre_base_exp = Flambda.saturation(B)
    Fibre, _ = Fibre.saturation(irrelevant)
    Radical = Fibre.radical()
    radical_square = I + R4.ideal([
        left*right for left in Radical.gens() for right in Radical.gens()
    ])
    square_contained = bool(radical_square <= Fibre)
    fibre_hp = Fibre.hilbert_polynomial()
    radical_hp = Radical.hilbert_polynomial()
    result["third_fibre"] = {{
      "base_saturation_exponent": int(fibre_base_exp),
      "dimension": int(Fibre.dimension()),
      "hilbert_polynomial": str(fibre_hp),
      "radical_hilbert_polynomial": str(radical_hp),
      "radical_square_contained_in_fibre_ideal": square_contained,
      "fibre_ideal_contained_in_radical": bool(Fibre <= Radical),
      "length": int(fibre_hp),
      "radical_length": int(radical_hp),
      "excess_length": int(fibre_hp-radical_hp),
    }}
    if int(fibre_hp) == 23 and int(radical_hp) == 15 and square_contained:
        result["third_fibre"].update({{
          "multiplicity_partition_geometric": [2]*8 + [1]*7,
          "ramification_contribution": 8,
          "passport_certified": True,
        }})
print("RESULT_JSON=" + json.dumps(result, sort_keys=True), flush=True)
'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=Path, default=HERE / "exact-f13-map.json"
    )
    parser.add_argument(
        "--image", default="2401ffa8e9fc", help="pinned SageMath 10.9 image id"
    )
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument(
        "--program-output", type=Path, default=HERE / "sage-passport-probe-generated.py"
    )
    parser.add_argument(
        "--output", type=Path, default=HERE / "sage-passport-probe.json"
    )
    parser.add_argument(
        "--stdout", type=Path, default=HERE / "sage-passport-probe.stdout"
    )
    parser.add_argument(
        "--stderr", type=Path, default=HERE / "sage-passport-probe.stderr"
    )
    args = parser.parse_args()

    receipt = json.loads(args.input.read_text())
    if receipt["status"] != "exact_f13_degree23_candidate_constructed":
        raise ValueError(f"map input is not accepted: {receipt['status']}")
    program = sage_program(receipt)
    args.program_output.write_text(program)
    completed = subprocess.run(
        [
            "docker", "--context", "pc", "run", "--rm", "-i",
            "--entrypoint", "sage", args.image, "-python", "-",
        ],
        input=program,
        text=True,
        capture_output=True,
        timeout=args.timeout,
        check=False,
    )
    args.stdout.write_text(completed.stdout)
    args.stderr.write_text(completed.stderr)
    marker = "RESULT_JSON="
    lines = [line for line in completed.stdout.splitlines() if line.startswith(marker)]
    result: dict[str, Any] = {
        "status": "sage_probe_failed",
        "runner": {
            "returncode": completed.returncode,
            "image_id": args.image,
            "timeout_seconds": args.timeout,
            "input_path": str(args.input.resolve()),
            "input_sha256": sha256(args.input),
            "generated_program_path": str(args.program_output.resolve()),
            "generated_program_sha256": sha256(args.program_output),
            "stdout_path": str(args.stdout.resolve()),
            "stderr_path": str(args.stderr.resolve()),
        },
    }
    if completed.returncode == 0 and len(lines) == 1:
        result.update(json.loads(lines[0][len(marker):]))
        result["runner"] = {
            **result.get("runner", {}),
            **{
                "returncode": completed.returncode,
                "image_id": args.image,
                "timeout_seconds": args.timeout,
                "input_path": str(args.input.resolve()),
                "input_sha256": sha256(args.input),
                "generated_program_path": str(args.program_output.resolve()),
                "generated_program_sha256": sha256(args.program_output),
                "stdout_path": str(args.stdout.resolve()),
                "stderr_path": str(args.stderr.resolve()),
            },
        }
    else:
        result["failure"] = {
            "result_marker_count": len(lines),
            "stderr_tail": completed.stderr[-4000:],
            "stdout_tail": completed.stdout[-4000:],
        }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(args.output)
    print("returncode", completed.returncode, "markers", len(lines))
    print(completed.stdout[-4000:])
    if completed.stderr:
        print(completed.stderr[-4000:])


if __name__ == "__main__":
    main()
