#!/usr/bin/env python3
"""Extract the exact degree-23 function-field polynomial for the F_13^6 map.

This is a bounded, deterministic companion to ``sage_passport_probe.py``.
It sends a self-contained exact Sage program to the already pinned SageMath
10.9 Docker image on the user-provided ``pc`` context.  The program eliminates
the canonical complete intersection on the affine chart ``x3=1``, removes the
fixed base/chart factors by factorisation over ``k(t)``, and records every
factor.  No repository path is mounted remotely.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from sage_passport_probe import polynomial


HERE = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sage_program(receipt: dict[str, Any]) -> str:
    q = polynomial(receipt["reduced_model"]["quadric"])
    c = polynomial(receipt["reduced_model"]["cubic"])
    n = polynomial(receipt["map"]["numerator"])
    d = polynomial(receipt["map"]["denominator"])
    return f'''from sage.all import *
import json

Fp = GF(13)
Ru = PolynomialRing(Fp, "u")
u = Ru.gen()
g = u**6-u**5-3*u**4-3*u**3+u**2+5*u+4
k = GF(13**6, name="a", modulus=g)
a = k.gen()

def fv(value):
    coeffs = list(k(value).polynomial())
    coeffs += [Fp(0)] * (6-len(coeffs))
    return [int(coeffs[5-i]) for i in range(6)]

R4 = PolynomialRing(k, names=("x0","x1","x2","x3"))
x0,x1,x2,x3 = R4.gens()
Q = R4({q})
C = R4({c})
N = R4({n})
D = R4({d})

Sx = PolynomialRing(k, names=("xx","y","z"))
xx,yx,zx = Sx.gens()
aff = R4.hom([xx,yx,zx,Sx(1)], Sx)
Qaff = aff(Q)
S = PolynomialRing(k, names=("y","z"))
y,z = S.gens()
drop = Sx.hom([S(0),y,z], S)
B = drop(Qaff.subs({{xx:Sx(0)}}))
A = drop(Qaff.subs({{xx:Sx(1)}}) - Qaff.subs({{xx:Sx(0)}}))
assert Qaff.degree(xx) == 1 and A != 0
FS = S.fraction_field()
sub = R4.hom([FS(-B/A),FS(y),FS(z),FS(1)], FS)

def substitute_clear(poly, power):
    value = sub(poly)
    return S(value * FS(A**power))

# The exact support shows degree_x0(C)<=2 and degree_x0(N),degree_x0(D)<=5.
assert C.degree(x0) <= 2 and N.degree(x0) <= 5 and D.degree(x0) <= 5
F = substitute_clear(C, 2)
Np = substitute_clear(N, 5)
Dp = substitute_clear(D, 5)
assert F != 0 and Np != 0 and Dp != 0

R = PolynomialRing(k, names=("tt","yy","zz"))
tt,yy,zz = R.gens()
FR = R(F(y=yy,z=zz))
NR = R(Np(y=yy,z=zz))
DR = R(Dp(y=yy,z=zz))
res = FR.resultant(NR-tt*DR, yy)
assert res != 0

Kt = FunctionField(k, "t")
t = Kt.gen()
P = PolynomialRing(Kt, "w")
w = P.gen()

def to_Kt(poly_t):
    # Coefficients are in k and tt is mapped to the function-field generator.
    return Kt(sum(k(coef) * t**exp[0] for exp,coef in poly_t.dict().items()))

resP = P([to_Kt(res.polynomial(zz)[i]) for i in range(res.degree(zz)+1)])
fac = resP.factor()
factors = []
for h,m in fac:
    hm = h.monic()
    factors.append({{
        "degree": int(hm.degree()),
        "multiplicity": int(m),
        "coefficients": [
            {{
                "numerator": [fv(v) for v in c.numerator().list()],
                "denominator": [fv(v) for v in c.denominator().list()],
            }}
            for c in hm.list()
        ],
    }})

out = {{
    "status": "exact_function_field_elimination",
    "field": {{"characteristic":13,"degree":6,"modulus":[1,12,10,10,1,5,4]}},
    "chart": {{
        "quadric_linear_coefficient_degree": [int(A.degree(y)),int(A.degree(z))],
        "quadric_remainder_degree": [int(B.degree(y)),int(B.degree(z))],
        "curve_degree": [int(F.degree(y)),int(F.degree(z))],
        "numerator_degree": [int(Np.degree(y)),int(Np.degree(z))],
        "denominator_degree": [int(Dp.degree(y)),int(Dp.degree(z))],
    }},
    "raw_resultant": {{
        "degree_t": int(res.degree(tt)),
        "degree_source": int(res.degree(zz)),
        "factor_count": len(factors),
    }},
    "factors": factors,
}}
out["degree23_factor_count"] = sum(1 for h in factors if h["degree"] == 23)
out["accepted"] = out["degree23_factor_count"] == 1
print("RESULT_JSON=" + json.dumps(out, sort_keys=True))
'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=HERE / "exact-f13-map.json")
    parser.add_argument("--image", default="2401ffa8e9fc")
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument(
        "--program-output", type=Path,
        default=HERE / "sage-function-field-polynomial-generated.py",
    )
    parser.add_argument(
        "--output", type=Path, default=HERE / "function-field-polynomial.json"
    )
    parser.add_argument(
        "--stdout", type=Path,
        default=HERE / "sage-function-field-polynomial.stdout",
    )
    parser.add_argument(
        "--stderr", type=Path,
        default=HERE / "sage-function-field-polynomial.stderr",
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
    wrapper: dict[str, Any] = {
        "status": "sage_function_field_extraction_failed",
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
        wrapper.update(json.loads(lines[0][len(marker):]))
    else:
        wrapper["diagnostic"] = {
            "result_marker_count": len(lines),
            "stdout_tail": completed.stdout[-4000:],
            "stderr_tail": completed.stderr[-4000:],
        }
    args.output.write_text(json.dumps(wrapper, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": wrapper["status"],
        "accepted": wrapper.get("accepted", False),
        "output": str(args.output),
    }, indent=2, sort_keys=True))
    if wrapper.get("accepted") is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
