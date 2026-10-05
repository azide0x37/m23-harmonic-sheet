#!/usr/bin/env python3
"""Emit the exact canonical model over L = F(sqrt(-23)) as a GP file.

An element of L is written  Mod(Mod(m(y),f) + Mod(n(y),f)*x, x^2+23)/d  with
f = y^6-y^5-3y^4-3y^3+y^2+5y+4 and x = sqrt(-23).  Series work uses the variable z,
which has lower priority than x and y.   Usage: python3 emit_L.py K [out.gp]
"""
import sys, json
sys.path.insert(0, "/home/claude/m23/relq"); sys.path.insert(0, "/home/claude/m23/mon13")
from lift13 import REL

MARGIN_MIN = 20.0


def gp(r):
    if all(t == 0 for t in r["m"]) and all(t == 0 for t in r["n"]):
        return "L0"
    m = "+".join(f"({c})*y^{i}" for i, c in enumerate(r["m"]) if c) or "0"
    n = "+".join(f"({c})*y^{i}" for i, c in enumerate(r["n"]) if c) or "0"
    body = f"Mod({m},f)" if n == "0" else f"Mod({m},f)+Mod({n},f)*x"
    return f"Mod({body},x^2+23)" + (f"/({r['d']})" if r["d"] != 1 else "")


def main(K, out):
    rec = {r["label"]: r for r in json.load(open(f"{REL}/reconCanon-{K}.json"))["items"]}
    st = json.load(open(f"{REL}/canon-{K}.json"))
    need = [f"q{a}{b}" for a in range(4) for b in range(a, 4)]
    cub = ["C" + "".join(map(str, m)) for m in st["cubic_monomials"]]
    need += cub + [f"Pm{j}" for j in range(4)] + ["J_plus"]
    bad = [lab for lab in need if rec[lab]["margin"] < MARGIN_MIN]
    assert not bad, f"insufficient margin for {bad}"
    lines = [
        "\\\\ the exact canonical model of the M23 cover over L = F(sqrt(-23)),",
        f"\\\\ reconstructed from the 13-adic lift at precision 13^{K}.",
        "f = y^6 - y^5 - 3*y^4 - 3*y^3 + y^2 + 5*y + 4;",
        "L0 = Mod(Mod(0,f),x^2+23); L1 = Mod(Mod(1,f),x^2+23);",
        "SQM = Mod(Mod(1,f)*x, x^2+23);   \\\\ sqrt(-23)",
        f"KPREC = {K};",
    ]
    lines.append("QUAD = matrix(4,4);")
    for a in range(4):
        for b in range(a, 4):
            e = gp(rec[f"q{a}{b}"])
            lines.append(f"QUAD[{a+1},{b+1}] = {e};" + (f" QUAD[{b+1},{a+1}] = QUAD[{a+1},{b+1}];" if a != b else ""))
    lines.append("CUBM = " + str([list(m) for m in st["cubic_monomials"]]).replace(" ", "") + ";")
    lines.append("CUBC = [" + ", ".join(gp(rec[c]) for c in cub) + "];")
    lines.append("PM = [" + ", ".join(gp(rec[f"Pm{j}"]) for j in range(4)) + "];")
    lines.append("JPLUS_recognised = " + gp(rec["J_plus"]) + ";")
    lines.append("QUINTM = " + str([list(m) for m in st["quintic_monomials"]]).replace(" ", "") + ";")
    lines.append("STDQ = " + str(st["standard_quintics"]).replace(" ", "") + ";")
    open(out, "w").write("\n".join(lines) + "\n")
    print(f"written {out}: {len(lines)} lines")


if __name__ == "__main__":
    K = int(sys.argv[1])
    main(K, sys.argv[2] if len(sys.argv) > 2 else f"{REL}/model_hat-{K}.gp")
