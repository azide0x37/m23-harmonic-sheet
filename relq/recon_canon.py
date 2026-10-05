#!/usr/bin/env python3
"""LLL reconstruction over L = F(sqrt(-23)) of the canonical model written by canon13.py.

For each element alpha of W_K the lattice {(d, m_0..m_5, n_0..n_5) : d*alpha = sum m_i xhat^i
+ shat * sum n_i xhat^i  mod 13^K} has covolume 13^{6K} in dimension 13, so a random element
has shortest vector of size about 13^{6K/13}; a recognised element of height 13^h is found
with margin 6K/13 - h.   Usage: python3 recon_canon.py K
"""
import sys, json, math
sys.path.insert(0, "/home/claude/m23/relq"); sys.path.insert(0, "/home/claude/m23/mon13")
from wk import Wk
from recon import shat, reconstruct, embed_check
from lift13 import log, REL


def rec(W, K, s, asc, label, out):
    d, m, n, _ = reconstruct(asc, K, s, W=W)
    assert embed_check(W, d, m, n, s, W.from_asc(asc)), f"{label}: embedding check failed"
    nrm2 = d * d + sum(t * t for t in m) + sum(t * t for t in n)
    h = (math.log(nrm2 >> max(0, nrm2.bit_length() - 60)) +
         max(0, nrm2.bit_length() - 60) * math.log(2)) / (2 * math.log(13))
    out.append({"label": label, "d": d, "m": m, "n": n, "h": h, "margin": 6 * K / 13 - h})
    return out[-1]


def main(K):
    W = Wk(K); s = shat(K)
    st = json.load(open(f"{REL}/canon-{K}.json"))
    items = []
    for key, v in sorted(st["qhat"].items()):
        items.append((f"q{key.replace(',', '')}", v))
    cub = [tuple(m) for m in st["cubic_monomials"]]
    for m, v in zip(cub, st["Chat"]):
        items.append(("C" + "".join(map(str, m)), v))
    quint = [tuple(m) for m in st["quintic_monomials"]]
    nonpiv = st["standard_quintics"]
    for i, v in zip(nonpiv, st["Nhat"]):
        items.append(("N" + "".join(map(str, quint[i])), v))
    for i, v in zip(nonpiv, st["Dhat"]):
        items.append(("D" + "".join(map(str, quint[i])), v))
    for j, v in enumerate(st["Pminus"]):
        items.append((f"Pm{j}", v))
    items.append(("kappa_hat", st["kappa_hat"]))
    items.append(("mu", st["mu"]))
    items.append(("J_plus", st["J_plus"]))
    out = []
    groups = {}
    for label, asc in items:
        r = rec(W, K, s, asc, label, out)
        g = label[0] if label[0] in "qCND" else label
        gg = groups.setdefault(g, [1e9, -1e9, None])
        if r["margin"] < gg[0]: gg[0] = r["margin"]; gg[2] = label
        gg[1] = max(gg[1], r["h"])
    log(f"expected random shortest vector: 13^{6*K/13:.1f}")
    for g, (worst, maxh, lab) in sorted(groups.items()):
        log(f"  {g:>10}: max height 13^{maxh:.1f} ({int(maxh*math.log10(13))+1} digits), "
            f"worst margin 13^{worst:.1f} at {lab}")
    json.dump({"K": K, "shat": s, "items": out}, open(f"{REL}/reconCanon-{K}.json", "w"))
    log(f"written reconCanon-{K}.json")


if __name__ == "__main__":
    main(int(sys.argv[1]))
