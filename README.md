# The Harmonic Sheet and the Sextic Orbit

Artifacts for **The Harmonic Sheet and the Sextic Orbit** (Alexander West
Templeton, October 2026), the second paper on the seven M23 Hurwitz points.
Paper DOI: *to be filled at Zenodo release*. Artifacts DOI: *to be filled*.
Part one: [The Galois action on the seven M23 Hurwitz points](https://doi.org/10.5281/zenodo.22073397),
artifacts [azide0x37/m23-seven-points](https://github.com/azide0x37/m23-seven-points).

## What the paper establishes

Huang, Jackson, Lee, Poonen, Pries and Zhang (arXiv:2608.08538) realize M23
over Q by a rational point in the seven-point Hurwitz fiber of covers of type
(2A, 23A, 23B). This paper determines the arithmetic of the fiber and the
geometry of the rational point.

- **At 13** (good reduction): an exact genus-4 cover over F_{13^6} with
  monodromy M23 and field of moduli exactly F_{13^6}; the tame Hurwitz
  problem extends to a finite étale scheme of rank 7 over Z_13. Hence the
  six nonrational points form one Galois orbit, a sextic field in which 13 is
  inert.
- **At 23** (bad reduction): the stable reduction of the rational cover is
  harmonic — tame point, two wild points and the attachment point of the
  unique new tail at 0, ±1, ∞ — with tails of conductors 7 and 15 and
  deformation datum c·z^6 dz/(z^22 − 1) on z^11 = t. One nonsquare cyclotomic
  character reverses the line and exchanges 23A, 23B; the harmonic datum has
  one lift per ordered class vector.
- **Over F(√−23)**, F = Q[x]/(x^6 − x^5 − 3x^4 − 3x^3 + x^2 + 5x + 4): an
  exact canonical model of one of the six covers, verified over that field,
  **including an exact certificate that its degree-23 map has exactly three
  branch points** (the eight-point test, new in this version), gives a point
  of the Hurwitz scheme at which the descended Petri invariant has minimal
  polynomial the recognized sextic. Hence E ≅ F: Galois closure S6,
  discriminant 2^4 · 11 · 23^4.

## Layout

| Path | Contents | Paper |
| --- | --- | --- |
| `paper/` | LaTeX source, ancillary table of the 49 coefficients of H', compiled PDFs | — |
| `relq/` | the exact canonical model over L (`model-hat.gp`, SHA-256 `7423ca42…`), the exact verification (V1)–(V8) (`run_verify.gp`, `verify_L.gp`, `check_V1_f13.py`, logs), and the 13-adic lift / recognition scripts that found it | §12.2, Thm 12.2 Steps 1–3, 6–8 |
| `eightpoint/` | the eight-point certificate (V9)–(V10): the plane relation H_0(τ₁, U₀) = 0 and the exact factorization H_0(b, U) = q₇ d₈² | Thm 12.2 Steps 4–5, Remark 12.5 |
| `petri-f13/` | the exact characteristic-13 model, the degree-23 map, passport, and the Petri value modulo 13 | §4, §12.1 (binding at 13) |
| `m23-sextic/` | the sextic Φ, the embedding r, Magma referee ledgers (copied from part one, MIT) | §13 |
| `mon13/`, `harmonic/`, `tails/`, `heptad/` | **being restored from archive** — the heptad certificate (Θ = Ψ₇₇Ψ₁₇₆), orientation, harmonic datum, Newton chambers, tails; see below | §5–§10 |

### Status of this repository

`relq/run_verify.gp` + `relq/check_V1_f13.py`, `eightpoint/`, `petri-f13/`
and `m23-sextic/` are complete and replay: every exact statement (V1)–(V10)
that Theorem 12.2 rests on is reproducible from this checkout. Still to be
restored from the 2026-09-29 archive (pre-release until then):

- the Part I/II receipt directories of Appendix D — `mon13/` (`heptad13.py`,
  `heptad-certificate.json` SHA-256 `0341e61a…`, `orient13.py`,
  `mon13-orientation.json` `83ceaf81…`, `hprime-canon.txt` `14bc240d…`,
  `m13a_canon.gp`, `elim_H.sing`, `canon_gaps.py`, `canon_2k.py`,
  `canon_E.sing`, `canon_L9.sing`, `cyclic_design.py`), `harmonic/`
  (`cartier.py`, `newton_chambers.py`, `inf_montes.gp`), `tails/`, `heptad/`;
- the outputs of the 13-adic lift that found the model, `relq/liftb-400.json`,
  `canon-400.json`, `reconCanon-400.json`, and the helper modules the lift
  scripts import (`wk.py`, `lift13.py`, `recon.py`, `petri13.py`,
  `petri13L.py`). These document how the model was found (Remark 12.4); the
  proof does not use them. `relq/check_V1.py` needs them; the replayable
  form of (V1) is `relq/check_V1_f13.py`, which compares `model-hat.gp`
  reduced at the place above 13 with `petri-f13/exact-f13-map.json`.

## Replay the exact model and the eight-point certificate

    pip install cypari2                      # PARI 2.17.x; gp for run_verify.gp and stage 1
    cd relq && gp -q --default parisizemax=4000000000 run_verify.gp   # (V2)–(V8), 3–7 min
    python3 check_V1_f13.py                   # (V1): reduction to the characteristic-13 model
    cd ../eightpoint && ./replay.sh           # (V9)–(V10), ~18 min on 2 cores

Expected: `run_verify.gp` prints `==> g_U = Phi : 1`; `check_V1_f13.py` prints
`(V1) verified`; `eightpoint/stage3.log` ends with `EIGHT-POINT TEST: PASS`
(exit code 1 otherwise).

## Verification layers

- **Exact over L** (the proof of Theorem 12.2): `relq/run_verify.gp` and
  `eightpoint/`. No height bound, no second prime, no floating point. The
  13-adic lift and the lattice recognition that *found* the model and the
  witness play no role in the proof.
- **Exact over F_{13^6}**: `petri-f13/` (Sage, python-flint, Magma calculator
  ledgers) and the heptad/orientation certificates (`mon13/`, pending
  restore).
- **Exact over F_23 / Q_23**: `harmonic/`, `tails/`, `heptad/` (pending
  restore).

## AI assistance

This project made extensive use of AI systems (Claude, Anthropic; ChatGPT,
OpenAI) for exploration, code, computation, analysis and drafting, under the
author's direction; the eight-point certificate and its library were written
and run in a Claude session on 2026-10-05. AI-generated suggestions were not
treated as mathematical evidence; every computational claim names a receipt
here.

## License

Code and receipts: MIT (see `LICENSE`). The paper (`paper/`) is CC BY 4.0.
`m23-sextic/` is copied from `azide0x37/m23-seven-points` (MIT).
