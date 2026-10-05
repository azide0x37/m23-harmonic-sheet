#!/bin/bash
# Replay the eight-point certificate (V9), (V10). Needs gp (PARI 2.17) and python3 + cypari2.
set -euo pipefail
cd "$(dirname "$0")"
gp -q --default parisizemax=6000000000 stage1-ND.gp > stage1.log
NPP2=160  python3 stage2.py
NPREC=2500 python3 stage3.py
tail -2 stage3.log
