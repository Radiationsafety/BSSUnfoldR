#!/usr/bin/env python
"""Check one method against every parity fixture and every synthetic spectrum.

Usage:
    python tools/parity/check_all.py unfold_maxed [--tol 1e-6] [--quick]

Runs the method through the R package and bssunfold on 4 fixtures
(main PTB, GSF, a deliberately non-uniform coarse grid, JINR) x 4 spectral
shapes. This is the anti-overfit gate: a port is only faithful if it holds on
grids the author did not tune against.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = {
    "main": "/tmp/bsscmp",
    "gsf": "/tmp/alt_gsf",
    "coarse": "/tmp/alt_coarse",
    "jinr": "/tmp/alt_jinr",
}
CASES = ["moderated", "fission", "evap", "split"]


# Seeded mode is the default: a stochastic port cannot be compared cell by
# cell against an unseeded reference, and unseeded Python-vs-Python already
# differs by rel_l2 ~0.8.  Methods whose R entry point has no random_state
# formal fall back to --no-rand for the whole run (cached in SEEDED).
SEEDED = {}


def _cell(method, fxdir, case, tol, seed):
    env = dict(os.environ, PARITY_OUT=fxdir)
    cmd = [sys.executable, str(ROOT / "tools" / "parity" / "one.py"),
           method, case, "--tol", str(tol)]
    if seed is None:
        cmd.append("--no-rand")
    else:
        cmd.extend(["--seed", str(seed)])
    p = subprocess.run(cmd, capture_output=True, text=True,
                       env=env, cwd=str(ROOT))
    out = p.stdout
    rel = None
    m = re.search(r"rel_l2\s*=\s*([-+0-9.eE]+)", out)
    if m:
        try:
            rel = float(m.group(1))
        except ValueError:
            rel = None
    if p.returncode == 2:
        side = "R" if "R  ERROR" in out else ("PY" if "PY ERROR" in out else "?")
        msg = [l for l in out.splitlines() if "ERROR" in l]
        return None, f"{side}-err: {(msg[0] if msg else '')[:70]}", out
    return rel, ("OK" if p.returncode == 0 else "DIVERGED"), out


def run_case(method, fxdir, case, tol, quick, seed=42):
    use_seed = SEEDED.get(method, seed)
    rel, status, out = _cell(method, fxdir, case, tol, use_seed)
    if (rel is None and status.startswith("R-err") and use_seed is not None
            and "random_state" in out):
        # the method does not accept a seed on the R side; compare unseeded
        SEEDED[method] = None
        rel, status, _ = _cell(method, fxdir, case, tol, None)
    SEEDED.setdefault(method, use_seed)
    return rel, status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("method")
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--quick", action="store_true",
                    help="only main fixture, moderated+fission")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    fixtures = {"main": FIXTURES["main"]} if args.quick else FIXTURES
    cases = ["moderated", "fission"] if args.quick else CASES

    worst = 0.0
    fails = 0
    print(f"{args.method}  (tol {args.tol:g})")
    for fname, fdir in fixtures.items():
        if not Path(fdir, "fixture.json").exists():
            print(f"  {fname}: fixture missing")
            continue
        for case in cases:
            rel, status = run_case(args.method, fdir, case, args.tol,
                                   args.quick, args.seed)
            if rel is not None:
                worst = max(worst, rel)
            if status not in ("OK",):
                fails += 1
            rels = "   n/a" if rel is None else f"{rel:.2e}"
            print(f"  {fname:<7} {case:<10} rel_l2={rels:>9}  {status}")
    print(f"=> worst rel_l2 = {worst:.3e}   failing cells = {fails}")
    sys.exit(0 if fails == 0 else 1)


if __name__ == "__main__":
    main()
