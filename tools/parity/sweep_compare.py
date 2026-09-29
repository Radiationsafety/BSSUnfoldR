#!/usr/bin/env python
"""Aggregate parity across every fixture and every synthetic spectrum.

Reads {fixture}/{py,r}/<case>.json produced by py_runner.py and r_runner.R and
writes one CSV with a rel_l2 per (fixture, case, method) cell.

Usage: python tools/parity/sweep_compare.py --out /tmp/parity_matrix.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


def rel_l2(xr, xp):
    xr = np.asarray(xr, float)
    xp = np.asarray(xp, float)
    n = min(len(xr), len(xp))
    xr, xp = xr[:n], xp[:n]
    d = np.linalg.norm(xp)
    if not np.isfinite(xr).all() or not np.isfinite(xp).all():
        return None, "nonfinite"
    if d == 0:
        return None, "py-zero"
    return float(np.linalg.norm(xr - xp) / d), "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixtures", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=float, default=1e-6)
    args = ap.parse_args()

    rows = []
    for fdir in args.fixtures:
        fdir = Path(fdir)
        tag = fdir.name
        pyd, rd = fdir / "py", fdir / "r"
        if not pyd.is_dir():
            print("missing py results for", tag, file=sys.stderr)
            continue
        # Only the known synthetic spectra; per-method scratch dumps an agent
        # may have left in these directories must not be read as cases.
        cases = [c for c in ("moderated", "fission", "evap", "split")
                 if (pyd / f"{c}.json").exists()]
        for case in cases:
            case_file = pyd / f"{case}.json"
            py = json.loads(case_file.read_text())
            rf = rd / case_file.name
            rr = json.loads(rf.read_text()) if rf.exists() else {}
            for m in sorted(set(py) | set(rr)):
                p, r = py.get(m), rr.get(m)
                # jsonlite writes an empty R list as [], not {}
                if not isinstance(p, dict):
                    p = {"error": "empty"} if p is not None else p
                if not isinstance(r, dict):
                    r = {"error": "empty"} if r is not None else r
                if p is None:
                    status = "missing-in-py"
                elif r is None:
                    status = "missing-in-r"
                elif "error" in p and "error" in r:
                    status = "err-both"
                elif "error" in p:
                    status = "py-error"
                elif "error" in r:
                    status = "r-error"
                else:
                    val, note = rel_l2(r["spectrum"], p["spectrum"])
                    status = note if note != "ok" else (
                        "ok" if val is not None and val <= args.tol else "diverged")
                    p = {**p, "rel_l2": val}
                rows.append({"fixture": tag, "case": case, "method": m,
                             "status": status,
                             "rel_l2": (p or {}).get("rel_l2", "")})

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["fixture", "case", "method",
                                           "status", "rel_l2"])
        w.writeheader()
        w.writerows(rows)

    # summary per method
    methods = sorted({r["method"] for r in rows})
    ncells = {}
    for m in methods:
        sub = [r for r in rows if r["method"] == m]
        ncells[m] = (sum(1 for r in sub if r["status"] == "ok"), len(sub))

    from collections import Counter
    print(Counter(r["status"] for r in rows).most_common())
    print(f"\n{'method':<36} {'ok/total':>9}")
    print("-" * 48)
    for m in methods:
        ok, tot = ncells[m]
        mark = "  " if ok == tot else "X "
        print(f"{mark}{m:<34} {ok:>4}/{tot:<4}")
    print("-" * 48)
    total = sum(t for _, t in ncells.values())
    good = sum(o for o, _ in ncells.values())
    print(f"cells ok: {good}/{total}   fully-clean methods: "
          f"{sum(1 for o, t in ncells.values() if o == t)}/{len(ncells)}")
    print("csv:", args.out, file=sys.stderr)


if __name__ == "__main__":
    main()
