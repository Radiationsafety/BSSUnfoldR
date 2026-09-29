#!/usr/bin/env python
"""Compare R vs bssunfold spectra method-by-method.

Reports, per method and case:
  status  both / r-error / py-error
  rel_l2  ||x_r - x_py|| / ||x_py||
  cos     cosine similarity
  int_rat sum(x_r) / sum(x_py)
Exit code lists methods above --tol (default 1e-3).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def metrics(xr, xp):
    xr = np.asarray(xr, float)
    xp = np.asarray(xp, float)
    n = min(len(xr), len(xp))
    same_len = len(xr) == len(xp)
    xr, xp = xr[:n], xp[:n]
    denom = np.linalg.norm(xp)
    rel = np.linalg.norm(xr - xp) / denom if denom else float("nan")
    nr = np.linalg.norm(xr)
    cos = float(xr @ xp / (nr * denom)) if nr and denom else float("nan")
    sr = float(xr.sum())
    sp = float(xp.sum())
    # pointwise worst relative deviation on bins carrying real signal
    mask = xp > 1e-6 * (sp / n if sp else 1.0)
    pw = float(np.max(np.abs(xr[mask] - xp[mask]) / xp[mask])) if mask.any() else float("nan")
    return {
        "len_ok": same_len,
        "rel_l2": float(rel),
        "cos": cos,
        "int_ratio": sr / sp if sp else float("nan"),
        "pointwise": pw,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--py", default="/tmp/bsscmp/py")
    ap.add_argument("--r", default="/tmp/bsscmp/r")
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--show", choices=["bad", "all", "errors"], default="bad")
    ap.add_argument("--csv", default="")
    args = ap.parse_args()

    rows = []
    for case_file in sorted(Path(args.py).glob("*.json")):
        case = case_file.stem
        py = json.loads(case_file.read_text())
        rp = Path(args.r) / case_file.name
        rr = json.loads(rp.read_text()) if rp.exists() else {}
        for m in sorted(set(py) | set(rr)):
            p, r = py.get(m), rr.get(m)
            row = {"case": case, "method": m}
            if p is None:
                row["status"] = "missing-in-py"
            elif r is None:
                row["status"] = "missing-in-r"
            elif "error" in p and "error" in r:
                row["status"] = "err-both"
                row["py_err"] = p["error"][:80]
                row["r_err"] = r["error"][:80]
            elif "error" in p:
                row["status"] = "py-error"
                row["py_err"] = p["error"][:80]
            elif "error" in r:
                row["status"] = "r-error"
                row["r_err"] = r["error"][:80]
            else:
                row["status"] = "both"
                row.update(metrics(r["spectrum"], p["spectrum"]))
            rows.append(row)

    hdr = ["case", "method", "status", "rel_l2", "cos", "int_ratio",
           "pointwise", "len_ok", "r_err", "py_err"]
    print(f"{'case':<10} {'method':<34} {'status':<13} {'rel_l2':>9} "
          f"{'cos':>8} {'int_rat':>8} {'ptw':>9}")
    print("-" * 104)
    for row in rows:
        st = row["status"]
        if args.show == "bad" and st == "both" and row.get("rel_l2", 0) <= args.tol:
            continue
        if args.show == "errors" and not st.endswith("error") and st not in (
                "missing-in-py", "missing-in-r", "err-both", "r-error"):
            if st == "both" and row.get("rel_l2", 0) <= args.tol:
                continue
        f = lambda k: (f"{row[k]:.4g}" if k in row and isinstance(row[k], float) else "-")
        print(f"{row['case']:<10} {row['method']:<34} {st:<13} "
              f"{f('rel_l2'):>9} {f('cos'):>8} {f('int_ratio'):>8} {f('pointwise'):>9}")

    both = [r for r in rows if r["status"] == "both"]
    print("-" * 104)
    print(f"compared={len(both)}  within_tol={sum(1 for r in both if r['rel_l2'] <= args.tol)}  "
          f"diverged={sum(1 for r in both if r['rel_l2'] > args.tol)}  "
          f"r_error={sum(1 for r in rows if r['status'] in ('r-error','err-both','missing-in-r'))}")
    if args.csv:
        import csv
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=hdr, extrasaction="ignore")
            w.writeheader()
            for row in rows:
                w.writerow(row)
        print("csv:", args.csv, file=sys.stderr)


if __name__ == "__main__":
    main()
