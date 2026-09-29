#!/usr/bin/env python
"""Re-score the cached parity sweep with rank-aware verdicts.

`sweep.sh` compares raw spectra and calls anything above `rel_l2 > tol` a
divergence.  That over-reports: spectrum-level rel_l2 is only a meaningful
parity gate when the response matrix has full column rank, and several
fixtures violate that (the JINR spheres in /tmp/alt_jinr give six *identical*
response columns for the lowest bins, so any allocation of flux among them
reproduces the readings exactly and the vertex choice is arbitrary).

This script reads the spectra both languages already cached under
`<fixture>/r/<case>.json` and `<fixture>/py/<case>.json` -- it does not re-run
anything -- and re-scores each cell on three independent quantities:

  spectrum    rel_l2 between the two exported spectra
  collapsed   rel_l2 after summing each group of numerically identical
              response columns into one bin, i.e. the part of the difference
              that no measurement can constrain
  fit         ||A x - b|| / (||A||_2 ||b||), recomputed here from each side's
              exported spectrum, so neither side is grading its own homework

Verdicts: OK (collapsed within --tol), DEGENERATE (collapsed within --tol but
the raw spectra differ, i.e. a pure null-space tie), NOISE (both fits and both
collapsed spectra agree to within solver tolerance), DIVERGED (a real
disagreement), plus the runnability states the sweep already recorded.

Usage:
    python tools/parity/reclassify.py [--fixtures ...] [--tol 1e-4]
                                      [--out /tmp/parity_degen.csv]
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

CASES = ("moderated", "fission", "evap", "split")


def rel(a, b):
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    n = min(len(a), len(b))
    if n == 0:
        return float("nan")
    return float(np.linalg.norm(a[:n] - b[:n]) / max(np.linalg.norm(b[:n]), 1e-300))


def identical_column_groups(M, tol=1e-10):
    """Columns of the design matrix that no measurement can tell apart."""
    norm = np.linalg.norm(M, axis=0)
    parent = list(range(M.shape[1]))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(M.shape[1]):
        for j in range(i + 1, M.shape[1]):
            if norm[i] == 0 or norm[j] == 0:
                continue
            if np.linalg.norm(M[:, i] - M[:, j]) / max(norm[i], norm[j]) <= tol:
                parent[find(i)] = find(j)
    groups: dict[int, list[int]] = {}
    for i in range(M.shape[1]):
        groups.setdefault(find(i), []).append(i)
    return [g for g in groups.values() if len(g) > 1]


def collapse(x, groups):
    x = np.asarray(x, float).ravel().copy()
    for g in groups:
        total = float(sum(x[i] for i in g))
        x[g[0]] = total
        for i in g[1:]:
            x[i] = 0.0
    return x


def design_matrix(fixture: Path):
    """The solver's design matrix, (n_detectors, n_bins), ln_steps included.

    Taken from bssunfold's own Detector so the lethargy convention cannot
    drift from what both runners actually handed their solvers.
    """
    from bssunfold import Detector

    fx = json.loads((fixture / "fixture.json").read_text())
    rf = {"E_MeV": fx["E_MeV"]}
    for name, sens in zip(fx["detector_names"], fx["sensitivities"]):
        rf[name] = sens
    det = Detector(response_functions=rf)
    return np.asarray(det.Amat, float).T, fx


def load_side(fixture: Path, side: str):
    out = {}
    for case in CASES:
        p = fixture / side / f"{case}.json"
        if not p.exists():
            continue
        out[case] = json.loads(p.read_text())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixtures", nargs="+",
                    default=["/tmp/bsscmp", "/tmp/alt_gsf",
                             "/tmp/alt_coarse", "/tmp/alt_jinr"])
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--out", default="/tmp/parity_degen.csv")
    args = ap.parse_args()

    rows = []
    for fixture in map(Path, args.fixtures):
        M, fx = design_matrix(fixture)
        groups = identical_column_groups(M)
        scale = max(np.linalg.norm(M, 2) *
                    max(np.linalg.norm(np.asarray(fx["cases"][c]["readings"], float))
                        for c in CASES), 1e-300)
        rside, pyside = load_side(fixture, "r"), load_side(fixture, "py")
        methods = sorted(rside.get("moderated", {}))
        for method in methods:
            for case in CASES:
                rec_r = rside.get(case, {}).get(method)
                rec_p = pyside.get(case, {}).get(method)
                b = np.asarray(fx["cases"][case]["readings"], float)
                r_err = rec_r is None or bool(rec_r.get("error"))
                p_err = rec_p is None or bool(rec_p.get("error"))
                if r_err and p_err:
                    # Both languages refuse the same input: that is parity.
                    status = "err-both"
                    detail = str((rec_r or {}).get("error", "missing"))[:80]
                elif r_err:
                    status = "r-error"
                    detail = str((rec_r or {}).get("error", "missing"))[:120]
                elif p_err:
                    status = "py-error" if rec_p is not None else "py-missing"
                    detail = str((rec_p or {}).get("error", ""))[:120]
                else:
                    xp = np.asarray(rec_p["spectrum"], float)
                    xr = np.asarray(rec_r["spectrum"], float)
                    if np.allclose(xp, 0):
                        status, detail = "py-zero", ""
                    else:
                        spectrum = rel(xr, xp)
                        collapsed = rel(collapse(xr, groups), collapse(xp, groups))
                        n = min(M.shape[1], len(xp), len(xr))
                        fp = np.linalg.norm(M[:, :n] @ xp[:n] - b) / scale
                        fr = np.linalg.norm(M[:, :n] @ xr[:n] - b) / scale
                        dose = float("nan")
                        try:
                            dk = list(rec_p["doserates"])
                            dose = rel([rec_r["doserates"][k] for k in dk],
                                       [rec_p["doserates"][k] for k in dk])
                        except Exception:  # noqa: BLE001
                            pass
                        if collapsed <= args.tol:
                            status = ("DEGENERATE" if spectrum > args.tol
                                      else "NOISE" if spectrum > 1e-5 else "OK")
                            if status == "NOISE" and max(fp, fr) > 1e-6:
                                status = "DIVERGED"
                        else:
                            status = "DIVERGED"
                        detail = (f"spectrum={spectrum:.3e} collapsed={collapsed:.3e} "
                                  f"fit(py={fp:.1e},r={fr:.1e}) dose={dose:.3e}")
                rows.append({"fixture": fixture.name, "case": case,
                             "method": method, "status": status,
                             "detail": detail})

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=("fixture", "case", "method",
                                           "status", "detail"))
        w.writeheader()
        w.writerows(rows)

    tally: dict[str, int] = {}
    for r in rows:
        tally[r["status"]] = tally.get(r["status"], 0) + 1
    print("cell tally:", "  ".join(f"{k}={v}" for k, v in sorted(tally.items())))
    print("degenerate response columns: " + ", ".join(
        f"{Path(f).name}={identical_column_groups(design_matrix(Path(f))[0])}"
        for f in args.fixtures))

    per: dict[str, dict[str, int]] = {}
    for r in rows:
        per.setdefault(r["method"], {})[r["status"]] = \
            per.setdefault(r["method"], {}).get(r["status"], 0) + 1
    bad = []
    for method, t in sorted(per.items()):
        if t.get("DIVERGED") or t.get("r-error"):
            bad.append((t.get("DIVERGED", 0), t.get("r-error", 0), method, t))
    bad.sort(reverse=True)
    print(f"\nmethods needing work: {len(bad)} of {len(per)}")
    for ndiv, nerr, method, t in bad:
        print(f"  {method:34s} diverged={ndiv:2d} r-error={nerr:2d}  " +
              " ".join(f"{k}:{v}" for k, v in sorted(t.items())))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
