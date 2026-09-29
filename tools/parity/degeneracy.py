#!/usr/bin/env python
"""Degeneracy-aware parity classifier.

Spectrum-level rel_l2 is only a meaningful parity gate when the response
matrix has full column rank.  Some fixtures violate that: the JINR spheres in
/tmp/alt_jinr have *identical* sensitivities for the lowest energy bins, so
any allocation of flux among those bins reproduces the readings exactly and
the solver's vertex choice is arbitrary.  A port can therefore be correct and
still show rel_l2 ~ 1.

For every (method, case) this script reports three independent numbers:

  spectrum   rel_l2 between the R and Python spectra
  fit        ||A x - b|| / ||A|| ||b|| computed here, in Python, from the
             spectrum each side actually exported (so it is not self-reported)
  shape      rel_l2 after collapsing the columns of A that are numerically
             identical, i.e. the part of the difference that no measurement
             can constrain

Verdict DIVERGED only when the collapsed-shape difference is large as well;
DEGENERATE means the two spectra are the same solution up to an unresolvable
tie, and NOISE means both are within solver tolerance.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CASES = ("moderated", "fission", "evap", "split")

_spec = importlib.util.spec_from_file_location("one", ROOT / "tools" / "parity" / "one.py")
one = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(one)


def rel(a, b):
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    n = min(len(a), len(b))
    return float(np.linalg.norm(a[:n] - b[:n]) / max(np.linalg.norm(b[:n]), 1e-300))


def column_groups(A, tol=1e-10):
    """Partition the columns of A that are numerically identical."""
    norm = np.linalg.norm(A, axis=0)
    parent = list(range(A.shape[1]))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(A.shape[1]):
        for j in range(i + 1, A.shape[1]):
            if norm[i] == 0 or norm[j] == 0:
                continue
            d = np.linalg.norm(A[:, i] - A[:, j]) / max(norm[i], norm[j])
            if d <= tol:
                parent[find(i)] = find(j)
    groups: dict[int, list[int]] = {}
    for i in range(A.shape[1]):
        groups.setdefault(find(i), []).append(i)
    return [g for g in groups.values() if len(g) > 1]


def collapse(x, groups):
    """Sum each degenerate group into its first member; keep the rest."""
    x = np.asarray(x, float).ravel().copy()
    out = list(x)
    for g in groups:
        total = sum(x[i] for i in g)
        out[g[0]] = total
        for i in g[1:]:
            out[i] = 0.0
    return np.array(out[: len(x)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("methods", nargs="+")
    ap.add_argument("--fixtures", nargs="+",
                    default=["/tmp/bsscmp", "/tmp/alt_gsf",
                             "/tmp/alt_coarse", "/tmp/alt_jinr"])
    ap.add_argument("--tol", type=float, default=1e-4)
    args = ap.parse_args()

    for fixture in args.fixtures:
        os.environ["PARITY_OUT"] = fixture
        one.HERE = Path(fixture)
        det, fx = one.load_det()
        # Amat is (n_bins, n_det) and already carries the ln_steps lethargy
        # weighting; transposed it is the solver's design matrix, whose
        # columns are the energy bins.
        M = np.asarray(det.Amat, float).T
        groups = column_groups(M)
        tag = Path(fixture).name
        if groups:
            print(f"# {tag}: {len(groups)} degenerate column group(s): {groups}")
        for method in args.methods:
            for case in CASES:
                readings = dict(zip(fx["detector_names"],
                                    fx["cases"][case]["readings"]))
                b = np.array([readings[n] for n in fx["detector_names"]], float)
                try:
                    xp = np.asarray(one.py_spectrum(method, case, {}, None)[0], float)
                    xr = np.asarray(one.r_spectrum(method, case, {}, None)[0], float)
                except SystemExit:
                    raise
                except Exception as exc:  # noqa: BLE001
                    print(f"{tag:11s} {method:26s} {case:10s} SKIP "
                          f"{type(exc).__name__}: {str(exc)[:70]}")
                    continue
                n = min(len(xp), len(xr), M.shape[1])
                xp, xr = xp[:n], xr[:n]
                scale = max(np.linalg.norm(M, 2) * np.linalg.norm(b), 1e-300)
                fp = np.linalg.norm(M @ xp - b) / scale
                fr = np.linalg.norm(M @ xr - b) / scale
                sp = rel(collapse(xr, groups), collapse(xp, groups))
                spectrum = rel(xr, xp)
                if max(fp, fr) <= 1e-8 and sp <= args.tol:
                    verdict = "DEGENERATE" if spectrum > args.tol else "OK"
                elif sp <= args.tol:
                    verdict = "NOISE"
                else:
                    verdict = "OK" if sp <= args.tol else "DIVERGED"
                print(f"{tag:11s} {method:26s} {case:10s} "
                      f"spectrum={spectrum:9.2e} fit(py={fp:8.1e},r={fr:8.1e}) "
                      f"collapsed={sp:9.2e}  {verdict}")


if __name__ == "__main__":
    main()
