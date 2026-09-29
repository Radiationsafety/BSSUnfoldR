#!/usr/bin/env python
"""Null-space diagnostic for degenerate cells.

Two spectra that differ only in the component of the solution lying in the
null space of the response matrix fit the measurements equally well, so a
large spectrum-level rel_l2 does not by itself mean the R port is wrong.
This script reports, for each (method, case), the spectral rel_l2 together
with the physically meaningful invariants: the detector residual norms and
the derived dose rates.  When the invariants agree but the spectra do not,
the cell is a null-space tie rather than a porting defect.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("one", ROOT / "tools" / "parity" / "one.py")
one = importlib.util.module_from_spec(spec)
spec.loader.exec_module(one)

PAIRS = [
    ("unfold_admm", "evap"),
    ("unfold_admm", "split"),
    ("unfold_ferdor", "evap"),
    ("unfold_ferdor", "split"),
    ("unfold_frank_wolfe", "moderated"),
    ("unfold_frank_wolfe", "evap"),
    ("unfold_hybrid_gmres", "split"),
    ("unfold_reconst", "fission"),
    ("unfold_reconst", "split"),
]


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


def main():
    fixture = Path(sys.argv[1]) if len(sys.argv) > 1 else "/tmp/alt_jinr"
    one.HERE = fixture
    for method, case in PAIRS:
        try:
            xp, pres = one.py_spectrum(method, case, {}, None)
            xr, rres = one.r_spectrum(method, case, {}, None)
        except Exception as exc:  # noqa: BLE001
            print(f"{method:24s} {case:10s} SKIP {type(exc).__name__}: {exc}")
            continue
        rp, rr = pres.get("residual_norm"), rres.get("residual_norm")
        dkey = None
        try:
            pyd = np.asarray(list(pres["doserates"].values()), float)
            rd_ = np.asarray(list(rres["doserates"].values()), float)
            dkey = rel(rd_, pyd)
        except Exception:  # noqa: BLE001
            pass
        print(
            f"{method:24s} {case:10s} "
            f"spectrum={rel(xr, xp):9.3e} "
            f"residual(py={rp if rp is None else format(rp, '.3e')}, "
            f"r={rr if rr is None else format(rr, '.3e')}) "
            f"doserate={dkey if dkey is None else format(dkey, '.3e')}"
        )


if __name__ == "__main__":
    main()
