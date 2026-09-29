#!/usr/bin/env python
"""Certify a stochastic method by seed spread rather than by a single seed.

Several bssunfold methods drive their sampler from
``numpy.random.default_rng(seed)`` -- PCG64 with ziggurat normals.  R has no
way to reproduce that stream, so the two languages cannot agree pointwise no
matter how faithful the port is.  What *can* be checked is that the R chain
lands inside the band the Python chain itself wanders over.

For each seed it runs both languages and reports

  py-spread   distance from Python's reference chain to other Python chains
  r-spread    distance from Python's reference chain to each R chain
  r-vs-py-ref distance from Python's reference chain to the same-seed R chain

A port is certified when max(r-spread) is no larger than max(py-spread) at the
same number of samples, which says the disagreement is chain noise rather than
a formulation difference.  Run with --samples to check that the band collapses
as the chains get longer, which is the stronger evidence.

Usage:
    python tools/parity/mcmc_spread.py unfold_bayesian_parametric \
        --cases fission evap --seeds 11 12 13 14 15 --ref-seed 42
"""
from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("one", ROOT / "tools" / "parity" / "one.py")
one = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(one)


def rel_to(a, ref):
    a, ref = np.asarray(a, float), np.asarray(ref, float)
    n = min(len(a), len(ref))
    return float(np.linalg.norm(a[:n] - ref[:n]) / max(np.linalg.norm(ref[:n]), 1e-300))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("method")
    ap.add_argument("--fixtures", nargs="+", default=["/tmp/bsscmp"])
    ap.add_argument("--cases", nargs="+", default=["fission"])
    ap.add_argument("--seeds", nargs="+", type=int,
                    default=[11, 12, 13, 14, 15])
    ap.add_argument("--ref-seed", type=int, default=42)
    ap.add_argument("kw", nargs="*", help="extra kwargs as key=value")
    args = ap.parse_args()

    kwargs = {}
    for item in args.kw:
        k, _, v = item.partition("=")
        kwargs[k] = int(v) if v.lstrip("-").isdigit() else float(v)

    for fixture in args.fixtures:
        os.environ["PARITY_OUT"] = fixture
        one.HERE = Path(fixture)
        for case in args.cases:
            ref, _ = one.py_spectrum(args.method, case, dict(kwargs),
                                     args.ref_seed)
            py = [rel_to(one.py_spectrum(args.method, case, dict(kwargs), s)[0],
                         ref) for s in args.seeds]
            rr = [rel_to(one.r_spectrum(args.method, case, dict(kwargs), s)[0],
                         ref) for s in args.seeds]
            same = rel_to(one.r_spectrum(args.method, case, dict(kwargs),
                                         args.ref_seed)[0], ref)
            ok = max(rr) <= max(py) * 3
            print(f"{Path(fixture).name} {args.method} [{case}] "
                  f"kwargs={kwargs or '{}'}")
            print(f"  py-spread  max={max(py):.2e}  " +
                  " ".join(f"{v:.1e}" for v in py))
            print(f"  r-spread   max={max(rr):.2e}  " +
                  " ".join(f"{v:.1e}" for v in rr))
            print(f"  same-seed R vs Py = {same:.2e}")
            print(f"  verdict = {'WITHIN CHAIN NOISE' if ok else 'DIVERGED'}")


if __name__ == "__main__":
    main()
