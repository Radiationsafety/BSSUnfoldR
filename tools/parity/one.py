#!/usr/bin/env python
"""Per-method parity check: compare one R method against bssunfold.

Usage:
    python tools/parity/one.py <method> [case] [kwarg=value ...]

Runs the Python method in-process and the R method via tools/parity/one.R,
then prints rel_l2 / cosine / integral ratio and the worst bins. Exit status
is 0 when rel_l2 <= --tol (default 1e-6), 1 otherwise, 2 on a hard error.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(os.environ.get("PARITY_OUT", "/tmp/bsscmp"))
ROOT = Path(__file__).resolve().parents[2]


def load_det():
    from bssunfold import Detector

    fx = json.loads((HERE / "fixture.json").read_text())
    rf = {"E_MeV": fx["E_MeV"]}
    for name, sens in zip(fx["detector_names"], fx["sensitivities"]):
        rf[name] = sens
    return Detector(response_functions=rf), fx


def py_spectrum(method, case, kwargs, seed):
    det, fx = load_det()
    readings = dict(zip(fx["detector_names"], fx["cases"][case]["readings"]))
    fn = getattr(det, method, None)
    if fn is None:
        raise RuntimeError(f"{method} not present in bssunfold Detector")
    import inspect

    if "random_state" in inspect.signature(fn).parameters:
        kwargs.setdefault("random_state", seed)
    res = fn(readings, **kwargs)
    return np.asarray(res["spectrum"], dtype=float), res


def r_spectrum(method, case, kwargs, seed):
    cmd = ["Rscript", str(ROOT / "tools" / "parity" / "one.R"), method, case]
    for k, v in kwargs.items():
        cmd.append(f"{k}={v}")
    if seed is not None:
        cmd.append(f"random_state={seed}")
    env = dict(os.environ, PARITY_OUT=str(HERE))
    p = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=str(ROOT))
    if p.returncode != 0:
        raise RuntimeError(p.stderr.strip().splitlines()[-1] if p.stderr.strip()
                           else "Rscript failed")
    payload = json.loads(p.stdout.strip().splitlines()[-1])
    if payload.get("error"):
        raise RuntimeError(payload["error"])
    return np.asarray(payload["spectrum"], dtype=float), payload


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("method")
    ap.add_argument("case", nargs="?", default="moderated")
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--seed", type=int,
                    default=int(os.environ.get("PARITY_SEED", "42")))
    ap.add_argument("--no-rand", action="store_true",
                    help="do not pass random_state (stochastic-free compare)")
    ap.add_argument("kw", nargs="*", help="extra kwargs as key=value")
    args = ap.parse_args()

    kwargs = {}
    for item in args.kw:
        k, _, v = item.partition("=")
        try:
            kwargs[k] = int(v) if v.lstrip("-").isdigit() else float(v)
        except ValueError:
            kwargs[k] = {"TRUE": True, "FALSE": False}.get(v, v)

    seed = None if args.no_rand else args.seed
    try:
        xp, pres = py_spectrum(args.method, args.case, dict(kwargs), seed)
    except Exception as exc:  # noqa: BLE001
        print(f"PY ERROR: {type(exc).__name__}: {exc}")
        sys.exit(2)
    try:
        xr, rres = r_spectrum(args.method, args.case, dict(kwargs), seed)
    except Exception as exc:  # noqa: BLE001
        print(f"R  ERROR: {exc}")
        sys.exit(2)

    n = min(len(xr), len(xp))
    xr, xp = xr[:n], xp[:n]
    rel = float(np.linalg.norm(xr - xp) / max(np.linalg.norm(xp), 1e-300))
    cos = float(xr @ xp / (max(np.linalg.norm(xr), 1e-300) *
                           max(np.linalg.norm(xp), 1e-300)))
    ratio = float(xr.sum() / xp.sum()) if xp.sum() else float("nan")
    print(f"{args.method} [{args.case}] kwargs={kwargs or '{}'}")
    print(f"  rel_l2   = {rel:.3e}")
    print(f"  cos      = {cos:.8f}")
    print(f"  int R/Py = {ratio:.6f}")
    print(f"  sum Py   = {xp.sum():.6f}   sum R = {xr.sum():.6f}")
    sig = xp > 1e-8 * max(xp.max(), 1e-300)
    if sig.any():
        worst = np.argsort(-np.abs(xr[sig] - xp[sig]) / xp[sig])[:5]
        e = np.abs(xr[sig] - xp[sig]) / xp[sig]
        print("  worst pointwise rel dev:",
              ", ".join(f"{v:.3g}" for v in e[worst][:5]))
    print(f"  verdict  = {'OK' if rel <= args.tol else 'DIVERGED'} (tol {args.tol:g})")
    sys.exit(0 if rel <= args.tol else 1)


if __name__ == "__main__":
    main()
