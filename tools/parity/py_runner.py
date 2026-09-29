#!/usr/bin/env python
"""Run every bssunfold unfolding method on the shared parity fixture.

Writes {method: {spectrum, residual_norm, iterations, converged}} per case to
results-python.json. Methods that raise are recorded with an `error` string so
the comparison table shows coverage rather than silently dropping them.
"""
from __future__ import annotations

import inspect
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

HERE = Path(os.environ.get("PARITY_OUT", "/tmp/bsscmp"))
OUT_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "py"


def _jsonable(v):
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, np.ndarray):
        return [_jsonable(x) for x in v.ravel().tolist()]
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, (bool, int, float, str)) or v is None:
        return v
    return str(v)


def build_detector(fx):
    from bssunfold import Detector

    rf = {"E_MeV": fx["E_MeV"]}
    for name, sens in zip(fx["detector_names"], fx["sensitivities"]):
        rf[name] = sens
    return Detector(response_functions=rf)


def dump(det, fx, case_name, methods, seed):
    readings = dict(zip(fx["detector_names"], fx["cases"][case_name]["readings"]))
    out = {}
    for mname in methods:
        fn = getattr(det, mname, None)
        if fn is None:
            out[mname] = {"error": "missing"}
            continue
        kwargs = {}
        try:
            params = inspect.signature(fn).parameters
            if "random_state" in params:
                kwargs["random_state"] = seed
            if "verbose" in params:
                kwargs["verbose"] = False
        except (TypeError, ValueError):
            pass
        rec = {}
        try:
            res = fn(readings, **kwargs)
            spec = np.asarray(res.get("spectrum"), dtype=float)
            rec["spectrum"] = [float(x) for x in spec]
            for key in ("residual_norm", "iterations", "converged", "chi2",
                        "reduced_chi2", "n_iterations"):
                if key in res and isinstance(res[key], (int, float, bool, np.number)):
                    rec[key] = _jsonable(res[key])
            if "doserates" in res:
                rec["doserates"] = _jsonable(res["doserates"])
        except Exception as exc:  # noqa: BLE001 - harness records all failures
            rec["error"] = f"{type(exc).__name__}: {str(exc)[:220]}"
        out[mname] = rec
    return out


def main():
    fx = json.loads((HERE / "fixture.json").read_text())
    det = build_detector(fx)
    methods = sorted(m for m in dir(det) if m.startswith("unfold_"))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    seed = int(os.environ.get("PARITY_SEED", "42"))
    for case in fx["cases"]:
        res = dump(det, fx, case, methods, seed)
        (OUT_DIR / f"{case}.json").write_text(json.dumps(res))
        ok = sum(1 for v in res.values() if "spectrum" in v)
        print(f"{case}: {ok}/{len(methods)} methods produced a spectrum",
              file=sys.stderr)


if __name__ == "__main__":
    main()
