"""Dump the full _reg1/omega/delta trajectory of Python's reconst for one cell.

Usage: python dump_traj_py.py <fixture> <case>
Writes /tmp/traj_py_<fixture>_<case>.json
Also, if the env var FIX_AB is set, reads A/b from /tmp/ab_all_<fixture>.json
for the given case and calls solve_reconst directly (no Detector involvement).
"""
import importlib
import json
import os

import numpy as np

fixture = os.environ.get("TRAJ_FX", "alt_gsf")
case = os.environ["TRAJ_CASE"]
U = importlib.import_module("bssunfold.core.unfold_reconst")

LOG = {"reg1": [], "omega": [], "delta": []}
o_reg1, o_om, o_de = U._reg1, U._compute_omega, U._compute_delta


def wreg1(B, OMO, A_vec, n, alpha, beta, ich):
    LOG["reg1"].append([float(alpha), float(beta), int(ich)])
    return o_reg1(B, OMO, A_vec, n, alpha, beta, ich)


def wom(OMO, D_inv, FI, n, alpha):
    v = o_om(OMO, D_inv, FI, n, alpha)
    LOG["omega"].append(float(v))
    return v


def wde(B, D_inv, FI, A_vec, F, S, n, m, beta):
    v = o_de(B, D_inv, FI, A_vec, F, S, n, m, beta)
    LOG["delta"].append(float(v))
    return v


U._reg1, U._compute_omega, U._compute_delta = wreg1, wom, wde

if os.environ.get("FIX_AB"):
    d = json.load(open("/tmp/ab_all_%s.json" % fixture))
    A = np.asarray(d[case]["A"], float)
    b = np.asarray(d[case]["b"], float)
    spec = np.asarray(U.solve_reconst(A, b), float)
else:
    fx = json.load(open(os.path.join("/tmp", fixture, "fixture.json")))
    rf = {"E_MeV": fx["E_MeV"]}
    for name, sens in zip(fx["detector_names"], fx["sensitivities"]):
        rf[name] = sens
    from bssunfold import Detector
    det = Detector(response_functions=rf)
    readings = dict(zip(fx["detector_names"], fx["cases"][case]["readings"]))
    spec = np.asarray(det.unfold_reconst(readings)["spectrum"], float)

json.dump({"reg1": LOG["reg1"], "omega": LOG["omega"], "delta": LOG["delta"],
           "spectrum": spec.tolist()},
          open("/tmp/traj_py_%s_%s.json" % (fixture, case), "w"))
print(fixture, case, "n_reg1", len(LOG["reg1"]), "omega", len(LOG["omega"]),
      "delta", len(LOG["delta"]))
