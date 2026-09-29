"""Is the R reconst port inside Python's own 1-ulp noise cone?

For every parity cell: build 12 replicas of (A, b) where each entry is
perturbed by at most one ulp (relative 1e-16), unfold each with Python's
solve_reconst, and compare the resulting spread with the distance of the R
port's answer (computed from the same unperturbed A, b) to the reference.
"""
import json

import numpy as np
import importlib

U = importlib.import_module("bssunfold.core.unfold_reconst")
EPS = 1e-16
NREP = 12
rng = np.random.default_rng(20260930)
print("cell         py_spread(max)  py_spread(med)  R_dev   R_dev/max_spread")
for fixture in ["bsscmp", "alt_gsf", "alt_coarse", "alt_jinr"]:
    d = json.load(open("/tmp/ab_all_%s.json" % fixture))
    for case in ["moderated", "fission", "evap", "split"]:
        A = np.asarray(d[case]["A"], float)
        b = np.asarray(d[case]["b"], float)
        s0 = np.asarray(U.solve_reconst(A, b), float)
        reps = []
        for _ in range(NREP):
            A2 = A * (1.0 + EPS * rng.uniform(-1, 1, size=A.shape))
            b2 = b * (1.0 + EPS * rng.uniform(-1, 1, size=b.shape))
            reps.append(np.asarray(U.solve_reconst(A2, b2), float))
        R = json.load(open("/tmp/rec_r_%s_%s.json" % (fixture, case)))["spectrum"]
        R = np.asarray(R, float)
        nrm = np.linalg.norm(s0)
        dev = [np.linalg.norm(r - s0) / nrm for r in reps]
        cross = [np.linalg.norm(r - R) / np.linalg.norm(R) for r in reps]
        rdev = np.linalg.norm(R - s0) / nrm
        print("%-7s %-9s %.3e      %.3e         %.3e   %s" % (
            fixture, case, max(dev), np.median(dev), rdev,
            "INSIDE" if rdev <= max(dev) else "outside"))
        json.dump({"py_spread": dev, "py_cross_R": cross, "r_dev": rdev},
                  open("/tmp/noisecone_%s_%s.json" % (fixture, case), "w"))
