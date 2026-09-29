"""How sensitive is Python's own reconst result to a 1-ulp change in A?"""
import importlib
import json
import os

import numpy as np

U = importlib.import_module("bssunfold.core.unfold_reconst")

for fixture in ["bsscmp", "alt_gsf", "alt_coarse", "alt_jinr"]:
    d = json.load(open("/tmp/ab_all_%s.json" % fixture))
    for case in ["moderated", "fission", "evap", "split"]:
        A = np.asarray(d[case]["A"], float)
        b = np.asarray(d[case]["b"], float)
        s0 = np.asarray(U.solve_reconst(A, b), float)
        A2 = A.copy()
        A2[0, 0] = np.nextafter(A2[0, 0], 2 * A2[0, 0])
        s1 = np.asarray(U.solve_reconst(A2, b), float)
        A3 = A.copy()
        A3 *= (1.0 + 1e-12)
        s2 = np.asarray(U.solve_reconst(A3, b), float)
        r1 = np.linalg.norm(s1 - s0) / np.linalg.norm(s0)
        r2 = np.linalg.norm(s2 - s0) / np.linalg.norm(s0)
        print("%-10s %-10s rel(1ulp)=%.3e rel(1e-12)=%.3e" %
              (fixture, case, r1, r2))
