"""Record np.linalg.cond(D) for every system matrix built by reconst."""
import importlib
import json
import sys

import numpy as np

U = importlib.import_module("bssunfold.core.unfold_reconst")
COND = []
orig = U._build_system_matrix


def wrapped(B, OMO, n, alpha, beta):
    D = orig(B, OMO, n, alpha, beta)
    COND.append([float(alpha), float(beta),
                 float(np.linalg.cond(D)), int(D.shape[0])])
    return D


U._build_system_matrix = wrapped
for fixture, case in [("alt_gsf", "moderated"), ("alt_gsf", "fission"),
                      ("bsscmp", "split"), ("alt_jinr", "moderated"),
                      ("alt_jinr", "fission"), ("alt_coarse", "moderated")]:
    d = json.load(open("/tmp/ab_all_%s.json" % fixture))
    A = np.asarray(d[case]["A"], float)
    b = np.asarray(d[case]["b"], float)
    COND.clear()
    s = np.asarray(U.solve_reconst(A, b), float)
    C = np.array(COND)
    near = np.abs(np.log10(C[:, 2]) - 12.0)
    print("%-8s %-9s n_inv=%3d cond range [%.3e, %.3e]  #within 0.3dec of 1e12: %d  min|log10-12|=%.3f"
          % (fixture, case, len(C), C[:, 2].min(), C[:, 2].max(),
             int((near < 0.3).sum()), near.min()))
    json.dump(COND, open("/tmp/cond_py_%s_%s.json" % (fixture, case), "w"))
