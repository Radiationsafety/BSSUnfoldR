"""Compare D, cond, D_inv, FI, omega, delta for one fixed (alpha,beta) point."""
import importlib
import json
import os
import sys

import numpy as np

U = importlib.import_module("bssunfold.core.unfold_reconst")
fixture, case = sys.argv[1], sys.argv[2]
alpha, beta = float(sys.argv[3]), float(sys.argv[4])

d = json.load(open("/tmp/ab_all_%s.json" % fixture))
A = np.asarray(d[case]["A"], float)
b = np.asarray(d[case]["b"], float)
M, N = A.shape
F = b.astype(float)
S = np.sqrt(np.maximum(F, 1e-10))
sa = np.exp(np.mean(np.log(np.maximum(S, 1e-300))))
S_norm = S / sa
W = A / S_norm[:, None]
B = W.T @ W
A_vec = A.T @ (F / S_norm ** 2)
OMO = U._build_omo_matrix(N, 1e-3)
Omega = U._omo_to_full(OMO, N)
D = beta * B + alpha * Omega
cond = np.linalg.cond(D)
D_inv = U._invert_system(D)
FI = D_inv @ A_vec * beta
om = U._compute_omega(OMO, D_inv, FI, N, alpha)
de = U._compute_delta(B, D_inv, FI, A_vec, F, S_norm, N, M, beta)
out = {"cond": float(cond), "D": D.tolist(), "D_inv": D_inv.tolist(),
       "FI": FI.tolist(), "omega": float(om), "delta": float(de),
       "d1": float(np.trace(B @ D_inv)), "d2": float(FI @ B @ FI),
       "d3": float(A_vec @ FI), "d4": float(np.sum((F / S_norm) ** 2)),
       "mOverBeta": float(M) / beta, "nOverAlpha": float(N) / alpha,
       "Bdiag": np.diag(B).tolist(), "Avec": A_vec.tolist(),
       "B": B.tolist(), "Omega": Omega.tolist(),
       "OMO": OMO.tolist(), "sa": float(sa),
       "Snorm": S_norm.tolist()}
json.dump(out, open("/tmp/pt_py_%s_%s.json" % (fixture, case), "w"))
print("cond %.6e omega %.10g delta %.10g d1 %.10g d2 %.10g d3 %.10g d4 %.10g" % (
    cond, om, de, out["d1"], out["d2"], out["d3"], out["d4"]))
