import importlib
import json

import numpy as np

U = importlib.import_module("bssunfold.core.unfold_reconst")
fixture, case = "alt_gsf", "moderated"
d = json.load(open("/tmp/ab_all_%s.json" % fixture))
A = np.asarray(d[case]["A"], float)
b = np.asarray(d[case]["b"], float)
M, N = A.shape
S = np.sqrt(np.maximum(b, 1e-10))
sa = float(np.exp(np.mean(np.log(np.maximum(S, 1e-300)))))
Sn = S / sa
W = A / Sn[:, None]
B = W.T @ W
Avec = A.T @ (b / Sn**2)
OMO = U._build_omo_matrix(N, 1e-3)
Om = U._omo_to_full(OMO, N)
alpha = 4.3543330946655486e-09
beta = 7.616246847424449
e = np.asarray(d[case]["spectrum"], float)


def fi_of(al, be, inv=np.linalg.inv):
    D = be * B + al * Om
    X = inv(D)
    return np.maximum(X @ Avec * be, 0), D, X


FI, D, X = fi_of(alpha, beta)
print("exact-internal-beta FI rel vs auto spectrum: %.3e" %
      (np.linalg.norm(FI - e) / np.linalg.norm(e)))
json.dump({"D_inv": X.tolist(), "D": D.tolist(), "B": B.tolist(),
           "Om": Om.tolist(), "Avec": Avec.tolist(), "sa": sa,
           "alpha": alpha, "beta": beta, "FI": FI.tolist(),
           "spectrum": e.tolist()},
          open("/tmp/final_alt_gsf_moderated.json", "w"))
for lab, f in [("beta+1e-7", (1 + 1e-7, 1.0)), ("beta+1e-5", (1 + 1e-5, 1.0)),
               ("alpha+1e-7", (1.0, 1 + 1e-7)), ("alpha+1e-5", (1.0, 1 + 1e-5))]:
    g, _, _ = fi_of(alpha * f[1], beta * f[0])
    print("  %-10s FI rel vs base %.3e" % (lab, np.linalg.norm(g - FI) / np.linalg.norm(FI)))
# solve-based (no explicit inverse)
g2 = np.maximum(np.linalg.solve(D, Avec) * beta, 0)
print("  solve vs inv FI rel %.3e" % (np.linalg.norm(g2 - FI) / np.linalg.norm(FI)))
print("  cond %.3e" % np.linalg.cond(D))
print("  FI[:6]", np.round(FI[:6], 6), " spec[:6]", np.round(e[:6], 6))
