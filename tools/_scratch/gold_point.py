import json

import numpy as np

p = json.load(open("/tmp/pt_py_alt_gsf_moderated.json"))
r = json.load(open("/tmp/pt_r_alt_gsf_moderated.json"))
D = np.array(p["D"], float)
Avec = np.array(p["Avec"], float)
beta = 0.017742987176815786
rhs = Avec * beta
x_py = np.array(p["FI"], float)
x_r = np.array(r["FI"], float)

# iterative refinement in float64 with longdouble residuals
Xl = np.zeros(len(rhs), dtype=np.longdouble)
Dl = D.astype(np.longdouble)
rhsl = rhs.astype(np.longdouble)
x = np.linalg.solve(D, rhs)
for it in range(6):
    res = rhsl - Dl @ x.astype(np.longdouble)
    dx = np.linalg.solve(D, res.astype(float))
    x = x + dx
    Xl = x.astype(np.longdouble)
gold = x
print("py vs gold rel  %.3e" % (np.linalg.norm(x_py - gold) / np.linalg.norm(gold)))
print("r  vs gold rel  %.3e" % (np.linalg.norm(x_r - gold) / np.linalg.norm(gold)))
print("py vs r  rel    %.3e" % (np.linalg.norm(x_py - x_r) / np.linalg.norm(gold)))
print("cond %.3e" % np.linalg.cond(D))
# delta with gold FI and py/r D_inv
B = np.array(p["B"], float)
d4 = p["d4"]
mb = p["mOverBeta"]


def delta_of(Xinv, FI):
    d1 = np.trace(B @ Xinv)
    d2 = FI @ B @ FI
    d3 = Avec @ FI
    return mb - (d1 + d2 - 2 * d3 + d4), d1, d2, d3


# gold D_inv: solve D X = I with refinement
I = np.eye(len(D))
Xg = np.linalg.solve(D, I)
for it in range(6):
    res = (I - D @ Xg).astype(np.longdouble)
    Xg = Xg + np.linalg.solve(D, res.astype(float))
FIg = Xg @ rhs
for nm, X, FI in [("py", np.array(p["D_inv"], float), x_py),
                  ("r", np.array(r["D_inv"], float), x_r),
                  ("gold", Xg, FIg)]:
    de, d1, d2, d3 = delta_of(X, FI)
    print("%-5s delta=%.6e d1=%.10f d2=%.10g d3=%.10g" % (nm, de, d1, d2, d3))
print("py delta stored %.6e ; r delta stored %.6e" % (p["delta"], r["delta"]))
