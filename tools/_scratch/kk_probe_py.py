"""Statistical parity probe for randomized Kaczmarz: seed ensembles."""
import importlib
import json

import numpy as np

K = importlib.import_module("bssunfold.core.unfold_randomized_kaczmarz")
fixture = "alt_jinr"
d = json.load(open("/tmp/ab_all_%s.json" % fixture))
out = {}
for case in ["moderated", "fission", "evap", "split"]:
    A = np.asarray(d[case]["A"], float)
    b = np.asarray(d[case]["b"], float)
    x0 = np.zeros(A.shape[1])
    specs = []
    iters = []
    for seed in range(1, 25):
        x, it, conv = K.solve_randomized_kaczmarz(A, b, x0, random_state=seed)
        specs.append(np.asarray(x, float))
        iters.append(int(it))
    S = np.array(specs)
    mean = S.mean(axis=0)
    rel = np.linalg.norm(S - mean, axis=1) / np.linalg.norm(mean)
    # pairwise
    pw = []
    for i in range(0, 8):
        for j in range(i + 1, 8):
            pw.append(np.linalg.norm(S[i] - S[j]) / np.linalg.norm(mean))
    out[case] = {"mean_sum": float(mean.sum()), "sums": [float(s.sum()) for s in S],
                 "rel_vs_mean": [float(v) for v in rel],
                 "pairwise": [float(v) for v in pw],
                 "iterations": iters}
    print("%-10s mean sum %.6g | sum spread rel %.3e | rel-to-mean median %.3e "
          "min %.3e max %.3e | pairwise median %.3e | iters %s" % (
              case, mean.sum(), np.std([s.sum() for s in S]) / mean.sum(),
              np.median(rel), rel.min(), rel.max(), np.median(pw),
              sorted(set(iters))))
json.dump(out, open("/tmp/kk_py_%s.json" % fixture, "w"))
