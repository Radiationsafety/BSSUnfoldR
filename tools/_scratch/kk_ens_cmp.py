"""Ensemble comparison of randomised-Kaczmarz runs: mean spectra, spread,
integral agreement.  Usage: python kk_ens_cmp.py <fixture_dir_name>"""
import json
import sys

import numpy as np

fixture = sys.argv[1] if len(sys.argv) > 1 else "alt_jinr"
d = json.load(open("/tmp/ab_all_%s.json" % fixture))
PY = json.load(open("/tmp/kk_ens_py_%s.json" % fixture))
MAXIT = max(max(PY[c]["iterations"]) for c in PY)
print("fixture=%s  max_iterations=%d  seeds=%d" % (fixture, MAXIT, len(PY["moderated"]["specs"])))
hdr = ("case        rel(mean)   cos(mean)   sum Py      sum R      int ratio  "
       "sd/|mean| Py  sd/|mean| R   per-seed vs Py-mean")
print(hdr)
for case in ["moderated", "fission", "evap", "split"]:
    Sp = np.squeeze(np.asarray(PY[case]["specs"], float)).reshape(len(PY[case]["specs"]), -1)
    import os
    per = "/tmp/kk_ens_r_%s_%s.json" % (fixture, case)
    comb = "/tmp/kk_ens_r_%s.json" % fixture
    r = json.load(open(per if os.path.exists(per) else comb))[case]
    Sr = np.squeeze(np.asarray(r["specs"], float)).reshape(len(r["specs"]), -1)
    mp, mr = Sp.mean(axis=0), Sr.mean(axis=0)
    rel = np.linalg.norm(mp - mr) / np.linalg.norm(mp)
    cos = mr @ mp / np.linalg.norm(mr) / np.linalg.norm(mp)
    sd_py = np.mean(np.linalg.norm(Sp - mp, axis=1) / np.linalg.norm(mp))
    sd_r = np.mean(np.linalg.norm(Sr - mr, axis=1) / np.linalg.norm(mr))
    cross = np.mean(np.linalg.norm(Sr - mp, axis=1) / np.linalg.norm(mp))
    # python-vs-python reference spread (same stream, different seeds)
    print("%-11s %.3e   %.6f   %.6g   %.6g   %.5f      %.3e      %.3e       %.3e" % (
        case, rel, cos, mp.sum(), mr.sum(), mr.sum() / mp.sum(),
        sd_py, sd_r, cross))
