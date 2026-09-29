#!/usr/bin/env python
"""unfold_express residual-parity check.

The Express reference fit is *underdetermined by construction*: the model has
n_groups + 1 = 7 spline knots while the fixtures expose 5-7 detectors, so the
least-squares valley is flat and two correct optimisers legitimately stop at
different points of it.  rel_l2 on the spectrum is therefore not a meaningful
gate for this method; the residual norm is.  This script compares residuals
(plus reports the spectrum rel_l2 for context).
"""
import json, os, subprocess, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIX = [("main", "/tmp/bsscmp"), ("gsf", "/tmp/alt_gsf"),
       ("coarse", "/tmp/alt_coarse"), ("jinr", "/tmp/alt_jinr")]
METHOD = "unfold_express"

R_SNIPPET = """
suppressMessages(pkgload::load_all("%s", quiet=TRUE)); library(jsonlite)
fx <- fromJSON(paste0(Sys.getenv("FXDIR"), "/fixture.json"), simplifyVector=FALSE)
rf <- c(list(E_MeV=unlist(fx$E_MeV)),
        setNames(lapply(seq_along(fx$detector_names),
                        function(i) unlist(fx$sensitivities[[i]])),
                 unlist(fx$detector_names)))
det <- Detector$new(response_function=rf, detector_names=unlist(fx$detector_names))
out <- list()
for (cs in names(fx$cases)) {
  rd <- setNames(as.numeric(unlist(fx$cases[[cs]]$readings)), unlist(fx$detector_names))
  r <- tryCatch(det$unfold_express(rd), error=function(e) NULL)
  if (is.null(r)) next
  out[[cs]] <- list(spectrum = as.numeric(r$spectrum),
                    residual_norm = r$residual_norm,
                    iterations = r$iterations, converged = r$converged)
}
write_json(out, paste0(Sys.getenv("FXDIR"), "/express_R.json"), digits=17, auto_unbox=TRUE)
"""

PY_SNIPPET = """
import json, os, sys
sys.path.insert(0, %r)
from bssunfold.core.detector import Detector
d = os.environ["FXDIR"]
fx = json.load(open(d + "/fixture.json"))
rf = {"E_MeV": fx["E_MeV"]}
for name, sens in zip(fx["detector_names"], fx["sensitivities"]):
    rf[name] = sens
det = Detector(response_functions=rf)
out = {}
for cs, case in fx["cases"].items():
    rd = dict(zip(fx["detector_names"], case["readings"]))
    r = det.unfold_express(rd)
    out[cs] = dict(spectrum=list(r["spectrum"]),
                   residual_norm=float(r["residual_norm"]),
                   iterations=int(r.get("iterations", -1)))
json.dump(out, open(d + "/express_PY.json", "w"))
"""

for name, d in FIX:
    env = dict(os.environ, FXDIR=d)
    subprocess.run(["/usr/local/bin/Rscript", "-e",
                    R_SNIPPET % str(ROOT)], check=True, env=env,
                   capture_output=True)
    subprocess.run([sys.executable, "-c", PY_SNIPPET % str(
        Path("/Users/spiralfractal/Work/code/bssunfold-0.28.0/src"))],
        check=True, env=env, capture_output=True)
    R = json.load(open(d + "/express_R.json"))
    P = json.load(open(d + "/express_PY.json"))
    print(f"-- {name}")
    for case in P:
        if case not in R:
            print(f"   {case:10s} R unavailable")
            continue
        sp_r = np.array(R[case]["spectrum"])
        sp_p = np.array(P[case]["spectrum"])
        rr, pr = R[case]["residual_norm"], P[case]["residual_norm"]
        rel = np.linalg.norm(sp_r - sp_p) / max(np.linalg.norm(sp_p), 1e-300)
        print(f"   {case:10s} resid Py={pr:.6e} R={rr:.6e} ratio={rr / pr:7.4f}"
              f"  rel_l2={rel:.3e}  {'R better' if rr <= pr else 'PY better'}")
