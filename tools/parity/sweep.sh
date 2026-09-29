#!/usr/bin/env bash
# Run both languages over every parity fixture, then aggregate.
set -uo pipefail
cd "$(dirname "$0")/../.."
PY=${PY:-/tmp/bssvenv/bin/python}
FIXTURES=${FIXTURES:-"/tmp/bsscmp /tmp/alt_gsf /tmp/alt_coarse /tmp/alt_jinr"}

for F in $FIXTURES; do
  if [ ! -f "$F/fixture.json" ]; then
    echo "SKIP $F (no fixture.json)" >&2
    continue
  fi
  echo "=== R   $F ===" >&2
  PARITY_OUT=$F Rscript tools/parity/r_runner.R "$F/r" 2>&1 | tail -6 >&2
  echo "=== PY  $F ===" >&2
  PARITY_OUT=$F "$PY" tools/parity/py_runner.py "$F/py" 2>&1 | tail -6 >&2
done

echo "=== aggregate ===" >&2
"$PY" tools/parity/sweep_compare.py --fixtures $FIXTURES \
     --out "${OUT:-/tmp/parity_matrix.csv}"
