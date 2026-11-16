#!/bin/sh
# Level A reproduction (docs/REPRODUCE.md): no quantum chemistry, a few minutes.
# Needs Python 3.12 with environment/paper-requirements.txt (set PYTHON=... if it is not python3), xz, sha256sum,
# and the CIF of CCDC 835950 as structure/BTPC.cif (or BTPC_CIF=/path/to/BTPC.cif). pandoc is optional (PANDOC=...
# or on the PATH); without it the pandoc part of make_tables.py --check is skipped.
set -eu
cd "$(dirname "$0")/.."
PY=${PYTHON:-python3}
BTPC_CIF=${BTPC_CIF:-structure/BTPC.cif}
export BTPC_CIF

echo "== 1. checksums of the deposited files"
sha256sum -c --quiet SHA256SUMS
echo "all files OK"

echo "== 2. unpack the grids"
sh scripts/unpack_grids.sh > /dev/null
echo "22 grids unpacked and checked"

echo "== 3. the CIF and the ion pair"
echo "17065c7ac118c838b386de8448cb7d631073c64e02cd2cd4d7e929e34c5349b5  $BTPC_CIF" | sha256sum -c
"$PY" scripts/dft/build_cluster.py "$BTPC_CIF" btpc --title Benzotriazolium-picrate
echo "dfa71f6f1a4b15cf9b9f8f4cf4240f88ef284b073fef5971cfe9c3e5a0f95476  structure/btpc_cluster.xyz" | sha256sum -c
"$PY" scripts/dft/curate_cases.py > /dev/null
cmp structure/btpc_gas.json cases/btpc_gas.json
echo "structure/btpc_gas.json: identical to cases/btpc_gas.json"

echo "== 4. the facts register"
REF=$(mktemp)
cp analysis/facts.json "$REF"
"$PY" scripts/analysis/build_facts.py
"$PY" scripts/tools/compare_facts.py "$REF" analysis/facts.json
rm -f "$REF"
"$PY" scripts/analysis/build_facts_md.py

echo "== 5. the tables"
if [ -n "${PANDOC:-}" ] || command -v pandoc > /dev/null 2>&1; then
  "$PY" scripts/tables/make_tables.py --check
else
  echo "pandoc not found: tables built without the pandoc checks"
  "$PY" scripts/tables/make_tables.py
fi
"$PY" scripts/tables/tables_to_csv.py > /dev/null

echo "== 6. every rebuilt file must be byte-identical to the deposit"
sha256sum -c --quiet SHA256SUMS
echo "SHA256SUMS: all files OK"
