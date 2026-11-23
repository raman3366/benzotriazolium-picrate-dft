#!/bin/sh
# Level A reproduction (docs/REPRODUCE.md): no quantum chemistry, a few minutes.
# Needs Python 3.12 with environment/paper-requirements.txt (set PYTHON=... if it is not python3), xz, sha256sum,
# and the CIF of CCDC 835950 as structure/BTPC.cif (or BTPC_CIF=/path/to/BTPC.cif; keep the file name BTPC.cif).
# pandoc 3.11 or newer is optional (PANDOC=... or on the PATH); with an older pandoc, or none, the tables are built
# without the pandoc part of make_tables.py --check.
set -eu
cd "$(dirname "$0")/.."
PY=${PYTHON:-python3}
BTPC_CIF=${BTPC_CIF:-structure/BTPC.cif}
export BTPC_CIF
PANDOC_MIN=3.11

version_ge() {  # true if dotted version $1 >= $2
  awk -v a="$1" -v b="$2" 'BEGIN { n = split(a, x, "."); m = split(b, y, "."); k = (n > m) ? n : m
    for (i = 1; i <= k; i++) { if (x[i] + 0 > y[i] + 0) exit 0; if (x[i] + 0 < y[i] + 0) exit 1 }
    exit 0 }'
}

echo "== 1. checksums of the deposited files"
sha256sum -c --quiet SHA256SUMS
echo "all files OK"

echo "== 2. unpack the grids"
sh scripts/unpack_grids.sh > /dev/null
echo "44 grids unpacked and checked"

echo "== 3. the CIF and the ion pair"
# The CCDC writes the download date into every CIF it serves (_audit_update_record: "... downloaded from the CCDC.")
# and the copy used here has CRLF line ends, so a fresh download has a different SHA-256. The CIF comparison is
# therefore only reported; the checks that count are the rebuilt ion pair and the case file (structure/README.md).
if [ ! -f "$BTPC_CIF" ]; then
  echo "CIF not found: $BTPC_CIF (download CCDC 835950; structure/README.md)" >&2
  exit 1
fi
CIF_SHA=17065c7ac118c838b386de8448cb7d631073c64e02cd2cd4d7e929e34c5349b5
CIF_NORM_SHA=9d93fe1c427457780670219100894a7ff2b48af9050038d552bc9e82104c67d9
if [ "$(sha256sum < "$BTPC_CIF" | cut -c1-64)" = "$CIF_SHA" ]; then
  echo "$BTPC_CIF: identical to the copy the calculations used"
elif [ "$(tr -d '\r' < "$BTPC_CIF" | grep -v 'downloaded from the CCDC' | sha256sum | cut -c1-64)" = "$CIF_NORM_SHA" ]; then
  echo "$BTPC_CIF: same content as the copy the calculations used (download date or line ends differ)"
else
  echo "WARNING: $BTPC_CIF differs from the copy the calculations used beyond the download date and line ends;"
  echo "         continuing: the rebuilt ion pair below decides (structure/README.md)"
fi
"$PY" scripts/dft/build_cluster.py "$BTPC_CIF" btpc --title Benzotriazolium-picrate
if ! echo "dfa71f6f1a4b15cf9b9f8f4cf4240f88ef284b073fef5971cfe9c3e5a0f95476  structure/btpc_cluster.xyz" | sha256sum -c; then
  echo "the rebuilt ion pair differs from the one the calculations used: check the CIF (structure/README.md)" >&2
  exit 1
fi
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
# make_tables.py --check looks for the table id pandoc writes from the caption ({#tbl:<stem>}); older pandoc
# versions (3.1.3, for example) write <table> without it, so the check runs only with pandoc >= PANDOC_MIN.
PANDOC=${PANDOC:-$(command -v pandoc || true)}
PANDOC_VERSION=
if [ -n "$PANDOC" ]; then
  PANDOC_VERSION=$("$PANDOC" --version 2> /dev/null | sed -n '1s/^pandoc[^0-9]*\([0-9][0-9.]*\).*/\1/p' || true)
fi
if [ -n "$PANDOC_VERSION" ] && version_ge "$PANDOC_VERSION" "$PANDOC_MIN"; then
  export PANDOC
  echo "pandoc $PANDOC_VERSION ($PANDOC)"
  "$PY" scripts/tables/make_tables.py --check
else
  if [ -z "$PANDOC" ]; then
    why="pandoc not found"
  elif [ -z "$PANDOC_VERSION" ]; then
    why="cannot read the version of $PANDOC"
  else
    why="pandoc $PANDOC_VERSION ($PANDOC) is older than $PANDOC_MIN"
  fi
  echo "NOTE: $why: tables built without the pandoc checks (they need pandoc $PANDOC_MIN or newer; set PANDOC=...)"
  "$PY" scripts/tables/make_tables.py
fi
# the eleven tables of the revision (made from facts.json like the others; fix_widths.py sets their column widths)
for t in rev_proton_tables rev_struct_tables rev_elec_tables fix_widths; do
  "$PY" scripts/tables/$t.py > /dev/null
done
echo "11 revision tables rebuilt"
"$PY" scripts/tables/tables_to_csv.py > /dev/null

echo "== 6. every rebuilt file must be byte-identical to the deposit"
sha256sum -c --quiet SHA256SUMS
echo "SHA256SUMS: all files OK"
