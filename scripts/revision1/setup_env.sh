#!/usr/bin/env bash
# Pinned, idempotent environment installer for the BTPC revision calculations (no sudo).
#
#   bash dft/revision/setup_env.sh [ROOT] [--multiwfn] [--lock FILE]
#
# ROOT (default ~/r007) receives:
#   ROOT/bin/micromamba        micromamba 2.10.0-0, SHA-256 checked
#   ROOT/mamba_root/           micromamba package cache
#   ROOT/env/                  the conda environment (python, PySCF, geomeTRIC, pyscf-dispersion, ...)
#   ROOT/env.lock.txt          explicit lock of what was installed (URLs + md5), written after every run
#   ROOT/env.versions.txt      versions and build strings of the packages that matter, plus the CPU level
#   ROOT/multiwfn/             Multiwfn noGUI (only with --multiwfn), SHA-256 in ROOT/multiwfn/SHA256
#
# Code of record (papers/btpc/data/facts.json, docs/environment/python_env.txt): PySCF 2.14.0 conda-forge build
# with rattler-build timestamp 1789238536 (= pyscf-2.14.0-py311hc37a911_1402, the x86-64 microarch-level-4 build of
# 13 Sep 2026), geomeTRIC 1.1.1, Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1, h5py built 19 Aug 2026 (= 3.16.0).
# conda-forge has PySCF 2.14.0 only for Python 3.10 and 3.11, so the environment is Python 3.11 (the record), not
# 3.12.  On an x86-64 CPU without AVX-512 (e.g. AMD Milan, GCP t2d) the solver takes the level-3 build of the same
# release (py311h990ede3_1302); on ARM64 it takes py311h6b91cb6_102.  The validation (dft/revision/validate.py)
# must be rerun on any new host.
#
# --lock FILE installs exactly the packages of an explicit lock written by an earlier run (same platform and CPU
# level only), instead of solving from the pins below.
set -euo pipefail

ROOT="$HOME/r007"; WANT_MWFN=0; LOCK=""
while [ $# -gt 0 ]; do
  case "$1" in
    --multiwfn) WANT_MWFN=1 ;;
    --lock) LOCK="$2"; shift ;;
    -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
    *) ROOT="$1" ;;
  esac
  shift
done
mkdir -p "$ROOT/bin"
ROOT="$(cd "$ROOT" && pwd)"

# ---------------------------------------------------------------- pins
MICROMAMBA_VERSION=2.10.0-0
declare -A MICROMAMBA_SHA256=(
  [linux-64]=701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd
  [linux-aarch64]=3bb56786e21b9b24024f599fe9185b4a59dc3f9f59937a42311ee583adbbb23c
)
SPECS=(
  "python=3.11.16"
  "pyscf=2.14.0=*02"            # the 13 Sep 2026 builds (_1402/_1302/_1102 x86-64, _102 aarch64); _1402 is the record
  "geometric=1.1.1"
  "numpy=2.4.6"
  "scipy=1.17.1"
  "h5py=3.16.0=nompi_*"
  "pyscf-dispersion=1.5.0"      # conda-forge build of 11 Sep 2026, contemporary with PySCF 2.14.0 (D3(BJ) via simple-dftd3)
  "pyyaml"
  "psutil"
)
# Multiwfn: 2026.9.13 is the version of record; the server keeps only the latest, so fall back to 2026.11.13.
MWFN_VERSIONS=(2026.9.13 2026.11.13)

case "$(uname -m)" in
  x86_64)  PLAT=linux-64 ;;
  aarch64) PLAT=linux-aarch64 ;;
  *) echo "unsupported machine $(uname -m)" >&2; exit 1 ;;
esac
export MAMBA_ROOT_PREFIX="$ROOT/mamba_root"
MM="$ROOT/bin/micromamba"
log() { echo "[setup_env $(date '+%Y-%m-%d %H:%M:%S %Z')] $*"; }

# ---------------------------------------------------------------- 1. micromamba
if [ -x "$MM" ] && [ "$("$MM" --version 2>/dev/null)" = "${MICROMAMBA_VERSION%-*}" ] \
   && [ "$(sha256sum "$MM" | cut -d' ' -f1)" = "${MICROMAMBA_SHA256[$PLAT]}" ]; then
  log "micromamba $MICROMAMBA_VERSION present"
else
  log "downloading micromamba $MICROMAMBA_VERSION ($PLAT)"
  tmp="$(mktemp "$ROOT/bin/micromamba.XXXXXX")"
  curl -fsSL -o "$tmp" \
    "https://github.com/mamba-org/micromamba-releases/releases/download/$MICROMAMBA_VERSION/micromamba-$PLAT"
  got="$(sha256sum "$tmp" | cut -d' ' -f1)"
  if [ "$got" != "${MICROMAMBA_SHA256[$PLAT]}" ]; then
    rm -f "$tmp"; echo "micromamba SHA-256 mismatch: $got" >&2; exit 1
  fi
  chmod +x "$tmp"; mv "$tmp" "$MM"
fi

# ---------------------------------------------------------------- 2. environment
check_env() {
  "$ROOT/env/bin/python" - <<'PY'
import sys, pyscf, geometric, numpy, scipy, h5py
assert pyscf.__version__ == '2.14.0', pyscf.__version__
assert geometric.__version__.startswith('1.1.1'), geometric.__version__
assert numpy.__version__ == '2.4.6', numpy.__version__
assert scipy.__version__ == '1.17.1', scipy.__version__
import pyscf.dispersion  # noqa: F401  (pyscf-dispersion)
print('ok', sys.version.split()[0])
PY
}
if [ -x "$ROOT/env/bin/python" ] && check_env >/dev/null 2>&1; then
  log "environment $ROOT/env present and matches the pins"
elif [ -n "$LOCK" ]; then
  log "creating $ROOT/env from explicit lock $LOCK"
  "$MM" create -y -q -p "$ROOT/env" --file "$LOCK"
else
  log "creating $ROOT/env from the pinned specs: ${SPECS[*]}"
  "$MM" create -y -q -p "$ROOT/env" -c conda-forge --override-channels "${SPECS[@]}"
fi
check_env

"$MM" env export -p "$ROOT/env" --explicit --md5 > "$ROOT/env.lock.txt"
{
  echo "# written by dft/revision/setup_env.sh on $(hostname) at $(date -Is)"
  echo "machine: $(uname -m)  platform: $PLAT  cpu: $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ //')"
  echo "micromamba: $("$MM" --version)"
  "$MM" list -p "$ROOT/env" 2>/dev/null \
    | grep -E '^\s*(python|pyscf|pyscf-dispersion|geometric|numpy|scipy|h5py|libxc-c|libxc|libcint|xcfun|libblas|libopenblas|simple-dftd3|dftd4|_x86_64-microarch-level|llvm-openmp|_openmp_mutex)\s' \
    | awk '{printf "%-28s %-14s %s\n", $1, $2, $3}'
  "$ROOT/env/bin/python" - <<'PY'
import numpy, pyscf
from pyscf import lib
from pyscf.dft import libxc
print('pyscf', pyscf.__version__, 'from', pyscf.__file__)
print('libxc', libxc.libxc_version() if hasattr(libxc, 'libxc_version') else '?')
cfg = numpy.show_config(mode='dicts') if 'mode' in numpy.show_config.__code__.co_varnames else {}
bl = (cfg.get('Build Dependencies', {}) or {}).get('blas', {})
print('numpy blas', bl.get('name'), bl.get('version'))
print('pyscf max threads', lib.num_threads())
PY
} > "$ROOT/env.versions.txt"
cat "$ROOT/env.versions.txt"

# dispersion keywords used by the revision (B3LYP-D3(BJ), wB97X-D3(BJ)) must parse and give a D3 energy
"$ROOT/env/bin/python" - <<'PY'
from pyscf import gto, dft
mol = gto.M(atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587', basis='sto-3g', verbose=0)
for xc in ('b3lyp-d3bj', 'wb97x-d3bj'):
    mf = dft.RKS(mol, xc=xc)
    e = mf.kernel()
    print(f'dispersion check {xc:12s} E = {e:.8f}  disp = {getattr(mf, "disp", None)}  Edisp = {mf.scf_summary.get("dispersion", float("nan")):.3e}')
PY

# ---------------------------------------------------------------- 3. Multiwfn (optional, x86-64 only)
if [ "$WANT_MWFN" = 1 ]; then
  if [ "$PLAT" != linux-64 ]; then
    log "Multiwfn: no official binary for $PLAT; skipped (build from source if needed)"
  elif [ -x "$ROOT/multiwfn/Multiwfn_noGUI" ]; then
    log "Multiwfn present: $(cat "$ROOT/multiwfn/VERSION" 2>/dev/null)"
  else
    mkdir -p "$ROOT/downloads" "$ROOT/multiwfn"
    got=""
    for v in "${MWFN_VERSIONS[@]}"; do
      z="Multiwfn_${v}_bin_Linux_noGUI.zip"
      if [ -f "$ROOT/downloads/$z" ] || curl -fsSL -o "$ROOT/downloads/$z" "http://sobereva.com/multiwfn/misc/$z"; then
        got="$v"; break
      fi
      rm -f "$ROOT/downloads/$z"; log "Multiwfn $v not downloadable"
    done
    [ -n "$got" ] || { echo "no Multiwfn version could be downloaded" >&2; exit 1; }
    z="Multiwfn_${got}_bin_Linux_noGUI.zip"
    sha256sum "$ROOT/downloads/$z" | tee "$ROOT/multiwfn/SHA256"
    rm -rf "$ROOT/downloads/unpack"; mkdir -p "$ROOT/downloads/unpack"
    "$ROOT/env/bin/python" -I -c "import sys, zipfile; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" \
      "$ROOT/downloads/$z" "$ROOT/downloads/unpack"
    d="$(dirname "$(find "$ROOT/downloads/unpack" -name Multiwfn_noGUI -type f | head -1)")"
    cp -r "$d"/. "$ROOT/multiwfn/" && rm -rf "$ROOT/downloads/unpack"
    chmod +x "$ROOT/multiwfn/Multiwfn_noGUI"
    echo "$got (http://sobereva.com/multiwfn/misc/$z, downloaded $(date -Is))" > "$ROOT/multiwfn/VERSION"
    # the project's settings (as in setup_linux.sh): 2-D grid extension and NCI density screening
    S="$ROOT/multiwfn/settings.ini"
    sed -i -E "s/^( *aug2D *= *)[^ ]*/\1 6.0/; s/^( *RDG_maxrho *= *)[^ ]*/\1 0.1/" "$S"
    grep -E "^ *(nthreads|aug2D|RDG_maxrho) *=" "$S" | cut -c1-60
  fi
fi
log "done: python $ROOT/env/bin/python"
