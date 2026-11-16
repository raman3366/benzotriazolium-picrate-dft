#!/bin/bash
# Recreate the software environment on Linux x86_64 (Ubuntu). Everything lands in tools/ (git-ignored)
# except the handful of apt packages below. Run from the repository root:  ./scripts/setup_linux.sh
set -euo pipefail
cd "$(dirname "$0")/.."      # repository root
ROOT="$PWD"
NPROC_BUILD=${NPROC_BUILD:-4}

echo "== 0. system packages (unzip for Multiwfn, xvfb+GL for headless ParaView) =="
need=""
# bzip2 is needed to untar the micromamba download (tar -xj) and is NOT present on a minimal
# Ubuntu 24.04 cloud image -- without it the environment install fails with "bzip2: Cannot exec".
for p in bzip2 xz-utils unzip xvfb libgomp1 libgl1 libglu1-mesa libxrender1 libxcursor1 libxft2 libxinerama1 libsm6 libice6; do
  dpkg -s "$p" >/dev/null 2>&1 || need="$need $p"
done
[ -n "$need" ] && { echo "installing:$need"; sudo apt-get update -qq && sudo apt-get install -y -qq $need; } || echo "all present"

echo "== 1. micromamba + conda env (PySCF with OpenMP, geomeTRIC, RDKit, plotting, python-docx) =="
mkdir -p tools/bin
[ -x tools/bin/micromamba ] || (cd tools && curl -Ls "https://micro.mamba.pm/api/micromamba/linux-64/${MICROMAMBA_VERSION:-latest}" | tar -xj bin/micromamba)
export MAMBA_ROOT_PREFIX="$ROOT/tools/mamba_root"
[ -x tools/env/bin/python ] || tools/bin/micromamba create -y -q -p "$ROOT/tools/env" -c conda-forge \
    python=3.11.16 pyscf=2.14.0 geometric=1.1.1 numpy=2.4.6 scipy=1.17.1 matplotlib=3.11.1 \
    rdkit gemmi pypdf python-docx      # versions as recorded for the reported runs (environment/compute-env.yml)
tools/env/bin/python -c "import pyscf, geometric, gemmi, docx; from pyscf import lib; print('PySCF', pyscf.__version__, '| geomeTRIC', geometric.__version__, '| threads', lib.num_threads())"

echo "== 2. Multiwfn (official Linux binary, no compilation needed) =="
MWFN_ZIP=Multiwfn_2026.9.13_bin_Linux_noGUI.zip
if [ ! -x tools/multiwfn/Multiwfn_noGUI ]; then
  mkdir -p tools/multiwfn
  [ -f "tools/$MWFN_ZIP" ] || curl -Ls -o "tools/$MWFN_ZIP" "http://sobereva.com/multiwfn/misc/$MWFN_ZIP"
  unzip -q -o "tools/$MWFN_ZIP" -d tools/multiwfn_unpack
  d=$(find tools/multiwfn_unpack -maxdepth 2 -name 'Multiwfn*' -type d | head -1)
  [ -n "$d" ] || d=tools/multiwfn_unpack
  cp -r "$d"/* tools/multiwfn/ && rm -rf tools/multiwfn_unpack
  chmod +x tools/multiwfn/Multiwfn_noGUI 2>/dev/null || true
fi
# settings.ini: same hard-won values as the 3ATB run (threads, 2-D grid extension, NCI density screening)
S=tools/multiwfn/settings.ini
sed -i -E "s/^( *nthreads *= *).*/\1 ${MWFN_THREADS:-12}/"  "$S"
sed -i -E "s/^( *aug2D *= *).*/\1 6.0/"                     "$S"
sed -i -E "s/^( *RDG_maxrho *= *).*/\1 0.1/"                "$S"
grep -E "^ *(nthreads|aug2D|RDG_maxrho) *=" "$S"
echo "Multiwfn version: $(echo -e '\nq' | tools/multiwfn/Multiwfn_noGUI 2>/dev/null | grep -oE 'Version [0-9.]+[^ ,]*' | head -1)"

echo "== 3. ParaView 6.1.1 (same version as the 3ATB figures) + xvfb for headless rendering =="
PV=ParaView-6.1.1-MPI-Linux-Python3.12-x86_64
if [ ! -x "tools/$PV/bin/pvbatch" ]; then
  [ -f "tools/$PV.tar.gz" ] || curl -Ls -o "tools/$PV.tar.gz" \
    "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v6.1&type=binary&os=Linux&downloadFile=$PV.tar.gz"
  tar -xzf "tools/$PV.tar.gz" -C tools/
fi
cat > tools/bin/pvbatch <<PVEOF
#!/bin/bash
# headless pvbatch: software GL inside a virtual framebuffer
exec xvfb-run -a --server-args="-screen 0 1920x1440x24" "$ROOT/tools/$PV/bin/pvbatch" "\$@"
PVEOF
chmod +x tools/bin/pvbatch
tools/bin/pvbatch -c "from paraview.simple import *; print('ParaView', GetParaViewVersion(), 'headless OK')"

echo
echo "== done =="
echo "  python : $ROOT/tools/env/bin/python"
echo "  multiwfn: $ROOT/tools/multiwfn/Multiwfn_noGUI"
echo "  pvbatch : $ROOT/tools/bin/pvbatch"
