#!/bin/bash
# Launch one case through the pipeline with every thread pool pinned to the same budget.
#
# Why this matters: PySCF's own thread count (DFT_THREADS) and the BLAS thread count are set independently.
# Leaving BLAS unpinned lets each of N PySCF threads spawn its own full-width BLAS pool, which oversubscribes
# the machine and slows a run by an order of magnitude (observed: a 26-atom SCF that takes 3 min pinned did
# not finish in 25 min unpinned).  Always go through this script.
#
# Usage: ./scripts/run_case.sh <case> [threads] [-- extra pipeline args]   (results go to results/<case>/)
set -euo pipefail
cd "$(dirname "$0")/.."      # repository root
CASE=$1; T=${2:-4}; shift 2 || shift 1 || true
[ "${1:-}" = "--" ] && shift || true

export DFT_THREADS=$T
export OMP_NUM_THREADS=$T MKL_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T NUMEXPR_NUM_THREADS=$T VECLIB_MAXIMUM_THREADS=$T
export DFT_MAX_MEMORY_MB=${DFT_MAX_MEMORY_MB:-6000}

# PySCF spills density-fitting intermediates to scratch -- 12 GB for a 853-basis-function gradient.
# The default /tmp on the compute VM was a 16 GB tmpfs (RAM!), so leaving it there both fails with
# "Disk quota exceeded" and eats the memory the calculation needs.  Put scratch on the real disk.
# Scratch MUST land on the big fast volume: an analytic Hessian for BTP wrote 275 GB, far more than the
# 240 GB boot disk holds.  /mnt/scratch was the local NVMe of the 48-vCPU compute VM; set SCRATCH_ROOT for another volume.
SCRATCH_ROOT=${SCRATCH_ROOT:-/mnt/scratch}
export SCRATCH_ROOT
if [ -d "$SCRATCH_ROOT" ]; then
  export PYSCF_TMPDIR=${PYSCF_TMPDIR:-$SCRATCH_ROOT/$CASE}
else
  export PYSCF_TMPDIR=${PYSCF_TMPDIR:-$PWD/scratch/$CASE}
fi
export TMPDIR=$PYSCF_TMPDIR
mkdir -p "$PYSCF_TMPDIR"
trap 'rm -rf "$PYSCF_TMPDIR"' EXIT
# Integration grid and PCM cavity quality.  These are NOT the library defaults (level 3 / Lebedev 17):
# every case in this project was computed at level 4 (~Gaussian UltraFine) and Lebedev order 29, and they
# used to be exported by hand on each launch.  A re-run launched without them silently dropped to the
# defaults and produced an SCF energy 0.6 mHa away from the rest of the project -- i.e. a different level
# of theory in the same table.  Defaulting them here makes that impossible; override deliberately.
export GRID_LEVEL=${GRID_LEVEL:-4}
export PCM_LEBEDEV_ORDER=${PCM_LEBEDEV_ORDER:-29}
export MWFN_THREADS=${MWFN_THREADS:-$T}
export PVBATCH=${PVBATCH:-$PWD/tools/bin/pvbatch}
export MWFN_BIN=${MWFN_BIN:-$PWD/tools/multiwfn/Multiwfn_noGUI}
# Frequency stage: the worker pool sizes itself from the cores that are actually idle when the
# machine-wide lock is acquired (other cases may still be optimising).  Set FD_WORKERS to override.
export FD_WORKER_THREADS=${FD_WORKER_THREADS:-2}
export FD_RESERVE_CORES=${FD_RESERVE_CORES:-2}

mkdir -p logs
echo "[run_case] $CASE  threads=$T  mem=${DFT_MAX_MEMORY_MB}MB  fd=${FD_WORKERS:-auto}x${FD_WORKER_THREADS}  $(date -Is)"
exec tools/env/bin/python scripts/dft/pipeline.py "$CASE" "$@"
