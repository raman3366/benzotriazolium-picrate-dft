# Environment of the revision calculations

One folder per compute VM (`docs/revision.md`, "How the jobs were run"):

| File | Contents |
|---|---|
| `env.lock.txt` | Explicit conda lock of the environment (package URLs with MD5), written by `scripts/revision1/setup_env.sh` after installing; `setup_env.sh <dir> --lock env.lock.txt` installs exactly these packages (Linux x86-64) |
| `env.versions.txt` | The packages that matter with their build strings (PySCF 2.14.0 `py311h990ede3_1302`, geomeTRIC 1.1.1, NumPy 2.4.6, SciPy 1.17.1, h5py 3.16.0, pyscf-dispersion 1.5.0, libxc 7.1.2, OpenBLAS 0.3.34, Python 3.11.16), the CPU, and the D3(BJ) check of the installer |
| `setup.log` | The installer's log, including the download and SHA-256 of Multiwfn 2026.11.13 (first VM) |

The second VM was installed from the first VM's lock, and its `env.lock.txt` is byte-identical to the first one's.
matplotlib 3.11.2 was added to the first VM's environment later (22 Nov 2026, for plots only) and is not in the
lock.
The locks of the two local machines on which the environment was first tested are in `scripts/revision1/locks/`.
In these files `<vm-1>`, `<vm-2>` stand for the VM host names and `<run-dir>` for the working folder on the VM.
