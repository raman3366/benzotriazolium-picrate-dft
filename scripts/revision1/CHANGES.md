# Changes made to the revision scripts for deposit

This folder is a copy of the project folder `dft/revision/` as it stood after the revision runs (23 Nov 2026), without
`__pycache__/`. In the files that ran and in the job records it is `dft/revision/`; the scripts expect that place in
the project layout (`docs/revision.md`, "Running one job again").

The edits below remove citations of internal notes, the cloud account name, VM host names and absolute paths of the
machines. They change comments, help text and default values only; no computation changes. Every rule was applied by
a script that counted the replacements. Nothing else in the folder differs from the project copy.

| Rule | File | What changed | Count |
|---|---|---|---:|
| S1 | `rev_common.py` | Module docstring: citation of an internal notes file replaced by "project notes, not deposited" | 1 |
| S2 | `setup_env.sh` | Comment: citation of an internal notes file replaced by the script that applies the same two Multiwfn settings (`setup_linux.sh`, deposited as `scripts/setup_linux.sh`) | 1 |
| S3 | `jobs_piserver.json` | `_about` line: citation of an internal notes file replaced by "project notes, not deposited" | 1 |
| S4 | `jobs_vm.json`, `jobs_vm2.json` | `_about` line: VM name and cloud account name replaced by `<vm-1>` / `<vm-2>`; the machine type is kept | 2 |
| S5 | `jobs_piserver.json`, `tests/validate_jobs.yaml`, `locks/linux-aarch64_piserver.versions.txt` | Absolute path of the working folder on the local ARM64 machine replaced by `<run-dir>` | 7 |
| S6 | `locks/linux-64_x86-64-v4_msi.versions.txt` | Absolute path of the working folder on the local x86-64 laptop (in its user's home folder) replaced by `<run-dir>` | 1 |
| S7 | `vm/create_vm.sh` | Default gcloud configuration (named after the cloud account) replaced by `default`; set `CFG` as before | 1 |
| S8 | `vm/create_vm.sh` | Default SSH key path written relative to `$HOME` instead of the controlling machine's root home folder (the same file when run as root) | 1 |
| S9 | `vm/create_vm.sh` | Comment: citation of an internal planning file replaced | 1 |

SHA-256 of the edited files before and after the edits. The job records (`calculations/revision1/<job>/result.json`,
`script_sha256`) give the values before the edits for the scripts that ran; `rev_common.py` is the only one of those
that was edited.

| File | Before (as run) | After (deposited) |
|---|---|---|
| `rev_common.py` | `fba929added9222d92adeca239e64fcb87eaaeefe6511ad436c2624eba72515a` | `940826d7774a6a4de1be6ad7deffd7bbaffac69cec45482238b953a4403df2ab` |
| `setup_env.sh` | `f192474b879eff0a08d5b6324279ab48047934e189335d856c80bde67eb593f5` | `1c0f6b3a4933294a7f22230e9ca0900a991204f2a686999788d9d64d9602576f` |
| `jobs_vm.json` | `7c3acf9bcd8b08f997714318e1203dacf0de1645da873f847e8eb72c021cfe60` | `a8f42be9c2c9888ea0d6d8b7ded1eb03790477f63a712796c740824750a5c57b` |
| `jobs_vm2.json` | `a4b11d3328373a7028132e998d80c2778a793b3fa9a07c0702e9c8c0feb44f79` | `211aea7def0e4e83fe0ec71405a496b105a101f26666c11e59e261b0010f6234` |
| `jobs_piserver.json` | `bb3243faff1c99d6616d2159ac65354bb0f5ec023cdab7bd64c90e333a0dfae1` | `3144849e8598c01ed49330be7e7247c33ea4a73a1cacd895537227f1305aa3ec` |
| `tests/validate_jobs.yaml` | `108a24491540ef5434d6b7b9ce3e3aca6c0b34f4845a8ffd720c23f1f14edc24` | `eff87d448571846204705fb215544d5077056967d4ef2c5cb29a05501b9f348c` |
| `vm/create_vm.sh` | `0a13a8496c999ef2e8dfde67ae3ed9974e0c7c2cf3fb6b75c4138d31fe2d2c07` | `4104e5881c3a9615fe98b2fe63d3690d90e991bdd7b44da051b9511ddcea3dff` |
| `locks/linux-64_x86-64-v4_msi.versions.txt` | `03889f08bb3ca2347519ddfcb9a45d7b6c30ec1bf5ec465e2bfa80c680099559` | `38d458527d7847d7c919bad0c8264374ad6d5b070f75f1a50c48e4e4827316de` |
| `locks/linux-aarch64_piserver.versions.txt` | `161cb0323f5f2879acab3dcfa372ff08b8b94bc172eb05ca208b81b5dbb1dc4a` | `9c4eb2e18a3549e19eebfbc802d56a25f785edf4b1afeb3adfa60cd16cfc9226` |

The edit of `rev_common.py` changed only line 5, the end of the docstring sentence that names the source of the
thread-pinning rule; the code is identical to the version recorded in the job records.

Kept as they are:

- `msi` and `piserver`, the names of the two local machines (a 4-core x86-64 laptop and a 4-core ARM64 board, a
  Raspberry Pi 5), in file names, `host_hint` and `est_hours_*` fields and comments; "the Pi" is the machine that
  controlled the runs.
- The working-folder default `~/r007` (`setup_env.sh`, `charges.py`, the job lists) and the resource names
  `r007-*` in `vm/create_vm.sh`: neutral names of the project, and changing them would change the defaults.
- Comments that mention plan items, review findings, the reviewers and "the user" (the person who ran the queue)
  refer to working notes of the project that are not deposited.
- `jobs_main.json`, `jobs_proton.json` and `jobs_piserver.json` are the earlier plans for the local machines;
  `jobs_vm.json` (first VM) and `jobs_vm2.json` (second VM) are the lists that ran.
