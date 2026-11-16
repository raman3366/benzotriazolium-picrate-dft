#!/bin/sh
# Unpack the xz-compressed grids (calculations/<case>/grids/*.xz) next to the .xz files, keeping the .xz files,
# and check every unpacked file against docs/grids.sha256 (SHA-256 of the uncompressed grids).
# The unpacked grids (about 420 MiB) are git-ignored.
set -eu
cd "$(dirname "$0")/.."
for f in calculations/*/grids/*.xz; do
  [ -f "${f%.xz}" ] || xz -dk "$f"
done
sha256sum -c docs/grids.sha256
