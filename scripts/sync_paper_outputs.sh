#!/bin/sh
# Copy the final paper outputs of the article from the project workspace into this repository.
# Run it once, after the revised manuscript and its supplementary information are final, then review and commit by
# hand. It never commits, never deletes a file here, and never copies a CIF.
#
#   sh scripts/sync_paper_outputs.sh /path/to/project-workspace
#
# Project path (under the workspace)                    -> here
#   papers/btpc/data/facts.json, facts.md                -> analysis/
#   papers/btpc/data/derived/*.csv                       -> analysis/derived/
#   papers/btpc/tables/tbl*.md                           -> tables/
#   papers/btpc/figures/*.png, *.pdf                     -> figures/
#   papers/btpc/si/BTPC_Supplementary_Information.pdf    -> si/
#   papers/btpc/data/build_facts.py                      -> scripts/analysis/build_facts.py
#   papers/btpc/data/rev/*.py                            -> scripts/analysis/rev/
#
# The copies deposited so far were adapted for this repository (docs/provenance.md): build_facts.py reads the
# deposited files through scripts/paths.py (rule P8) and writes no host names or internal file names into
# facts.json (rules F1-F16); facts.md and tables/INDEX.md carry rules M1-M5 and P9. A plain copy undoes those
# changes. After running this script:
#   1. git status; git diff --stat; git diff scripts/analysis/build_facts.py   (what changed, and what was undone)
#   2. carry the adaptations over to the new build_facts.py and rev/ modules, including a mapping of the revision
#      results (project path data/revision_runs/<queue>/out/<job>/ -> calculations/revision1/<job>/) in
#      scripts/paths.py, so that the script runs inside this repository;
#   3. rebuild analysis/facts.json and facts.md with the deposited scripts and compare with the copied register
#      (python scripts/tools/compare_facts.py <copied facts.json> analysis/facts.json); rebuild the tables, their
#      CSV copies and tables/INDEX.md (scripts/tables/make_tables.py --check, scripts/tables/tables_to_csv.py);
#   4. fix every hit of the scan printed below, and check the host names written in facts.json and facts.md (the
#      deposited register describes the compute hosts by their hardware, rules F1-F3); rerun the release gate
#      (docs/REPRODUCE.md), regenerate SHA256SUMS last, and commit. The pre-commit hook checks identity, file
#      types, sizes and wording again.
set -eu
W=${1:?usage: sh scripts/sync_paper_outputs.sh /path/to/project-workspace}
cd "$(dirname "$0")/.."            # repository root
P="$W/papers/btpc"
[ -f "$P/data/facts.json" ] || { echo "not a project workspace: $P/data/facts.json missing" >&2; exit 1; }

LIST=$(mktemp)
trap 'rm -f "$LIST"' EXIT
put() {   # put SOURCE DEST: copy one file, refusing CIF and other crystallographic files
  case "$1" in
    *.cif|*.CIF|*.hkl|*.fcf) echo "refused: $1" >&2; exit 1 ;;
  esac
  [ -f "$1" ] || { echo "missing: $1" >&2; exit 1; }
  mkdir -p "$(dirname "$2")"
  cp "$1" "$2"
  echo "$2" >> "$LIST"
}

put "$P/data/facts.json" analysis/facts.json
put "$P/data/facts.md" analysis/facts.md
for f in "$P"/data/derived/*.csv; do put "$f" "analysis/derived/$(basename "$f")"; done
for f in "$P"/tables/tbl*.md; do put "$f" "tables/$(basename "$f")"; done
for f in "$P"/figures/*.png "$P"/figures/*.pdf; do
  if [ -e "$f" ]; then put "$f" "figures/$(basename "$f")"; fi
done
put "$P/si/BTPC_Supplementary_Information.pdf" si/BTPC_Supplementary_Information.pdf
put "$P/data/build_facts.py" scripts/analysis/build_facts.py
for f in "$P"/data/rev/*.py; do put "$f" "scripts/analysis/rev/$(basename "$f")"; done

echo "copied $(wc -l < "$LIST") files"
echo "== scan of the copied files (absolute paths, cloud host names, internal file names); fix every hit:"
hits=0
while read -r f; do
  out=$(grep -aniE "(^|[[:space:]\"'(=])/(home|root|data|mnt|tmp|scratch)/|[a-z0-9-]+\\.internal|PROGRESS\\.md|HANDOFF|ANALYSIS_AND_PLAN|handover" "$f" || true)
  if [ -n "$out" ]; then
    printf '%s\n' "$out" | cut -c1-160 | sed "s|^|$f:|"
    hits=1
  fi
done < "$LIST"
if [ "$hits" = 0 ]; then echo "(no hits)"; fi
echo "== next: review with git diff, carry over the adaptations (see the header of this script), rebuild, check, commit"
