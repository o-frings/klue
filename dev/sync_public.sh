#!/bin/zsh
# =============================================================================
# Mirror the working repository into the public clone (klue_public/, the
# github.com/o-frings/klue working copy).
#
# Layout published there (established 2026-09-18): the package lives in klue/,
# the paper's reproduction materials at the root. That means install_github and
# r-universe both need the subdir "klue" -- o-frings.r-universe.dev/packages.json
# must carry "subdir": "klue" or the package feed cannot find DESCRIPTION.
#
# What is mirrored:
#   klue/                the package, in full
#   R/ dev/ studies/     the empirical adapters, orchestration scripts, drivers
#   output/              the result files the manuscript reports, with the
#                        subfolders in OUTPUT_SUBDIRS (dev/publish_filter.zsh)
#   README_REPRODUCE.md  the table -> script -> result-file map
#   renv.lock            package versions (dev/record_versions.R)
#   fig_*.pdf            the figures (dev/make_figures.R)
#   CITATION.cff         klue/CITATION.cff, so the root carries the package version
#
# What is not, and why:
#   test_data_files/     no restricted-access data are redistributed
#   output/ files that dev/publish_filter.zsh leaves out: run logs, resume
#                        part-files, development cross-checks, self-tests
#
# This script only writes into klue_public/; it never commits and never pushes.
# Usage:  dev/sync_public.sh [--dry-run]
# =============================================================================
set -e
cd "$(dirname "$0")/.."
PUB="klue_public"
git -C "$PUB" rev-parse --git-dir >/dev/null 2>&1 || { echo "$PUB is not a git clone; aborting" >&2; exit 1 }
DRY=""; [[ "$1" == "--dry-run" ]] && DRY="--dry-run"
RS=(rsync -a --delete $DRY --exclude '.DS_Store' --exclude '.Rproj.user')

# --- the package ------------------------------------------------------------
"${RS[@]}" klue/ "$PUB/klue/"

# --- reproduction code ------------------------------------------------------
"${RS[@]}" --include '*.R' --exclude '*' R/       "$PUB/R/"
"${RS[@]}" --include '*.R' --exclude '*' studies/ "$PUB/studies/"
"${RS[@]}" --exclude 'golden/' --exclude 'brief-*.md' dev/ "$PUB/dev/"
cp $DRY README_REPRODUCE.md "$PUB/" 2>/dev/null || rsync -a $DRY README_REPRODUCE.md "$PUB/"
rsync -a $DRY renv.lock fig_*.pdf "$PUB/"
rsync -a $DRY klue/CITATION.cff "$PUB/CITATION.cff"

# --- results ----------------------------------------------------------------
mkdir -p "$PUB/output"
source dev/publish_filter.zsh   # keep_output, OUTPUT_SUBDIRS
published=()
for f in output/*(.N); do
  b="${f:t}"
  if keep_output "$b"; then published+=("$b"); [[ -n "$DRY" ]] || cp "$f" "$PUB/output/$b"; fi
done
for d in $OUTPUT_SUBDIRS; do
  "${RS[@]}" --include '*/' --include '*.rds' --include '*.csv' --exclude '*' \
    "output/$d/" "$PUB/output/$d/"
done
# drop anything in the public copy that is no longer published
for f in $PUB/output/*(.N); do
  b="${f:t}"
  if [[ ${published[(Ie)$b]} -eq 0 ]]; then
    echo "  removing stale $b"; [[ -n "$DRY" ]] || rm "$f"
  fi
done

echo "synced -> $PUB  (package $(awk '/^Version:/{print $2}' klue/DESCRIPTION), ${#published} result files)"
[[ -n "$DRY" ]] && echo "(dry run: nothing written)"
echo "review with:  git -C $PUB status --short"
