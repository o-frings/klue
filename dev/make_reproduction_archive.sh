#!/bin/zsh
# =============================================================================
# Build the reproduction archive that accompanies the paper's release.
#
# Produces klue-reproduction-v<VERSION>.zip containing everything the data and
# code availability statement promises -- the simulation drivers, the
# orchestration scripts, the five empirical applications, the result files the
# manuscript reports, the figures, and the exact package tarball used -- so the
# archive can be attached to the GitHub release and deposited at Zenodo.
#
# No restricted-access data are redistributed: test_data_files/ is excluded, and
# the Data section of README_REPRODUCE.md says where each dataset comes from.
# The result files follow the same rule as the public repository
# (dev/publish_filter.zsh): run logs, resume part-files, completion markers and
# development cross-checks are left out, and the output/ subfolders the
# manuscript cites are included. The archive also holds renv.lock (package
# versions, dev/record_versions.R), CITATION.cff and the licence.
#
# Usage:  dev/make_reproduction_archive.sh [version]      # default: DESCRIPTION
# =============================================================================
set -e
cd "$(dirname "$0")/.."
VERSION="${1:-$(awk '/^Version:/{print $2}' klue/DESCRIPTION)}"
STAGE="$(mktemp -d)/klue-reproduction-v${VERSION}"
ZIP="klue-reproduction-v${VERSION}.zip"
mkdir -p "$STAGE"

# --- code -------------------------------------------------------------------
cp README_REPRODUCE.md "$STAGE/"
mkdir -p "$STAGE/studies" "$STAGE/dev" "$STAGE/R"
cp studies/*.R "$STAGE/studies/"
cp dev/*.R dev/*.sh dev/*.zsh "$STAGE/dev/" 2>/dev/null || true
for f in dev/*.md(N); do [[ "${f:t}" == brief-* ]] || cp "$f" "$STAGE/dev/"; done   # brief-*.md: notes of another project
cp -R dev/cpp_experiment "$STAGE/dev/" 2>/dev/null || true
cp R/*.R "$STAGE/R/"

# --- the exact package ------------------------------------------------------
R CMD build --no-build-vignettes klue > /dev/null 2>&1
mv "klue_${VERSION}.tar.gz" "$STAGE/"

# --- results the manuscript reports ----------------------------------------
mkdir -p "$STAGE/output"
source dev/publish_filter.zsh   # keep_output, OUTPUT_SUBDIRS
for f in output/*(.N); do
  if keep_output "${f:t}"; then cp "$f" "$STAGE/output/"; fi
done
for d in $OUTPUT_SUBDIRS; do
  rsync -a --include '*/' --include '*.rds' --include '*.csv' --exclude '*' \
    "output/$d/" "$STAGE/output/$d/"
done
cp fig_*.pdf renv.lock "$STAGE/"
cp klue/CITATION.cff "$STAGE/"
for l in klue_public/LICENSE klue/LICENSE; do
  [[ -f "$l" ]] && { cp "$l" "$STAGE/LICENSE"; break }
done

# --- package it -------------------------------------------------------------
rm -f "$ZIP"
(cd "$(dirname "$STAGE")" && zip -rq "$OLDPWD/$ZIP" "$(basename "$STAGE")" -x '*.DS_Store')
rm -rf "$(dirname "$STAGE")"
echo "built $ZIP  ($(du -h "$ZIP" | cut -f1))"
unzip -l "$ZIP" | tail -1
