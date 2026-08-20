#!/bin/bash
# Warning: this script is AI-generated.
#
# Replaces the old alignments bucket prefix with the new one, in every file whose
# name starts with `workpackage1_` in a given directory. Files are edited in
# place; a `.bak` copy of each modified file is kept.
#
# Usage: ./replace_bam_paths.sh <input_dir>
#
set -euo pipefail

INPUT_DIR=$1

OLD_STRING="gs://fc-secure-8f7d6a20-04ce-40d7-8c88-aececeac3e09/CCS/terra-f6671367/outputs/GRCh38/alignments/"
NEW_STRING="gs://prod-drc-broad/longreads/demo_group/UW/CCS/outputs/GRCh38/alignments/"

if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: not a directory: ${INPUT_DIR}" >&2
    exit 1
fi

shopt -s nullglob
FOUND="0"
for FILE in "${INPUT_DIR}"/workpackage1_*; do
    if [ ! -f "${FILE}" ]; then
        continue
    fi
    FOUND="1"
    if ! grep -qF -- "${OLD_STRING}" "${FILE}"; then
        echo "Unchanged (string not found): ${FILE}"
        continue
    fi
    # `|` as delimiter, since the strings contain `/`. Portable in-place edit
    # (GNU sed and BSD/macOS sed disagree on the syntax of `-i`).
    sed "s|${OLD_STRING}|${NEW_STRING}|g" "${FILE}" > "${FILE}.tmp"
    mv "${FILE}" "${FILE}.bak"
    mv "${FILE}.tmp" "${FILE}"
    echo "Updated: ${FILE}"
done

if [ "${FOUND}" = "0" ]; then
    echo "WARNING: no workpackage1_* files in ${INPUT_DIR}" >&2
fi
