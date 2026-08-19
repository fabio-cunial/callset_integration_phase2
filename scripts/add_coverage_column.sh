#!/bin/bash
# Warning: this script is AI-generated.
#
# For every file whose name starts with `workpackage1_` in a given directory,
# creates a new TSV with format `ID \t number \t ...other fields...`, i.e. the
# original record with a new second column inserted. `ID` is the first column of
# the input file, and `number` is taken from the record of a lookup TSV (with
# format `ID \t number`) that has the same ID. Input files are left untouched;
# output files keep the same basename and are written to a separate directory.
#
# IDs that do not occur in the lookup TSV get `NA` as their number, and are
# reported on stderr.
#
# Usage: ./add_coverage_column.sh <input_dir> <lookup_tsv> [output_dir]
#
set -euo pipefail

INPUT_DIR=$1
LOOKUP_TSV=$2
OUTPUT_DIR=${3:-"${INPUT_DIR}/mapped"}

if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: not a directory: ${INPUT_DIR}" >&2
    exit 1
fi
if [ ! -f "${LOOKUP_TSV}" ]; then
    echo "ERROR: not a file: ${LOOKUP_TSV}" >&2
    exit 1
fi
mkdir -p "${OUTPUT_DIR}"

shopt -s nullglob
FOUND="0"
for FILE in "${INPUT_DIR}"/workpackage1_*; do
    if [ ! -f "${FILE}" ]; then
        continue
    fi
    FOUND="1"
    OUTPUT_FILE="${OUTPUT_DIR}/$(basename "${FILE}")"
    # Pass 1 (FNR==NR) loads the lookup table; pass 2 emits one line per input
    # line, in the original order. `-F '\t'` on both, so that IDs containing
    # spaces are handled correctly.
    awk -F '\t' -v OFS='\t' -v FILENAME_IN="${FILE}" '
        FNR==NR { if (NF>=2) NUMBER[$1]=$2; next }
        {
            if ($1 in NUMBER) VALUE=NUMBER[$1];
            else { VALUE="NA"; MISSING++ }
            REST="";
            for (I=2; I<=NF; I++) REST=REST OFS $I;
            print $1, VALUE REST;
        }
        END {
            if (MISSING>0) printf("WARNING: %d ID(s) not found in the lookup table: %s\n", MISSING, FILENAME_IN) > "/dev/stderr";
        }
    ' "${LOOKUP_TSV}" "${FILE}" > "${OUTPUT_FILE}"
    echo "Created: ${OUTPUT_FILE}"
done

if [ "${FOUND}" = "0" ]; then
    echo "WARNING: no workpackage1_* files in ${INPUT_DIR}" >&2
fi
