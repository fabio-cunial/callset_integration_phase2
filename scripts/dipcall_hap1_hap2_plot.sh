#!/bin/bash
#
INPUT_TSV="hprc_y2_males.tsv"

set -euxo pipefail

rm -f *_hap*.csv
while read ROW; do
    ID=$(echo "${ROW}" | cut -f 1)
    
    # Coverage
    FIELD_TO_PLOT="6"
    X_VALUE=$(tail -n 1 ${ID}_x_hap1.txt | cut -f ${FIELD_TO_PLOT})
    Y_VALUE=$(tail -n 1 ${ID}_y_hap1.txt | cut -f ${FIELD_TO_PLOT})
    echo "${X_VALUE},${Y_VALUE}" >> coverage_hap1.csv
    X_VALUE=$(tail -n 1 ${ID}_x_hap2.txt | cut -f ${FIELD_TO_PLOT})
    Y_VALUE=$(tail -n 1 ${ID}_y_hap2.txt | cut -f ${FIELD_TO_PLOT})
    echo "${X_VALUE},${Y_VALUE}" >> coverage_hap2.csv
    
    # Meandepth
    FIELD_TO_PLOT="7"
    X_VALUE=$(tail -n 1 ${ID}_x_hap1.txt | cut -f ${FIELD_TO_PLOT})
    Y_VALUE=$(tail -n 1 ${ID}_y_hap1.txt | cut -f ${FIELD_TO_PLOT})
    echo "${X_VALUE},${Y_VALUE}" >> meandepth_hap1.csv
    X_VALUE=$(tail -n 1 ${ID}_x_hap2.txt | cut -f ${FIELD_TO_PLOT})
    Y_VALUE=$(tail -n 1 ${ID}_y_hap2.txt | cut -f ${FIELD_TO_PLOT})
    echo "${X_VALUE},${Y_VALUE}" >> meandepth_hap2.csv
done < ${INPUT_TSV}
