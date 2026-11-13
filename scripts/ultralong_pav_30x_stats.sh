#!/bin/bash
#
ULTRALONG_VCF_GZ=$1


set -euxo pipefail


function CountTsv() {
    local INPUT_TSV=$1

    java GetDiscoveryAnnotations ${INPUT_TSV} 1 0 1 histograms.tsv
    cut -f 2 histograms.tsv > counts.tsv
    rm -f histograms.tsv
    awk -F',' '{ gsub(/\r$/, "") }
    NR == 1 {
        ncols = NF
        for (i = 1; i <= ncols; i++) header[i] = $i
        next
    }
    NF > 0 {
        for (i = 1; i <= ncols; i++) {
            col_sum[i] += $i
        }
        for (i = 2; i <= ncols; i++) {
            total += $i
        }
        subset_sum += $5 + $6 + $7 + $8
    }
    END {
        for (i = 1; i <= ncols; i++) {
            printf "%d: %g\n", header[i], col_sum[i]
        }
        printf "Sum of columns 5,6,7,8    : %g\n", subset_sum
        printf "Sum of all data columns   : %g\n", total
        if (total > 0) {
            printf "Ratio (cols 5-8 / total)  : %.6f (%.2f%%)\n", subset_sum / total, (subset_sum / total) * 100
        } else {
            print "Ratio (cols 5-8 / total)  : N/A (total is 0)"
        }
    }' counts.tsv
    rm -f counts.tsv
}


# ------------------------------ Main program ----------------------------------

bcftools query --format '%ID[\t%SUPP_PBSV,%SUPP_SNIFFLES,%SUPP_PAV,0]\n' --include 'INFO/SVTYPE="DEL"' ${ULTRALONG_VCF_GZ} > del.tsv &
bcftools query --format '%ID[\t%SUPP_PBSV,%SUPP_SNIFFLES,%SUPP_PAV,0]\n' --include 'INFO/SVTYPE="INS" && INSDUP_SVLEN="."' ${ULTRALONG_VCF_GZ} > ins.tsv &
bcftools query --format '%ID[\t%SUPP_PBSV,%SUPP_SNIFFLES,%SUPP_PAV,0]\n' --include 'INFO/SVTYPE="INS" && INSDUP_SVLEN!="."' ${ULTRALONG_VCF_GZ} > insdup.tsv &
wait
CountTsv del.tsv
CountTsv ins.tsv
CountTsv insdup.tsv
