#!/bin/bash
#
INPUT_BCF=$1  # Assumed to contain just one chromosome


Continental groups!!!! Filter before everything else.



# Building a sparse matrix with just the necessary information and with integers
# instead of sample strings.
bcftools view --header-only ${INPUT_BCF} | tail -n 1 | tr '\t' '\n' | tail -n +10 > samples.txt
bcftools query --format '%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' --include 'GT="alt"' ${INPUT_BCF} | awk -F'\t' -v OFS='\t' '{
    svlen=length($3)-length($2)
    for (i=4; i<=NF; i++) {
        split($i, sample_gt, "=")
        n = split(sample_gt[2], alleles, "[/|]")
        alt_count = 0
        for (j = 1; j <= n; j++) {
            if (alleles[j] + 0 > 0) alt_count++
        }
        $i = sample_gt[1] "=" alt_count
    }
    out = $1 OFS svlen OFS
    for (i = 4; i <= NF; i++) out = out OFS $i
    print out
}' | python3 -c '
import re, sys
with open("samples.txt", "r") as f:
    mapping = {line.rstrip("\r\n"): str(i) for i, line in enumerate(f, 1) if line.rstrip("\r\n")}
sorted_keys = sorted(mapping.keys(), key=len, reverse=True)
pattern = re.compile("|".join(re.escape(k) for k in sorted_keys))
for line in sys.stdin:
    sys.stdout.write(pattern.sub(lambda m: mapping[m.group(0)], line))
' | gzip -1 -c > matrix.tsv.gz





 

