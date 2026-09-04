#!/bin/bash
#
# Classifies the BND records of a VCF/BCF that is assumed to contain only BNDs.
#
# Every BND record is a single breakend whose ALT encodes the mate position p
# and the orientation of the join, in one of four forms:
#
#   t[p[   3to5      ]p]t   5to3      t]p]   3to3      [p[t   5to5
#
# t[p[ and ]p]t are each other's mate; t]p] and [p[t are self-mated, i.e. both
# breakends of an inverted join use the same form (see symmetrize() in
# BndCanonize.java). Since canonization keeps the endpoint with the smaller
# (chrom,pos), an intra-chromosomal record of a canonized input always has
# p>POS, and its sub-type follows from the ALT form alone:
#
#   DEL-like   t[p[ with p>POS
#   DUP-like   ]p]t with p>POS
#   INV-like   t]p] with p>POS   (3to3, head-to-head)
#          or  [p[t with p>POS   (5to5, tail-to-tail)
#
# For generality each sub-type below also collects the non-canonical mate
# representation, which is the same form with p<POS for the two self-mated
# ones and the other form with p<POS for DEL/DUP:
#
#   DEL-like   t[p[ p>POS  or  ]p]t p<POS
#   DUP-like   ]p]t p>POS  or  t[p[ p<POS
#   INV-like   t]p] or [p[t, either side
#
# These three are mutually exclusive and, together with the degenerate p==POS
# case, cover every intra-chromosomal record; the reconciliation is asserted.
#
INPUT_VCF_GZ=$1
MAIN_CHR_LIST="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

set -euo pipefail

if [ -z "${INPUT_VCF_GZ}" ]; then
    echo "Usage: $0 <input.bcf|input.vcf.gz>"
    exit 1
fi

# %CHROM, %POS: coordinates of this breakend.
# %ALT: bracket notation from which the mate coordinates are parsed.
bcftools query --format '%CHROM\t%POS\t%ALT\n' ${INPUT_VCF_GZ} | awk -v mainChrList="${MAIN_CHR_LIST}" '
BEGIN {
    FS="\t"; OFS="\t";
    n=split(mainChrList,tmp," ");
    for (i=1; i<=n; i++) main[tmp[i]]=1;
}
{
    total++;
    chrom=$1; pos=$2+0; alt=$3;

    # BNDs are single-ALT by construction; if not, only the first ALT is used.
    if (split(alt,altList,",")>1) multiallelic++;
    alt=altList[1];

    # Splitting on brackets: "t[p[" -> ("t","p",""), "[p[t" -> ("","p","t").
    nFields=split(alt,fields,/[][]/);
    if (nFields!=3) { noBracket++; next; }
    spec=fields[2];
    if (!match(spec,/:[0-9]+$/)) { malformed++; next; }
    mateChrom=substr(spec,1,RSTART-1);
    matePos=substr(spec,RSTART+1)+0;

    # The bracket type and the side the base sits on give the orientation.
    match(alt,/[][]/); bracket=substr(alt,RSTART,1);
    if (fields[1]!="") { form=(bracket=="[") ? "t[p[" : "t]p]"; }
    else { form=(bracket=="[") ? "[p[t" : "]p]t"; }

    if (mateChrom in main && chrom in main) {
        bothMain++;
        if (mateChrom==chrom) bothMainIntra++; else bothMainInter++;
    }

    if (mateChrom!=chrom) { inter++; next; }
    intra++;
    if (matePos>pos) side="mateRight"; else if (matePos<pos) side="mateLeft"; else side="mateSame";
    byForm[form " " side]++;
}
END {
    del=byForm["t[p[ mateRight"]+byForm["]p]t mateLeft"];
    dup=byForm["]p]t mateRight"]+byForm["t[p[ mateLeft"];
    inv=byForm["t]p] mateRight"]+byForm["t]p] mateLeft"]+byForm["t]p] mateSame"] \
        +byForm["[p[t mateRight"]+byForm["[p[t mateLeft"]+byForm["[p[t mateSame"];
    # p==POS with a 3to5/5to3 join tells DEL from DUP only by the mate side,
    # so those records cannot be assigned to either.
    samePos=byForm["t[p[ mateSame"]+byForm["]p]t mateSame"];

    printf("total_bnd\t%d\n",total+0);
    printf("unparsable\t%d\n",noBracket+malformed+0);
    printf("inter_chromosomal\t%d\n",inter+0);
    printf("intra_chromosomal\t%d\n",intra+0);
    printf("intra_del_like\t%d\n",del+0);
    printf("intra_dup_like\t%d\n",dup+0);
    printf("intra_inv_like\t%d\n",inv+0);
    printf("intra_ambiguous_same_pos\t%d\n",samePos+0);
    printf("both_breakpoints_main_chr\t%d\n",bothMain+0);
    printf("both_breakpoints_main_chr_intra\t%d\n",bothMainIntra+0);
    printf("both_breakpoints_main_chr_inter\t%d\n",bothMainInter+0);

    # Which ALT forms the intra-chromosomal counts above actually came from: on
    # a canonized input the non-canonical (mateLeft) cells should all be zero.
    printf("\n# intra-chromosomal breakdown by ALT form and mate side\n");
    split("t[p[ ]p]t t]p] [p[t",forms," ");
    split("mateRight mateLeft mateSame",sides," ");
    for (i=1; i<=4; i++) for (j=1; j<=3; j++) {
        key=forms[i] " " sides[j];
        printf("%s\t%s\t%d\n",forms[i],sides[j],byForm[key]+0);
    }

    if (inter+intra+noBracket+malformed!=total+0)
        printf("\n# WARNING: inter+intra+unparsable (%d) != total_bnd (%d)\n",inter+intra+noBracket+malformed,total+0);
    if (del+dup+inv+samePos!=intra+0)
        printf("\n# WARNING: del+dup+inv+ambiguous (%d) != intra_chromosomal (%d)\n",del+dup+inv+samePos,intra+0);
    if (noBracket+0>0) printf("\n# WARNING: %d records with no bracket in ALT (single breakends or non-BND), not classified\n",noBracket);
    if (malformed+0>0) printf("# WARNING: %d records whose ALT mate spec is not <chrom>:<pos>, not classified\n",malformed);
    if (multiallelic+0>0) printf("# WARNING: %d multiallelic records, only the first ALT used\n",multiallelic);
}
'
