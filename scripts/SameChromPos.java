import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.io.*;

/**
 * Computes basic properties of distinct (CHROM,POS) pairs in a VCF file.
 * Used to investigate the `bcftools annotate -c CHROM,POS,~ID` bug.
 */
public class SameChromPos {
    /**
     * The program prints to STDOUT a CSV row with the following fields:
     * 
     * 0: total number of records;
     * 1: number of records with the same CHROM,POS as another record;
     * 
     * 2: number of distinct CHROM,POS pairs;
     * 3: number of CHROM,POS pairs with an INS and a DEL;
     * 4: number of CHROM,POS pairs with multiple INS;
     * 5: number of CHROM,POS pairs with multiple DEL;
     * 6: number of CHROM,POS pairs with >1 record.
     * 
     * It prints to STDERR a CSV row with a field for every INFO annotation,
     * where the value of the field is the number of distinct CHROM,POS pairs 
     * with non-constant value for that field.
     */
    public static int[] counts = new int[7];

    /**
     * Annotations
     */
    public static final String[] ANNOTATIONS = new String[] {"KS_1", "KS_2", "SQ", "GQ", "DP", "AD_NON_ALT", "AD_ALL", "GT_COUNT", "SUPP_PBSV", "SUPP_SNIFFLES", "SUPP_PAV"};
    public static double[] currentAnnotations = new double[ANNOTATIONS.length];
    public static boolean[] currentAnnotationsDiffer = new boolean[ANNOTATIONS.length];
    public static int[] blocksWithDifferentAnnotations = new int[ANNOTATIONS.length];

    /**
     * @param args 
     * 0: assumed to be sorted by CHROM,POS;
     * 3: a minimal output VCF with one record per distinct (CHROM,POS) pair.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String OUTPUT_TXT_CLUSTER_SIZE = args[1];
        final String OUTPUT_CSV_MAX_DELTAS = args[2];
        final String OUTPUT_VCF_CLUSTERS = args[3];

        int i;
        int svlen, pos, currentPos, currentRecords, currentIns, currentDel;
        double currentMin, currentMax, annotation;
        String str, svtype, currentChrom, value;
        BufferedReader br;
        BufferedWriter bwSize, bwDelta, bwVcf;
        String[] tokens;

        bwSize = new BufferedWriter(new FileWriter(OUTPUT_TXT_CLUSTER_SIZE));
        bwDelta = new BufferedWriter(new FileWriter(OUTPUT_CSV_MAX_DELTAS));
        bwVcf = new BufferedWriter(new FileWriter(OUTPUT_VCF_CLUSTERS));
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        currentChrom=""; currentPos=-1; currentRecords=0; currentIns=0; currentDel=0; currentMin=Integer.MAX_VALUE; currentMax=Integer.MIN_VALUE;
        Arrays.fill(currentAnnotations,Integer.MAX_VALUE);
        Arrays.fill(currentAnnotationsDiffer,false);
        Arrays.fill(blocksWithDifferentAnnotations,0);
        tokens=null; pos=-1;
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                bwVcf.write(str); bwVcf.newLine();
                str=br.readLine();
                continue;
            }
            counts[0]++;
            tokens=str.split("\t");
            pos=Integer.parseInt(tokens[1]);
            if (!tokens[0].equals(currentChrom) || pos!=currentPos) {
                // End of a block
                if (currentChrom.length()>0) {
                    counts[2]++;
                    if (currentRecords>1) {
                        counts[6]++;
                        counts[1]+=currentRecords;
                        if (currentIns>0 && currentDel>0) counts[3]++;
                        if (currentIns>1) counts[4]++;
                        if (currentDel>1) counts[5]++;
                        bwSize.write(currentRecords+"\n");
                        bwDelta.write((currentMax-currentMin)+","+((currentMax-currentMin)/currentMin)+"\n");
                        for (i=0; i<currentAnnotations.length; i++) {
                            if (currentAnnotationsDiffer[i]) blocksWithDifferentAnnotations[i]++;
                        }
                        bwVcf.write(currentChrom+"\t"+currentPos+"\t.\tN\tN\t.\t.\t.\tGT\t0/1\n");
                    }
                }
                currentChrom=tokens[0]; currentPos=pos;
                currentRecords=0; currentIns=0; currentDel=0;
                currentMin=Integer.MAX_VALUE; currentMax=Integer.MIN_VALUE;
                Arrays.fill(currentAnnotations,Integer.MAX_VALUE);
                Arrays.fill(currentAnnotationsDiffer,false);
            }
            currentRecords++;
            svtype=getInfoField(tokens[7],"SVTYPE");
            if (svtype.equalsIgnoreCase("INS")) currentIns++;
            else if (svtype.equalsIgnoreCase("DEL")) currentDel++;
            svlen=Integer.parseInt(getInfoField(tokens[7],"SVLEN"));
            if (svlen<0) svlen=-svlen;
            if (svlen<currentMin) currentMin=svlen;
            if (svlen>currentMax) currentMax=svlen;
            for (i=0; i<ANNOTATIONS.length; i++) {
                value=getInfoField(tokens[7],ANNOTATIONS[i]);
                if (value!=null) {
                    annotation=Double.parseDouble(value);
                    if (currentAnnotations[i]==Integer.MAX_VALUE) currentAnnotations[i]=annotation;
                    else if (currentAnnotations[i]!=annotation) currentAnnotationsDiffer[i]=true;
                }
            }
            str=br.readLine();
        }
        // End of last block
        if (currentChrom.length()>0) {
            counts[2]++;
            if (currentRecords>1) {
                counts[6]++;
                counts[1]+=currentRecords;
                if (currentIns>0 && currentDel>0) counts[3]++;
                if (currentIns>1) counts[4]++;
                if (currentDel>1) counts[5]++;
                bwSize.write(currentRecords+"\n");
                bwDelta.write((currentMax-currentMin)+","+((currentMax-currentMin)/currentMin)+"\n");
                for (i=0; i<currentAnnotations.length; i++) {
                    if (currentAnnotationsDiffer[i]) blocksWithDifferentAnnotations[i]++;
                }
                bwVcf.write(currentChrom+"\t"+currentPos+"\t.\tN\tN\t.\t.\t.\tGT\t0/1\n");
            }
        }
        br.close(); bwSize.close(); bwDelta.close(); bwVcf.close();

        // Outputting
        System.out.print(counts[0]+"");
        for (i=1; i<counts.length; i++) System.out.print(","+counts[i]);
        System.out.println();
        System.err.print(blocksWithDifferentAnnotations[0]+"");
        for (i=1; i<ANNOTATIONS.length; i++) System.err.print(","+blocksWithDifferentAnnotations[i]);
        System.err.println();
    }


    /**
	 * @return NULL if $field$ does not occur in $info$.
	 */
	private static final String getInfoField(String info, String field) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
        
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return null;
		q=info.indexOf(";",p+FIELD_LENGTH);
		return info.substring(p+FIELD_LENGTH,q<0?info.length():q);
	}

}