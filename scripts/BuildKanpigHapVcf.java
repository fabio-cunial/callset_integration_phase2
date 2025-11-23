import java.util.Arrays;
import java.io.*;


/**
 * Given all and only the VCF records that correspond to a kanpig window over 
 * multiple samples, the program outputs a VCF where the entire reference
 * haplotype in the window is the REF of every record, and every distinct non-
 * ref phased haplotype (from any sample) is the ALT of a corresponding record.
 *
 * Remark: every input VCF record is assumed to be phased. If a record is not
 * phased, it is assumed to be phased and it is taken in the given unphased
 * order.
 *
 * Remark: multiple active INS at the same position are taken in the order in
 * which they appear in the VCF.
 *
 * Remark: VCF records may (correctly) overlap by one position on the same hap, 
 * and this is handled by the program. If a VCF record is found to overlap with 
 * a previous record by more than one position on the same hap, it is discarded 
 * and the program continues (i.e. colliding records are removed greedily 
 * rather than optimally). This can happen even if the single-sample VCFs were
 * re-genotyped with kanpig (which disallows collisions), since truvari collapse
 * may create overlaps when choosing cluster representatives.
 *
 * Example of a collision:
 *
 * chr6	89548	id1	GTACATGGAGGGGAACAACACACACCAGGGCCTCTCAGCGGGACAGGGGGTAGGAGACCATCAGGACAAACACGTGGA	G	3462	PASS	NumCollapsed=59;NumConsolidated=199;SUPP_PAV=1;SUPP_PBSV=1;SUPP_SNIFFLES=1;SVTYPE=DEL;SVLEN=77	GT:FT:SQ:GQ:PS:NE:DP:AD:KS:SCORE:CALIBRATION_SENSITIVITY:SUPP_PBSV:SUPP_SNIFFLES:SUPP_PAV	0|1:0:100:25:.:88971:20:8,12:93,96:0.6514:0.9983:0:1:0
 * chr6	89579	id2	CTCTCAGCGGGACAGGGGGTAGGAGACCATCAGGACAAACACGTGGATACATGGAGGGGAACAACACACACCAGGGCCTCTCAGCGGGACAGGGGGTAGGAGACCATCAGGACAAACACGTGGGTACATGGAGGGGAACAACACACACCAGGGCCTCTCAGGGGGACAGGGGGTAGGAGACCATCAGGACAAACACGTGGGTACATGGAGGGGAACAACACACACCAGGGCCTCTCAGGGGGACAGGGGTAGGAGACCATCAGGACAAACACGTGGGTACATGGAGGGCAACAACACACACCAGGGCCTCTCAGGGGGACAGGGGGTAGGAGACCATCAGTACAAACACGTGGATACATGGAGGGGAACAGCACACACCAGGGCCTCTCAGCGGGACAGGGGTAGGAGACCATCAGGACAAACACGTGGGTACATGGAGGGGAACAACACACACCAGGGCC	C	196	PASS	SUPP_PAV=0;SUPP_PBSV=1;SUPP_SNIFFLES=0;SVTYPE=DEL;SVLEN=460;NumCollapsed=17;NumConsolidated=51	GT:FT:SQ:GQ:PS:NE:DP:AD:KS:SCORE:CALIBRATION_SENSITIVITY:SUPP_PBSV:SUPP_SNIFFLES:SUPP_PAV	0|1:0:100:25:.:88971:20:8,12:93,96:0.7083:0.9966:0:1:0
 *
 * Remark: if all the VCF records of a window are inactive in all samples, the
 * output VCF is created but it is empty.
 */
public class BuildKanpigHapVcf {
    
    /**
     * Example usage: BuildKanpigHapVcf list.txt chr6.fa 12680 ./out_dir
     * 
     * Output files: 
     * - `CHUNKID_records.vcf`: the input VCF, without header, FORMAT, SAMPLE;
     * - `CHUNKID_haps.vcf`: the headerless haplotypes VCF, where all POS,REF 
     *   fields are the same, and with no FORMAT or SAMPLE;
     * - `CHUNKID_map.csv`: maps each record of `CHUNKID_haps.vcf` to the 
     *   records of `CHUNKID_records.vcf` that compose it.
     *
     * @param args 
     * 0: a list of files (one per row), each being a headerless, uncompressed
     *    VCF that contains all the records in a truvari-collapsed VCF, limited
     *    to one kanpig window;
     * 2: total number of samples in the cohort;
     * 5: 1=tests the program: assumes that the input VCF contains only one 
     *    autosome and one sample, and prints GTs in the output VCF. The test
     *    must show that `input.vcf > BuildKanpigHapVcf > ConvertKanpigHapVcf > 
     *    output.vcf` is such that `input.vcf` and `output.vcf` have the same 
     *    records in the same order, and the same GTs modulo '/', '|', '.' and
     *    side of '|'.
     */
    public static void main(String[] args) throws IOException {
        final String CHUNKS_FILE = args[0];
        final String CHROMOSOME_FA = args[1];
        final int N_SAMPLES = Integer.parseInt(args[2]);
        final String OUTPUT_DIR = args[3];
        final boolean VERBOSE = args[4].equals("1");
        final boolean DEBUG_MODE = args[5].equals("1");
        
        boolean debugModeFlag;
        int i, j;
        int lastDistinctHap, idGenerator, pos, refLength, altLength, svtype, lastRecord, nRecords;
        int refFirst, refLast;  // Zero-based, inclusive.
        long nOverlaps, nOverlapsTolerated;
        String str1, str2, refChrom, chunkID, refSequence;
        StringBuilder chromosome, tmpBuffer;
        BufferedReader br1, br2;
        BufferedWriter bw1, bw2;
        Haplotype tmpHap;
        int[] lastType1, lastType2;  // 0=undefined, 1=DEL, 2=INS, 3=REPL.
        int[] lastPos1, lastPos2;  // Zero-based, inclusive.
        boolean[][] hap1records, hap2records, matrix; 
        String[] tokens;
        Haplotype[] distinctHaps;
        StringBuilder[] hap1, hap2;
        
        // Loading the reference
        System.err.print("Loading the reference... ");
        chromosome = new StringBuilder();
        br1 = new BufferedReader(new FileReader(CHROMOSOME_FA));
        str1=br1.readLine(); str1=br1.readLine();  // Skipping header
        while (str1!=null) {
            chromosome.append(str1);
            str1=br1.readLine();
        }
        br1.close();
        System.err.println("done");
        
        // Allocating memory
        System.err.print("Allocating memory for "+(2*N_SAMPLES)+" haplotypes... ");
        lastPos1 = new int[N_SAMPLES]; lastPos2 = new int[N_SAMPLES];
        lastType1 = new int[N_SAMPLES]; lastType2 = new int[N_SAMPLES];
        hap1 = new StringBuilder[N_SAMPLES];
        for (i=0; i<N_SAMPLES; i++) hap1[i] = new StringBuilder();
        hap2 = new StringBuilder[N_SAMPLES];
        for (i=0; i<N_SAMPLES; i++) hap2[i] = new StringBuilder();
        distinctHaps = new Haplotype[2*N_SAMPLES];
        for (i=0; i<distinctHaps.length; i++) distinctHaps[i] = new Haplotype(false,-1,null);
        hap1records=null; hap2records=null;
        System.err.println("done");
        
        // Translating chunks
        nOverlaps=0; nOverlapsTolerated=0;
        br1 = new BufferedReader(new FileReader(CHUNKS_FILE));
        str1=br1.readLine();
        while (str1!=null) {  // For each chunk
            chunkID=str1.substring(str1.lastIndexOf("/")+1,str1.lastIndexOf(".vcf"));
            
            // Computing `refChrom,refFirst,refLast`, and outputting the input
            // VCF without samples.
            nRecords=0; refChrom=""; refFirst=Integer.MAX_VALUE; refLast=-1;
            bw1 = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+chunkID+"_records.vcf"));
            br2 = new BufferedReader(new FileReader(str1));
            str2=br2.readLine();        
            while (str2!=null) {
                nRecords++;
                tokens=str2.split("\t");
                if (refChrom.length()==0) refChrom=tokens[0];
                pos=Integer.parseInt(tokens[1]);  // One-based
                if (pos-1<refFirst) refFirst=pos-1;
                refLength=tokens[3].length();
                if (pos-1+refLength-1>refLast) refLast=pos-1+refLength-1;
                bw1.write(tokens[0]);
                for (i=1; i<=7; i++) bw1.write("\t"+tokens[i]);
                bw1.newLine();
                str2=br2.readLine();
            }
            br2.close(); bw1.close();
            if (refFirst==Integer.MAX_VALUE) {
                System.err.println("WARNING: "+str1+" is empty.");
                str1=br1.readLine(); continue;
            }
            if (hap1records==null || hap1records.length<nRecords) hap1records = new boolean[nRecords][N_SAMPLES];
            if (hap2records==null || hap2records.length<nRecords) hap2records = new boolean[nRecords][N_SAMPLES];
            
            // Building the new REF and new distinct ALT haplotypes
            Arrays.fill(lastPos1,refFirst-1); Arrays.fill(lastPos2,refFirst-1);
            Arrays.fill(lastType1,0); Arrays.fill(lastType2,0);
            for (i=0; i<nRecords; i++) Arrays.fill(hap1records[i],false);
            for (i=0; i<nRecords; i++) Arrays.fill(hap2records[i],false);
            br2 = new BufferedReader(new FileReader(str1));
            str2=br2.readLine(); lastRecord=-1;
            while (str2!=null) {
                lastRecord++;
                tokens=str2.split("\t");
                pos=Integer.parseInt(tokens[1]);  // One-based
                refLength=tokens[3].length();
                altLength=tokens[4].length();
                if (refLength>1 && altLength==1) svtype=1;
                else if (refLength==1 && altLength>1) svtype=2;
                else svtype=3;
                for (i=0; i<N_SAMPLES; i++) {
                    if (tokens[9+i].charAt(0)=='1') {
                        if (pos-1<lastPos1[i]) {
                            nOverlaps++;
                            if (svtype==1 && svtype==lastType1[i]) {
                                // Merging two colliding DELs
                                lastPos1[i]=pos-1+refLength-1;
                                nOverlapsTolerated++;
                                hap1records[lastRecord][i]=true;
                            }
                            if (VERBOSE) overlappingError(str1,i,1);
                        }
                        else {
                            if (pos-1>=lastPos1[i]+2) hap1[i].append(chromosome.substring(lastPos1[i]+1,pos-1));
                            if (pos-1==lastPos1[i]) hap1[i].append(tokens[4].substring(1));
                            else hap1[i].append(tokens[4]);
                            lastPos1[i]=pos-1+refLength-1;
                            lastType1[i]=svtype;
                            hap1records[lastRecord][i]=true;
                        }
                    }
                    if (tokens[9+i].charAt(2)=='1') {
                        if (pos-1<lastPos2[i]) {
                            nOverlaps++;
                            if (svtype==1 && svtype==lastType2[i]) {
                                // Merging two colliding DELs
                                lastPos2[i]=pos-1+refLength-1;
                                nOverlapsTolerated++;
                                hap2records[lastRecord][i]=true;
                            }
                            if (VERBOSE) overlappingError(str1,i,2);
                        }
                        else {
                            if (pos-1>=lastPos2[i]+2) hap2[i].append(chromosome.substring(lastPos2[i]+1,pos-1));
                            if (pos-1==lastPos2[i]) hap2[i].append(tokens[4].substring(1));
                            else hap2[i].append(tokens[4]);
                            lastPos2[i]=pos-1+refLength-1;
                            lastType2[i]=svtype;
                            hap2records[lastRecord][i]=true;
                        }
                    }
                }
                str2=br2.readLine();
            }
            br2.close();
            for (i=0; i<N_SAMPLES; i++) {
                if (lastPos1[i]<refLast) hap1[i].append(chromosome.substring(lastPos1[i]+1,refLast+1));
                if (lastPos2[i]<refLast) hap2[i].append(chromosome.substring(lastPos2[i]+1,refLast+1));
            }
            System.err.println("All haps of all samples built");
            
            // Keeping only distinct haplotypes.
            // Remark: if a haplotype can be represented by multiple distinct
            // sets of records, only one such set is selected arbitrarily.
            j=-1;
            for (i=0; i<N_SAMPLES; i++) {
                j++;
                distinctHaps[j].matrixID=false;
                distinctHaps[j].sampleID=i;
                distinctHaps[j].sequence=hap1[i].toString();
                distinctHaps[j].count=1;
            }
            for (i=0; i<N_SAMPLES; i++) {
                j++;
                distinctHaps[j].matrixID=true;
                distinctHaps[j].sampleID=i;
                distinctHaps[j].sequence=hap2[i].toString();
                distinctHaps[j].count=1;
            }
            Arrays.parallelSort(distinctHaps);
            lastDistinctHap=0;
            for (i=1; i<2*N_SAMPLES; i++) {
                if (!distinctHaps[i].equals(distinctHaps[lastDistinctHap])) {
                    lastDistinctHap++;
                    tmpHap=distinctHaps[lastDistinctHap];
                    distinctHaps[lastDistinctHap]=distinctHaps[i];
                    distinctHaps[i]=tmpHap;
                }
                else distinctHaps[lastDistinctHap].count++;
            }
            System.err.println("Distinct haps: "+(lastDistinctHap+1));
            
            // Outputting the haps VCF and the hap->records map
            idGenerator=0;
            bw1 = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+chunkID+"_haps.vcf"));
            bw2 = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+chunkID+"_map.csv"));
            refSequence=chromosome.substring(refFirst,refLast+1);
            debugModeFlag=false;
            for (i=0; i<=lastDistinctHap; i++) {
                if (distinctHaps[i].sequence.equalsIgnoreCase(refSequence)) continue;
                bw1.write(refChrom+"\t"+(refFirst+1)+"\t"+chunkID+"-"+(idGenerator++)+"\t"+refSequence+"\t"+distinctHaps[i].sequence+"\t100\tPASS\t.");
                if (DEBUG_MODE) {
                    if (distinctHaps[i].count==2) bw1.write("\tGT\t1|1");
                    else {
                        if (!debugModeFlag) { bw1.write("\tGT\t0|1"); debugModeFlag=true; }
                        else bw1.write("\tGT\t1|0");
                    }
                }
                bw1.newLine();
                matrix=distinctHaps[i].matrixID?hap2records:hap1records;
                for (j=0; j<nRecords; j++) {
                    if (matrix[j][distinctHaps[i].sampleID]) bw2.write(j+",");
                }
                bw2.newLine();
            }
            bw1.close(); bw2.close();
            
            // Next chunk
            System.err.println("Chunk "+chunkID+" completed");
            for (i=0; i<N_SAMPLES; i++) hap1[i].delete(0,hap1[i].length());
            for (i=0; i<N_SAMPLES; i++) hap2[i].delete(0,hap2[i].length());
            str1=br1.readLine();
        }
        br1.close();
        System.err.println("Over all chunks: nOverlaps="+nOverlaps+" nOverlapsTolerated="+nOverlapsTolerated+" ("+((100.0*nOverlapsTolerated)/nOverlaps)+"%)");
    }
    
    
    /**
     * @param hap 1 or 2.
     */
    private static final void overlappingError(String chunkFile, int sampleID, int hap) throws IOException {
        int i;
        String str, out;
        BufferedReader br;
        String[] tokens;
        
        out="ERROR: overlapping records on the same hap"+hap+" of the "+sampleID+"-th sample of the VCF (zero-based).\n";
        br = new BufferedReader(new FileReader(chunkFile));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            if (tokens[9+sampleID].charAt(0)=='1' || tokens[9+sampleID].charAt(2)=='1') {
                for (i=0; i<9; i++) out+=tokens[i]+"\t";
                out+=tokens[9+sampleID]+"\n";
            }
            str=br.readLine();
        }
        br.close();
        System.err.println(out);
    }
    
    
    private static class Haplotype implements Comparable {
        public boolean matrixID;  // false=hap1records, true=hap2records.
        public int sampleID;
        public String sequence;
        public int count;  // Number of instances of the hap in the cohort
        
        public Haplotype(boolean matrixID, int sampleID, String sequence) {
            this.matrixID=matrixID;
            this.sampleID=sampleID;
            this.sequence=sequence;
        }
        
        public boolean equals(Object other) {
            Haplotype otherHaplotype = (Haplotype)other;
            return sequence.equalsIgnoreCase(otherHaplotype.sequence);
        }
        
        public int compareTo(Object other) {
            Haplotype otherHaplotype = (Haplotype)other;
            return sequence.compareTo(otherHaplotype.sequence);
        }
    }
    
}