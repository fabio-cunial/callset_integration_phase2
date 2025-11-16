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
 * Remark: the program handles VCF records that overlap by one position on the
 * same hap. It stops with an error if two VCF records overlap by more than one
 * position on the same hap.
 */
public class BuildKanpigHapVcf {
    
    /**
     * Example usage: BuildKanpigHapVcf list.txt chr6.fa 12680 ./out_dir
     * 
     * Remark: the program writes in output a headerless VCF where all POS,REF 
     * fields are the same, and with no FORMAT or SAMPLE columns.
     *
     * @param args 
     * 0: a list of files (one per row), each being a headerless, uncompressed
     *    VCF that contains all the records in a truvari-collapsed VCF, limited
     *    to one kanpig window;
     * 2: total number of samples in the cohort.
     */
    public static void main(String[] args) throws IOException {
        final String CHUNKS_FILE = args[0];
        final String CHROMOSOME_FA = args[1];
        final int N_SAMPLES = Integer.parseInt(args[2]);
        final String OUTPUT_DIR = args[3];
        
        int i, j;
        int lastDistinctHap, idGenerator, pos, refLength;
        int refFirst, refLast;  // Zero-based, inclusive.
        String str1, str2, refChrom, chunkID, refSequence;
        StringBuilder chromosome, tmpBuffer;
        BufferedReader br1, br2;
        BufferedWriter bw;
        int[] lastPos1, lastPos2;  // Zero-based, inclusive.
        String[] tokens, distinctHaps;
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
        hap1 = new StringBuilder[N_SAMPLES];
        for (i=0; i<N_SAMPLES; i++) hap1[i] = new StringBuilder();
        hap2 = new StringBuilder[N_SAMPLES];
        for (i=0; i<N_SAMPLES; i++) hap2[i] = new StringBuilder();
        distinctHaps = new String[2*N_SAMPLES];
        System.err.println("done");
        
        // Translating chunks
        br1 = new BufferedReader(new FileReader(CHUNKS_FILE));
        str1=br1.readLine();
        while (str1!=null) {  // For each chunk     
            // Computing `refChrom,refFirst,refLast`.
            refChrom=""; refFirst=Integer.MAX_VALUE; refLast=-1;
            br2 = new BufferedReader(new FileReader(str1));
            str2=br2.readLine();        
            while (str2!=null) {                
                tokens=str2.split("\t");
                if (refChrom.length()==0) refChrom=tokens[0];
                pos=Integer.parseInt(tokens[1]);  // One-based
                if (pos-1<refFirst) refFirst=pos-1;
                refLength=tokens[3].length();
                if (pos-1+refLength-1>refLast) refLast=pos-1+refLength-1;
                str2=br2.readLine();
            }
            br2.close();
            if (refFirst==Integer.MAX_VALUE) {
                System.err.println("WARNING: "+str1+" is empty.");
                str1=br1.readLine(); continue;
            }
            
            // Building new REF, and new ALT haplotypes for each sample.
            Arrays.fill(lastPos1,refFirst-1); Arrays.fill(lastPos2,refFirst-1);
System.err.println("VITTU> 0  lastPos1="+lastPos1[0]+" lastPos2="+lastPos2[0]);                        
            br2 = new BufferedReader(new FileReader(str1));
            str2=br2.readLine(); 
            while (str2!=null) {
                tokens=str2.split("\t");
                pos=Integer.parseInt(tokens[1]);  // One-based
                refLength=tokens[3].length();
                for (i=0; i<N_SAMPLES; i++) {
System.err.println("VITTU> 1  "+tokens[9+i].charAt(0)+" / "+tokens[9+i].charAt(2));                    
                    if (tokens[9+i].charAt(0)=='1') {
                        if (pos-1>=lastPos1[i]+2) hap1[i].append(chromosome.substring(lastPos1[i]+1,pos-1));
                        if (pos-1<lastPos1[i]) {
                            System.err.println("ERROR: overlapping records on the same hap1.");
                            System.exit(1);
                        }
                        else if (pos-1==lastPos1[i]) hap1[i].append(tokens[4].substring(1));
                        else hap1[i].append(tokens[4]);
                        lastPos1[i]=pos-1+refLength-1;
                    }
                    if (tokens[9+i].charAt(2)=='1') {
                        if (pos-1>=lastPos2[i]+2) hap2[i].append(chromosome.substring(lastPos2[i]+1,pos-1));
                        if (pos-1<lastPos2[i]) {
                            System.err.println("ERROR: overlapping records on the same hap2.");
                            System.exit(1);
                        }
                        else if (pos-1==lastPos2[i]) hap2[i].append(tokens[4].substring(1));
                        else hap2[i].append(tokens[4]);
                        lastPos2[i]=pos-1+refLength-1;
                    }
                }
                str2=br2.readLine();
            }
            br2.close();
System.err.println("VITTU> 2  hap1="+hap1[0]+" hap2="+hap2[0]);            
System.err.println("VITTU> 2  lastPos1="+lastPos1[0]+" lastPos2="+lastPos2[0]);
            for (i=0; i<N_SAMPLES; i++) {
                if (lastPos1[i]<refLast) hap1[i].append(chromosome.substring(lastPos1[i]+1,refLast+1));
                if (lastPos2[i]<refLast) hap2[i].append(chromosome.substring(lastPos2[i]+1,refLast+1));
            }
            System.err.println("All haps of all samples built");
System.err.println("VITTU> 3  hap1="+hap1[0]+" hap2="+hap2[0]);               
            
            // Keeping only distinct haplotypes
            j=-1;
            for (i=0; i<N_SAMPLES; i++) distinctHaps[++j]=hap1[i].toString();
            for (i=0; i<N_SAMPLES; i++) distinctHaps[++j]=hap2[i].toString();

System.err.println("VITTU> 3.5  distinctHaps before sorting:");
for (int x=0; x<=j; x++) System.err.println(distinctHaps[x]); 

            Arrays.parallelSort(distinctHaps);
            
System.err.println("VITTU> 4  sorted haps:");
for (int x=0; x<=j; x++) System.err.println(distinctHaps[x]);            
            
            lastDistinctHap=0;
            for (i=1; i<2*N_SAMPLES; i++) {
                if (!distinctHaps[i].equalsIgnoreCase(distinctHaps[lastDistinctHap])) distinctHaps[++lastDistinctHap]=distinctHaps[i];
            }
            
System.err.println("VITTU> 5  distinct haps:");
for (int x=0; x<=lastDistinctHap; x++) System.err.println(distinctHaps[x]);                        
            
            
            System.err.println("Distinct haps: "+(lastDistinctHap+1));            
            
            // Writing the output chunk
            idGenerator=0;
            chunkID=str1.substring(str1.lastIndexOf("/")+1,str1.lastIndexOf(".vcf"));
            bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+chunkID+"_haps.vcf"));
            refSequence=chromosome.substring(refFirst,refLast+1);
            for (i=0; i<=lastDistinctHap; i++) {
                if (!distinctHaps[i].equalsIgnoreCase(refSequence)) bw.write(refChrom+"\t"+(refFirst+1)+"\t"+chunkID+"-"+(idGenerator++)+"\t"+refSequence+"\t"+distinctHaps[i]+"\t100\tPASS\t.\n");
            }
            bw.close();
            
            // Next chunk
            System.err.println("Chunk "+chunkID+" completed");
            for (i=0; i<N_SAMPLES; i++) hap1[i].delete(0,hap1[i].length());
            for (i=0; i<N_SAMPLES; i++) hap2[i].delete(0,hap2[i].length());
            str1=br1.readLine();
        }
        br1.close();
    }
    
}