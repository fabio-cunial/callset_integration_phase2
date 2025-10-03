import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Uses the NE field printed by kanpig to build a BED of kanpig regions.
 */
public class GetComplexSvs {
    
    /**
     * 
     *
     * @param args 0 a VCF whose FORMAT field is:
     *
     * GT:FT:SQ:GQ:PS:NE:DP:AD:KS
     *
     */
    public static void main(String[] args) throws IOException {
        final String COHORT_VCF_GZ = args[0];
        
        final int QUANTUM = 10000;
        
        int lastCall, nRecords, region, currentRegion;
        String str;
        BufferedReader br;
        String[] tokens, tokensPrime;
        String[][] calls;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(COHORT_VCF_GZ))));
        str=br.readLine(); currentRegion=-1; nRecords=0;
        while (str!=null) { 
            if (str.charAt(0)=='#') { str=br.readLine(); continue; }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records...");
            tokens=str.split("\t");
            tokensPrime=tokens[9].split(":");
            region=Integer.parseInt(tokensPrime[5]);
            if (region!=currentRegion) {
                if (currentRegion!=-1) getComplexSvs(calls,lastCall);
                currentRegion=region; lastCall=0; calls[0]=tokens;
            }
            else {
                lastCall++;
                if (lastCall==calls.length) {
                    String[][] newCalls = new String[calls.length<<1][calls[0].length];
                    System.arraycopy(calls,0,newCalls,0,calls.length);
                }
                calls[lastCall]=tokens;
            }
            str=br.readLine();
        }
        if (currentRegion!=-1) getComplexSvs(calls,lastCall);
        br.close();
    }
    
    
    /**
     * @param maxDistance maximum distance between two calls to be considered
     * part of the same complex event;
     * @param minCalls (>=2) min number of calls in a complex event.
     */
    private static final void getComplexSvs(String[][] calls, int lastCall, int maxDistance, int minCalls) {
        final int N_COLUMNS = calls[0].length;
        final boolean IS_AUTOSOMAL = calls[0][9].indexOf(":")==3;
        
        final int CAPACITY = 20;  // Arbitrary
        
        int i, j;
        int last;
        int[] haplotype, tmpArray;
        
        haplotype = new int[CAPACITY];
        tmpArray = new int[2];
        for (j=10; j<N_COLUMNS; j++) {
            loadHaplotype(calls,lastCall,false,IS_AUTOSOMAL);
            getComplexSVs(calls,maxDistance,minCalls,tmpArray);
            if (IS_AUTOSOMAL) {
                loadHaplotype(calls,lastCall,true,true);
                getComplexSVs(calls,maxDistance,minCalls,tmpArray);
            }
            
            
        }
        
        
    }
    
    
    private static int[] haplotype;
    private static int haplotypeLast;
    
    
    /**
     * Remark: the procedure loads call IDs in global variable `haplotype`.
     *
     * @param which 0=left, 1=right.
     */
    private static final void loadHaplotype(String[][] calls, int lastCall, boolean which, boolean autosomal) {
        int i;
        
        haplotypeLast=-1;
        for (i=0; i<=lastCall; i++) {
            if (autosomal && calls[i][j].charAt(1)!='|') {
                // Skipping unphased calls
                continue;
            }
            if (calls[i][j].charAt(which?2:0)=='1') {
                haplotypeLast++;
                if (haplotypeLast==haplotype.length) {
                    int[] newArray = new int[haplotype.length<<1];
                    System.arraycopy(haplotype,0,newArray,0,haplotype.length);
                    haplotype=newArray;
                }
                haplotype[haplotypeLast]=i;
            }
        }
    }
    
    
    private static final void getComplexSVs(String[][] calls, int maxDistance, int minCalls, int[] tmpArray) {
        
        
        getInterval(calls,haplotype[0],tmpArray);
        currentFirst=0; currentLast=0; currentStart=tmpArray[0]; currentEnd=tmpArray[1];
        for (i=1; i<=haplotypeLast; i++) {
            getInterval(calls,haplotype[i],tmpArray);
            if (tmpArray[0]-currentEnd>maxDistance) {
                if (currentLast-currentFirst+1>=minCalls) {
                    // Print complex SV...
                    
                    // Increment its count to get AF in the end...
                    
                    // Next iteration
                    currentFirst=i; currentStart=tmpArray[0];
                }
                currentLast=i; currentEnd=tmpArray[1];
            }
            
        }
        
    }
    
    
    /**
     * @param out [start,end], 1-based, inclusive.
     */
    private static final void getInterval(String[][] calls, int id, int[] out) {
        final int POS = Integer.parseInt(calls[id][1]);
        final String REF = calls[id][3];
        final String ALT = calls[id][4];
        
        int start, end;
        if (REF.length()==1 && ALT.length()>1) {
            // INS
            start=pos;
            end=pos+1;
        }
        else if (REF.length()>1 && ALT.length()==1) {
            // DEL
            start=pos+1;
            end=pos+ref.length();
        }
        else {
            start=pos+1;
            end=pos+ref.length();
        }
        out[0]=start; out[1]=end;
    }
    
    
}