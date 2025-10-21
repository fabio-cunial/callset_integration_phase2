import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Finds clusters of nearby SVs that are phased on the same haplotype.
 */
public class GetCompositeSvs {
    /**
     * Temporary space to store the calls that occur consecutively on the same 
     * haplotype.
     */
    private static int[] haplotype;
    private static int haplotypeLast;
    
    
    /**
     * @param args 
     * 0: a sorted and phased inter-sample VCF. Only phased calls are
     *    considered;
     * 3: only clusters with at least one call of this length or longer are
     *    considered.
     */
    public static void main(String[] args) throws IOException {
        final String COHORT_VCF_GZ = args[0];
        final int MAX_DISTANCE = Integer.parseInt(args[1]);
        final int MIN_CALLS = Integer.parseInt(args[2]);
        final int MIN_SV_LENGTH = Integer.parseInt(args[3]);
        
        final int CAPACITY = 10;  // Arbitrary
        final int QUANTUM = 10000;  // Arbitrary
        
        int lastCall, nRecords, region, currentRegion;
        String str;
        BufferedReader br;
        String[] tokens, tokensPrime, columns;
        String[][] calls;
        
        calls = new String[CAPACITY][0]; lastCall=-1;
        columns=null;
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(COHORT_VCF_GZ))));
        str=br.readLine(); currentRegion=-1; nRecords=0;
        while (str!=null) { 
            if (str.charAt(0)=='#') { 
                if (str.substring(0,6).equalsIgnoreCase("#CHROM")) columns=str.split("\t");
                str=br.readLine();
                continue; 
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records...");
            tokens=str.split("\t");
            tokensPrime=tokens[9].split(":");
            region=Integer.parseInt(tokensPrime[5]);
            if (region!=currentRegion) {
                if (currentRegion!=-1) getCompositeSvs(calls,lastCall,MAX_DISTANCE,MIN_CALLS,MIN_SV_LENGTH,columns);
                currentRegion=region; lastCall=0; calls[0]=tokens;
            }
            else {
                lastCall++;
                if (lastCall==calls.length) {
                    String[][] newCalls = new String[calls.length<<1][calls[0].length];
                    System.arraycopy(calls,0,newCalls,0,calls.length);
                    calls=newCalls;
                }
                calls[lastCall]=tokens;
            }
            str=br.readLine();
        }
        if (currentRegion!=-1) getCompositeSvs(calls,lastCall,MAX_DISTANCE,MIN_CALLS,MIN_SV_LENGTH,columns);
        br.close();
    }
    
    
    /**
     * @param maxDistance maximum distance between two calls that are present on
     * the same haplotype to be considered part of the same composite event;
     * @param minCalls (>=2) min number of calls in a composite event;
     * @param minSvLength only clusters with at least one call of this length or
     * longer are printed in output.
     */
    private static final void getCompositeSvs(String[][] calls, int lastCall, int maxDistance, int minCalls, int minSvLength, String[] columns) throws IOException {
        final int N_COLUMNS = calls[0].length;
        final boolean IS_AUTOSOMAL = calls[0][9].indexOf(":")==3;
        
        final int CAPACITY = 20;  // Arbitrary
        
        boolean found;
        int i, j, k;
        int lastCompositeSv, callId, idGenerator, refLength, altLength;
        String key;
        BufferedWriter bw1, bw2, bw3;
        int[] tmpArray;
        String[] tokens;
        int[][] compositeSvs;
        Vector<String> value;
        HashMap<String,Vector<String>> compositeSv2samples;
        Iterator<Map.Entry<String,Vector<String>>> iterator;
        Map.Entry<String,Vector<String>> entry;
        
        // Building `compositeSvs` and `compositeSv2samples`.
        compositeSvs = new int[lastCall+1][2];
        compositeSv2samples = new HashMap<String,Vector<String>>();
        haplotype = new int[CAPACITY];
        tmpArray = new int[2];
        for (j=10; j<N_COLUMNS; j++) {
            loadHaplotype(j,false,IS_AUTOSOMAL,calls,lastCall);
            lastCompositeSv=getCompositeSvs(calls,maxDistance,minCalls,compositeSvs,tmpArray);
            found=false;
            for (k=0; k<=lastCompositeSv; k++) {
                key=compositeSv2key(compositeSvs[k][0],compositeSvs[k][1]);
                if (compositeSv2samples.containsKey(key)) {
                    value=compositeSv2samples.get(key);
                    value.add(columns[j]+"_0");
                }
                else {
                    value = new Vector<String>();
                    value.add(columns[j]+"_0");
                    compositeSv2samples.put(key,value);
                }
            }
            if (IS_AUTOSOMAL) {
                loadHaplotype(j,true,true,calls,lastCall);
                lastCompositeSv=getCompositeSvs(calls,maxDistance,minCalls,compositeSvs,tmpArray);
                for (k=0; k<=lastCompositeSv; k++) {
                    key=compositeSv2key(compositeSvs[k][0],compositeSvs[k][1]);
                    if (compositeSv2samples.containsKey(key)) {
                        value=compositeSv2samples.get(key);
                        value.add(columns[j]+"_1");
                    }
                    else {
                        value = new Vector<String>();
                        value.add(columns[j]+"_1");
                        compositeSv2samples.put(key,value);
                    }
                }
            }
        }
        
        // Printing `compositeSvs` and `compositeSv2samples`.
        idGenerator=0;
        iterator=compositeSv2samples.entrySet().iterator();
        while (iterator.hasNext()) {
            entry=iterator.next();
            key=entry.getKey(); value=entry.getValue();
            tokens=key.split("-"); found=false;
            for (i=0; i<tokens.length; i++) {
                callId=Integer.parseInt(tokens[i]);
                refLength=calls[callId][3].length();
                altLength=calls[callId][4].length();
                if (refLength-altLength>=minSvLength || altLength-refLength>=minSvLength) { found=true; break; }
            }
            if (!found) continue;
            idGenerator++;
            bw1 = new BufferedWriter(new FileWriter(idGenerator+"_graph.txt"));
            bw1.write("-> ");
            bw2 = new BufferedWriter(new FileWriter(idGenerator+"_refalt.txt"));
            for (i=0; i<tokens.length; i++) {
                callId=Integer.parseInt(tokens[i]);
                refLength=calls[callId][3].length();
                altLength=calls[callId][4].length();
                if (refLength==1 && altLength>1) {
                    bw1.write("INS_"+(altLength-refLength)+" -> ");
                    bw2.write(calls[callId][4]+"\n");
                }
                else if (refLength>1 && altLength==1) {
                    bw1.write("DEL_"+(refLength-altLength)+" -> ");
                    bw2.write(calls[callId][3]+"\n");
                }
                else {
                    bw1.write("REPL_"+refLength+" -> ");
                    bw2.write(calls[callId][3]+" -> "+calls[callId][4]+"\n");
                }
            }
            bw1.close(); bw2.close();
            bw3 = new BufferedWriter(new FileWriter(idGenerator+"_haps.txt"));
            for (i=0; i<value.size(); i++) bw3.write(value.elementAt(i)+"\n");
            bw3.close();
        }
    }
    
    
    /**
     * The procedure loads call IDs (i.e. indexes in `calls`) in global variable
     * `haplotype`. Only phased calls are loaded.
     *
     * @param which 0=left hap, 1=right hap.
     */
    private static final void loadHaplotype(int column, boolean which, boolean isAutosomal, String[][] calls, int lastCall) {
        int i, j;
        
        haplotypeLast=-1;
        for (i=0; i<=lastCall; i++) {
            if (isAutosomal && calls[i][column].charAt(1)!='|') {
                // Skipping unphased calls
                continue;
            }
            if (calls[i][column].charAt(which?2:0)=='1') {
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
    
    
    /**
     * Scans `haplotype` for blocks of >=`minCalls` consecutive calls at
     * distance <=`maxDistance` from each other, and stores `(i,j)` pairs in 
     * `out`, where `i` (resp. `j`) is the index of the first (resp. last) call 
     * in a block (indexes are WRT `haplotype`).
     *
     * @param tmpArray temporary space with >=2 cells;
     * @return the last element in `out`.
     */
    private static final int getCompositeSvs(String[][] calls, int maxDistance, int minCalls, int[][] out, int[] tmpArray) {
        int i;
        int lastCompositeEvent, currentFirst, currentLast, currentStart, currentEnd;
        
        lastCompositeEvent=-1;
        getInterval(haplotype[0],calls,tmpArray);
        currentFirst=0; currentLast=0; currentStart=tmpArray[0]; currentEnd=tmpArray[1];
        for (i=1; i<=haplotypeLast; i++) {
            getInterval(haplotype[i],calls,tmpArray);
            if (tmpArray[0]-currentEnd>maxDistance) {
                if (currentLast-currentFirst+1>=minCalls) {
                    lastCompositeEvent++;
                    out[lastCompositeEvent][0]=currentFirst;
                    out[lastCompositeEvent][1]=currentLast;
                }
                currentFirst=i; currentStart=tmpArray[0];
            }
            currentLast=i; currentEnd=tmpArray[1];
        }
        if (currentLast-currentFirst+1>=minCalls) {
            lastCompositeEvent++;
            out[lastCompositeEvent][0]=currentFirst;
            out[lastCompositeEvent][1]=currentLast;
        }
        return lastCompositeEvent;
    }
    
    
    /**
     * @param id position in `calls`;
     * @param out [start,end], 1-based, inclusive.
     */
    private static final void getInterval(int id, String[][] calls, int[] out) {
        final int POS = Integer.parseInt(calls[id][1]);
        final String REF = calls[id][3];
        final String ALT = calls[id][4];
        
        int start, end;
        
        if (REF.length()==1 && ALT.length()>1) {
            // INS
            start=POS;
            end=POS+1;
        }
        else if (REF.length()>1 && ALT.length()==1) {
            // DEL
            start=POS+1;
            end=POS+REF.length();
        }
        else {
            start=POS+1;
            end=POS+REF.length();
        }
        out[0]=start; out[1]=end;
    }
    
    
    /**
     * @para, first,last positions in `haplotype`.
     */
    private static final String compositeSv2key(int first, int last) {
        int i;
        String out;
        
        out=haplotype[first]+"";
        for (i=first+1; i<=last; i++) out+="-"+haplotype[i];
        return out;
    }
    
}