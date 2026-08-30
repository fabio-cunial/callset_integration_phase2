import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Finds clusters of nearby SVs on the same sample (not necessarily in phase)
 * from the main SV VCF. Prints a BED file (not necessarily sorted) with format:
 * 
 * CHROM \t START \t END \t N \t T \t N_DEL \t N_INS \t N_DUP \t N_INV \t N_BND \t L \t C \t (sample1,P1);...;(sampleX,PX)
 * 
 * where:
 * `N` is the number of calls involved in the composite event;
 * `T` is the number of distinct SVTYPEs in the composite event;
 * `N_DEL` is the number of deletions in the composite event;
 * `L` is the max length of a call in the composite event;
 * `C` is a concise description of the calls in the composite event;
 * `(sampleI,PI)` lists every sample with that composite event, where `PI=1` iff 
 * >=minCalls in sampleI lie on the same haplotype.
 * 
 * Remark: the program can be used on the ultralong VCF, on the BND VCF, and
 * even on the bcftools concat of ultralong and BND VCF.
 * 
 * Remark: this program could be easily made much faster.
 */
public class GetCompositeSvsPrime {

    private static final int CAPACITY = 20;  // Arbitrary

    /**
     * Sample names from the VCF header
     */
    private static String[] samples;
    private static int nSamples;

    /**
     * Temporary space to store all calls that form a cluster in the multi-
     * sample VCF.
     */
    private static String calls_chrom;
    private static String[][] calls_genotypes = new String[CAPACITY][]; 
    private static int[][] calls_start_end = new int[CAPACITY][2];  // Zero-based, inclusive.
    private static int[] calls_svlen = new int[CAPACITY];
    private static String[] calls_svtype = new String[CAPACITY];
    private static int[] calls_pos = new int[CAPACITY];
    private static int lastCall;

    /**
     * Temporary space to store the calls from the multi-sample cluster that
     * occur on the same sample.
     */
    private static int[] sample = new int[CAPACITY];
    private static String[] sample_gt = new String[CAPACITY];
    private static int sampleLast;

    /**
     * Temporary, reused space.
     */
    private static StringBuilder path = new StringBuilder();
    private static HashMap<String,ArrayList<String>> compositeSv2samples = new HashMap<String,ArrayList<String>>();
    private static int[][] compositeSvs = new int[CAPACITY][3];
    private static int lastCompositeSv;
    
    
    /**
     * @param args 
     * 0: a sorted inter-sample VCF;
     * 1: treat intra-chromosomal BNDs as intervals (0) or as isolated 
     *    breakpoints (1); in the first (resp. second) case, the VCF is assumed
     *    to contain only one BND record (resp. two symmetrical BND records) per
     *    event; inter-chromosomal BNDs are always assumed to be represented as
     *    pairs of symmetrical records;
     * 4: only clusters with at least one call >=this length are considered.
     */
    public static void main(String[] args) throws IOException {
        final String COHORT_VCF_GZ = args[0];
        final int BND_MODE = Integer.parseInt(args[1]);
        final int MAX_DISTANCE = Integer.parseInt(args[2]);
        final int MIN_CALLS = Integer.parseInt(args[3]);
        final int MIN_SV_LENGTH = Integer.parseInt(args[4]);
        final String OUTPUT_BED = args[5];
        
        final int QUANTUM = 10000;  // Arbitrary
        
        int i;
        int nRecords, currentFirst, currentLast;
        String str, chrom, currentChrom;
        BufferedReader br;
        BufferedWriter bw;
        int[] tmpArray = new int[2];
        String[] tokens;
        
        lastCall=-1;
        samples=null;
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(COHORT_VCF_GZ))));
        bw = new BufferedWriter(new FileWriter(OUTPUT_BED));
        str=br.readLine(); currentChrom=null; currentFirst=-1; currentLast=-1; nRecords=0; nSamples=0;
        while (str!=null) { 
            if (str.charAt(0)=='#') { 
                if (str.substring(0,6).equalsIgnoreCase("#CHROM")) {
                    tokens=str.split("\t");
                    nSamples=tokens.length-9;
                    samples = new String[nSamples];
                    System.arraycopy(tokens,9,samples,0,nSamples);
                }
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records...");
            tokens=str.split("\t");
            chrom=tokens[0];
            getInterval(tokens,BND_MODE,tmpArray);
            if (currentFirst==-1) { 
                currentChrom=chrom; currentFirst=tmpArray[0]; currentLast=tmpArray[1];
                for (i=0; i<=lastCall; i++) calls_genotypes[i]=null;
                calls_genotypes[0] = new String[nSamples];
                System.arraycopy(tokens,9,calls_genotypes[0],0,nSamples);
                calls_start_end[0][0]=tmpArray[0];
                calls_start_end[0][1]=tmpArray[1];
                calls_svtype[0]=getInfoField(tokens[7],"SVTYPE");
                if (calls_svtype[0].equals("BND")) calls_svlen[0]=BND_MODE==0?calls_start_end[0][1]-calls_start_end[0][0]:Integer.MAX_VALUE;
                else calls_svlen[0]=Math.abs(Integer.parseInt(getInfoField(tokens[7],"SVLEN")));
                calls_pos[0]=Integer.parseInt(tokens[1]);
                calls_chrom=chrom;
                lastCall=0;
            }
            else if (!chrom.equals(currentChrom) || tmpArray[0]>currentLast+MAX_DISTANCE) {
                getCompositeSvs(MAX_DISTANCE,MIN_CALLS,MIN_SV_LENGTH,bw);
                currentChrom=chrom; currentFirst=tmpArray[0]; currentLast=tmpArray[1];
                for (i=0; i<=lastCall; i++) calls_genotypes[i]=null;
                calls_genotypes[0] = new String[nSamples];
                System.arraycopy(tokens,9,calls_genotypes[0],0,nSamples);
                calls_start_end[0][0]=tmpArray[0];
                calls_start_end[0][1]=tmpArray[1];
                calls_svtype[0]=getInfoField(tokens[7],"SVTYPE");
                if (calls_svtype[0].equals("BND")) calls_svlen[0]=BND_MODE==0?calls_start_end[0][1]-calls_start_end[0][0]:Integer.MAX_VALUE;
                else calls_svlen[0]=Math.abs(Integer.parseInt(getInfoField(tokens[7],"SVLEN")));
                calls_pos[0]=Integer.parseInt(tokens[1]);
                calls_chrom=chrom;
                lastCall=0;
            }
            else {
                if (tmpArray[1]>currentLast) currentLast=tmpArray[1];
                lastCall++;
                if (lastCall==calls_genotypes.length) {
                    String[][] newCalls = new String[calls_genotypes.length<<1][];
                    System.arraycopy(calls_genotypes,0,newCalls,0,calls_genotypes.length);
                    calls_genotypes=newCalls;
                    int[][] newCallsStartEnd = new int[calls_start_end.length<<1][2];
                    System.arraycopy(calls_start_end,0,newCallsStartEnd,0,calls_start_end.length);
                    calls_start_end=newCallsStartEnd;
                    int[] newCallsSvlen = new int[calls_svlen.length<<1];
                    System.arraycopy(calls_svlen,0,newCallsSvlen,0,calls_svlen.length);
                    calls_svlen=newCallsSvlen;
                    String[] newCallsSvtype = new String[calls_svtype.length<<1];
                    System.arraycopy(calls_svtype,0,newCallsSvtype,0,calls_svtype.length);
                    calls_svtype=newCallsSvtype;
                    int[] newCallsPos = new int[calls_pos.length<<1];
                    System.arraycopy(calls_pos,0,newCallsPos,0,calls_pos.length);
                    calls_pos=newCallsPos;
                }
                calls_genotypes[lastCall] = new String[nSamples];
                System.arraycopy(tokens,9,calls_genotypes[lastCall],0,nSamples);
                calls_start_end[lastCall][0]=tmpArray[0];
                calls_start_end[lastCall][1]=tmpArray[1];
                calls_svtype[lastCall]=getInfoField(tokens[7],"SVTYPE");
                if (calls_svtype[lastCall].equals("BND")) calls_svlen[lastCall]=BND_MODE==0?calls_start_end[lastCall][1]-calls_start_end[lastCall][0]:Integer.MAX_VALUE;
                else calls_svlen[lastCall]=Math.abs(Integer.parseInt(getInfoField(tokens[7],"SVLEN")));
                calls_pos[lastCall]=Integer.parseInt(tokens[1]);
            }
            str=br.readLine();
        }
        getCompositeSvs(MAX_DISTANCE,MIN_CALLS,MIN_SV_LENGTH,bw);
        br.close(); bw.close();
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
    
    
    /**
     * @param maxDistance maximum distance between two calls that are present on
     * the same sample to be considered part of the same composite event;
     * @param minCalls (>=2) min number of calls in a composite event;
     * @param minSvLength only clusters with at least one call of this length or
     * longer are printed in output.
     */
    private static final void getCompositeSvs(int maxDistance, int minCalls, int minSvLength, BufferedWriter outputBed) throws IOException {
        boolean found;
        int i, j, k;
        int callId, start, end, nDel, nIns, nDup, nInv, nBnd, nTypes, maxLength;
        String key;
        ArrayList<String> value;
        Iterator<Map.Entry<String,ArrayList<String>>> iterator;
        Map.Entry<String,ArrayList<String>> entry;
        String[] tokens;

        // Early exit in trivial cases
        if (lastCall+1<minCalls) return;
        found=false;
        for (i=0; i<=lastCall; i++) {
            if (calls_svlen[i]>=minSvLength) { found=true; break; }
        }
        if (!found) return;
        
        // Building `compositeSvs` and `compositeSv2samples`.
        if (compositeSvs.length<lastCall+1) compositeSvs = new int[lastCall+1][3];
        compositeSv2samples.clear();
        if (sample.length<lastCall+1) {
            sample = new int[lastCall+1];
            sample_gt = new String[lastCall+1];
        }
        for (j=0; j<nSamples; j++) {
            loadSample(j);
            if (sampleLast==-1) continue;
            getCompositeSvsInSample(maxDistance,minCalls);
            for (k=0; k<=lastCompositeSv; k++) {
                key=compositeSv2key(compositeSvs[k][0],compositeSvs[k][1]);
                if (compositeSv2samples.containsKey(key)) {
                    value=compositeSv2samples.get(key);
                    value.add(samples[j]);
                    value.add(compositeSvs[k][2]+"");
                }
                else {
                    value = new ArrayList<String>();
                    value.add(samples[j]);
                    value.add(compositeSvs[k][2]+"");
                    compositeSv2samples.put(key,value);
                }
            }
        }
        
        // Printing the output BED
        iterator=compositeSv2samples.entrySet().iterator();
        while (iterator.hasNext()) {
            entry=iterator.next();
            key=entry.getKey(); value=entry.getValue();
            tokens=key.split("-"); found=false;
            start=Integer.MAX_VALUE; end=0; path.delete(0,path.length()); nDel=0; nIns=0; nDup=0; nInv=0; nBnd=0; maxLength=0;
            for (i=0; i<tokens.length; i++) {
                callId=Integer.parseInt(tokens[i]);
                if (calls_svlen[callId]>=minSvLength) found=true;
                if (calls_svlen[callId]!=Integer.MAX_VALUE) maxLength=Math.max(maxLength,calls_svlen[callId]);
                start=Math.min(start,calls_start_end[callId][0]);
                end=Math.max(end,calls_start_end[callId][1]);
                if (calls_svtype[callId].equals("DEL")) nDel++;
                else if (calls_svtype[callId].equals("INS")) nIns++;
                else if (calls_svtype[callId].equals("DUP")) nDup++;
                else if (calls_svtype[callId].equals("INV")) nInv++;
                else if (calls_svtype[callId].equals("BND")) nBnd++;
                path.append(path.length()==0?"":",").append(calls_svtype[callId]).append("_").append(calls_pos[callId]).append("_").append(calls_svlen[callId]!=Integer.MAX_VALUE?calls_svlen[callId]:"NA");
            }
            if (found) {
                nTypes=(nDel>0?1:0)+(nIns>0?1:0)+(nDup>0?1:0)+(nInv>0?1:0)+(nBnd>0?1:0);
                outputBed.write(calls_chrom+"\t"+start+"\t"+end+"\t"+tokens.length+"\t"+nTypes+"\t"+nDel+"\t"+nIns+"\t"+nDup+"\t"+nInv+"\t"+nBnd+"\t"+maxLength+"\t"+path.toString()+"\t");
                outputBed.write("("+value.get(0)+","+value.get(1)+")");
                for (i=2; i<value.size(); i+=2) { outputBed.write(";("+value.get(i)+","+value.get(i+1)+")"); }
                outputBed.newLine();
            }
        }
    }
    
    
    /**
     * The procedure loads call IDs (i.e. indexes in `calls_*`) in global
     * variable `sample`.
     */
    private static final void loadSample(int column) {
        int i, p;

        sampleLast=-1;
        for (i=0; i<=lastCall; i++) {
            if (calls_genotypes[i][column].charAt(0)=='1' || (calls_genotypes[i][column].length()>=3 && (calls_genotypes[i][column].charAt(1)=='/' || calls_genotypes[i][column].charAt(1)=='|') && calls_genotypes[i][column].charAt(2)=='1')) {
                sampleLast++;
                sample[sampleLast]=i;
                p=calls_genotypes[i][column].indexOf(':');
                if (p<0) p=calls_genotypes[i][column].length();
                sample_gt[sampleLast]=calls_genotypes[i][column].substring(0,p);
            }
        }
    }
    
    
    /**
     * Scans `sample` for blocks of >=`minCalls` consecutive calls at distance 
     * <=`maxDistance` from each other, and stores `(i,j,k)` pairs in 
     * `compositeSvs`, where `i` (resp. `j`) is the index of the first (resp. 
     * last) call in a block, and `k` is 1 iff `>=minCalls` in that block lie on
     * the same haplotype.
     */
    private static final void getCompositeSvsInSample(int maxDistance, int minCalls) {
        int i;
        int currentFirst, currentLast, currentEnd, nCallsLeft, nCallsRight;
        
        lastCompositeSv=-1;
        if (isPhased(sample_gt[0])) {
            nCallsLeft=sample_gt[0].charAt(0)=='1'?1:0;
            nCallsRight=(sample_gt[0].length()==3 && sample_gt[0].charAt(2)=='1')?1:0;
        }
        else { nCallsLeft=0; nCallsRight=0; }
        currentFirst=0; currentLast=0; currentEnd=calls_start_end[sample[0]][1];
        for (i=1; i<=sampleLast; i++) {
            if (calls_start_end[sample[i]][0]-currentEnd>maxDistance) {
                if (currentLast-currentFirst+1>=minCalls) {
                    lastCompositeSv++;
                    compositeSvs[lastCompositeSv][0]=currentFirst;
                    compositeSvs[lastCompositeSv][1]=currentLast;
                    compositeSvs[lastCompositeSv][2]=(nCallsLeft>=minCalls||nCallsRight>=minCalls)?1:0;
                }
                currentFirst=i; nCallsLeft=0; nCallsRight=0;
            }
            currentLast=i; currentEnd=Math.max(currentEnd,calls_start_end[sample[i]][1]);
            if (isPhased(sample_gt[i])) {
                nCallsLeft+=sample_gt[i].charAt(0)=='1'?1:0;
                nCallsRight+=(sample_gt[i].length()==3 && sample_gt[i].charAt(2)=='1')?1:0;
            }
        }
        if (currentLast-currentFirst+1>=minCalls) {
            lastCompositeSv++;
            compositeSvs[lastCompositeSv][0]=currentFirst;
            compositeSvs[lastCompositeSv][1]=currentLast;
            compositeSvs[lastCompositeSv][2]=(nCallsLeft>=minCalls||nCallsRight>=minCalls)?1:0;
        }
    }


    private static final boolean isPhased(String gt) {
        return gt.length()==1 || (
               gt.length()==3 && (
                 gt.charAt(1)=='|' || 
                 (gt.charAt(0)=='1' && gt.charAt(2)=='1')
               )
             );
    }
    
    
    /**
     * @param out [start,end], 0-based, inclusive.
     */
    private static final void getInterval(String[] call, int bndMode, int[] out) {        
        int from, to, pos, svlen, altPos;
        String svtype, altChr;

        pos=Integer.parseInt(call[1]);
        svtype=getInfoField(call[7],"SVTYPE");
        if (svtype.equals("INS")) {
            from=pos-1;  // Zero-based, inclusive.
            to=pos;  // Zero-based, inclusive.
        }
        else if (svtype.equals("BND")) {
            altChr=bndGetAltChr(call[4]);
            if (bndMode==0 && altChr.equals(call[0])) {
                altPos=bndGetAltPos(call[4]);
                if (altPos<pos) {
                    System.err.println("ERROR: BND with altPos<POS: "+altPos+" < "+pos+". Has the BND VCF been canonized?");
                    System.exit(1);
                }
                from=pos-1;  // Zero-based, inclusive.
                to=altPos-1;  // Zero-based, inclusive.
            }
            else {
                from=pos-1;  // Zero-based, inclusive.
                to=pos+1;  // Zero-based, inclusive.
            }
        }
        else {
            from=pos;  // Zero-based, inclusive.
            svlen=Math.abs(Integer.parseInt(getInfoField(call[7],"SVLEN")));
            to=pos+svlen-1;  // Zero-based, inclusive.
        }
        out[0]=from; out[1]=to;
    }


    private static String bndGetAltChr(String alt) {
        int p, q, r;

        p=alt.indexOf('['); q=alt.indexOf(']'); 
        if (p>=0) {
            r=alt.indexOf(':',p+1);
            return alt.substring(p+1,r);
        }
        else if (q>=0) {
            r=alt.indexOf(':',q+1);
            return alt.substring(q+1,r);
        }
        else {
            System.err.println("ERROR: unrecognized ALT: "+alt);
            System.exit(1);
        }
        return null;
    }


    private static int bndGetAltPos(String alt) {
        int p, q, r;

        p=alt.indexOf('['); q=alt.indexOf(']'); 
        if (p>=0) {
            r=alt.indexOf(':',p+1);
            p=alt.indexOf('[',r+1);
            return Integer.parseInt(alt.substring(r+1,p));
        }
        else if (q>=0) {
            r=alt.indexOf(':',q+1);
            q=alt.indexOf(']',r+1);
            return Integer.parseInt(alt.substring(r+1,q));
        }
        else {
            System.err.println("ERROR: unrecognized ALT: "+alt);
            System.exit(1);
        }
        return -1;
    }
    
    
    /**
     * @param first,last positions in `sample`.
     */
    private static final String compositeSv2key(int first, int last) {
        int i;
        String out;
        
        out=sample[first]+"";
        for (i=first+1; i<=last; i++) out+="-"+sample[i];
        return out;
    }
    
}