import java.util.*;
import java.io.*;


/**
 * Finds ultralong interval calls that overlap calls in the main SV VCF. Prints
 * a BED file (not necessarily sorted) with format:
 * 
 * container_CHROM \t container_START \t container_END \t container_SVTYPE \t N \t D \t T \t L \t C \t sample1,...,sampleX
 * 
 * where:
 * `N` is the number of contained calls;
 * `D` is the N divided by the span of the container (a proxy for density/
 *     complexity);
 * `T` is the number of distinct SVTYPEs of contained calls;
 * `L` is the max length of a contained call;
 * `C` is a concise description of the calls in the composite event;
 * `sampleI` lists every sample with that contained set.
 * 
 * Remark: this program could be easily made much faster.
 */
public class GetOverlappingSvsUltralongAndMain {
    
    /**
     * WARNING: the program deletes every file it reads.
     * 
     * @param args 
     * 1: >=1; only contained sets with >=this number of calls in the main VCF
     *    are reported;
     * 2: only contained sets with a call >=this length are reported;
     * 4: format: `CHROM,START,END,GT1,...,GTn` (comma-separated);
     * 5: directory with one VCF per record `j` of `CONTAINER_CHUNK_FILE`, with 
     *    name `i_j.vcf` where `i=CONTAINER_CHUNK_ID`; the order of the samples
     *    must be identical between the contained and the container VCF.
     */
    public static void main(String[] args) throws IOException {
        final String CONTAINER_SVTYPE = args[0];
        final int CONTAINED_MIN_CALLS = Integer.parseInt(args[1]);
        final int CONTAINED_MIN_SV_LENGTH = Integer.parseInt(args[2]);
        final String CONTAINER_CHUNK_ID = args[3];
        final String CONTAINER_CHUNK_FILE = args[4];
        final String CONTAINED_VCFS_DIR = args[5];
        final String OUTPUT_BED = args[6];

        final int CAPACITY = 100;  // Arbitrary
        final double SCALING_FACTOR = 1.2;  // Arbitrary
        
        boolean hasDel, hasIns, hasDup, hasInv, samplesLoaded;
        int i, j;
        int lastContained, maxLength, callId, nContained, nContainedTypes, containerId, nSamples, containerLength;
        String str1, str2, keyStr, prefix;
        StringBuilder key = new StringBuilder();
        StringBuilder description = new StringBuilder();
        BufferedReader br1, br2;
        BufferedWriter bw;
        int[] contained_pos, contained_svlen;
        String[] tokens, samples, container_tokens, contained_svtype;
        byte[][] contained_gts;  // 00, 01, 10, 11
        ArrayList<String> value;
        HashMap<String,ArrayList<String>> contained2samples;
        Iterator<Map.Entry<String,ArrayList<String>>> iterator;
        Map.Entry<String,ArrayList<String>> entry;

        // Allocating reused data structures
        samples=null;
        contained_gts = new byte[CAPACITY][];
        contained_svtype = new String[CAPACITY];
        contained_svlen = new int[CAPACITY];
        contained_pos = new int[CAPACITY];
        contained2samples = new HashMap<String,ArrayList<String>>();
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_BED));
        br1 = new BufferedReader(new InputStreamReader(new FileInputStream(CONTAINER_CHUNK_FILE)));
        str1=br1.readLine(); containerId=0; samplesLoaded=false;
        while (str1!=null) {
            // Loading container and contained tokens
            container_tokens=str1.split(",");
            prefix=container_tokens[0]+"\t"+container_tokens[1]+"\t"+container_tokens[2]+"\t"+CONTAINER_SVTYPE+"\t";
            containerLength=Integer.parseInt(container_tokens[2])-Integer.parseInt(container_tokens[1]);
            nSamples=container_tokens.length-3;
            if (!samplesLoaded) samples = new String[nSamples];
            br2 = new BufferedReader(new InputStreamReader(new FileInputStream(CONTAINED_VCFS_DIR+"/"+CONTAINER_CHUNK_ID+"_"+containerId+".vcf")));
            str2=br2.readLine(); lastContained=-1;
            while (str2!=null) {
                if (str2.charAt(0)=='#') { 
                    if (!samplesLoaded && str2.startsWith("#CHROM")) {
                        tokens=str2.split("\t");
                        System.arraycopy(tokens,9,samples,0,nSamples);
                        samplesLoaded=true;
                    }
                    str2=br2.readLine(); continue;
                }
                lastContained++;
                if (lastContained==contained_gts.length) {
                    byte[][] contained_gts_new = new byte[(int)(contained_gts.length*SCALING_FACTOR)][];
                    String[] contained_svtype_new = new String[(int)(contained_gts.length*SCALING_FACTOR)];
                    int[] contained_svlen_new = new int[(int)(contained_gts.length*SCALING_FACTOR)];
                    int[] contained_pos_new = new int[(int)(contained_gts.length*SCALING_FACTOR)];
                    System.arraycopy(contained_gts,0,contained_gts_new,0,contained_gts.length);
                    System.arraycopy(contained_svtype,0,contained_svtype_new,0,contained_svtype.length);
                    System.arraycopy(contained_svlen,0,contained_svlen_new,0,contained_svlen.length);
                    System.arraycopy(contained_pos,0,contained_pos_new,0,contained_pos.length);
                    contained_gts=contained_gts_new;
                    contained_svtype=contained_svtype_new;
                    contained_svlen=contained_svlen_new;
                    contained_pos=contained_pos_new;
                }
                tokens=str2.split("\t");
                if (contained_gts[lastContained]==null) contained_gts[lastContained] = new byte[nSamples];
                for (j=0; j<nSamples; j++) contained_gts[lastContained][j]=gt2byte(tokens[9+j]);
                contained_svtype[lastContained]=getInfoField(tokens[7],"SVTYPE");
                contained_svlen[lastContained]=Math.abs(Integer.parseInt(getInfoField(tokens[7],"SVLEN")));
                contained_pos[lastContained]=Integer.parseInt(tokens[1]);
                str2=br2.readLine();
            }
            br2.close();
            new File(CONTAINED_VCFS_DIR+"/"+CONTAINER_CHUNK_ID+"_"+containerId+".vcf").delete();

            // Early exit
            if (lastContained+1<CONTAINED_MIN_CALLS) { containerId++; str1=br1.readLine(); continue; }
            maxLength=0;
            for (i=0; i<=lastContained; i++) maxLength=Math.max(maxLength,contained_svlen[i]);
            if (maxLength<CONTAINED_MIN_SV_LENGTH) { containerId++; str1=br1.readLine(); continue; }

            // Checking containment in every sample
            for (j=0; j<nSamples; j++) {
                if (!isPresent(container_tokens[3+j])) continue;
                nContained=0; maxLength=0; key.delete(0,key.length());
                for (i=0; i<=lastContained; i++) {
                    if (contained_gts[i][j]!=0) { 
                        nContained++;
                        maxLength=Math.max(maxLength,contained_svlen[i]);
                        if (key.length()==0) key.append(i+"");
                        else key.append("-"+i);
                    }
                }
                if (nContained==0 || nContained<CONTAINED_MIN_CALLS || maxLength<CONTAINED_MIN_SV_LENGTH) continue;
                keyStr=key.toString();
                if (contained2samples.containsKey(keyStr)) {
                    value=contained2samples.get(keyStr);
                    value.add(samples[j]);
                }
                else {
                    value = new ArrayList<String>();
                    value.add(samples[j]);
                    contained2samples.put(keyStr,value);
                }
            }

            // Outputting
            iterator=contained2samples.entrySet().iterator();
            while (iterator.hasNext()) {
                entry=iterator.next();
                keyStr=entry.getKey(); value=entry.getValue();
                tokens=keyStr.split("-");
                description.delete(0,description.length());
                hasDel=false; hasIns=false; hasDup=false; hasInv=false; maxLength=0;
                for (i=0; i<tokens.length; i++) {
                    callId=Integer.parseInt(tokens[i]);
                    if (contained_svtype[callId].equals("DEL")) hasDel=true;
                    if (contained_svtype[callId].equals("INS")) hasIns=true;
                    if (contained_svtype[callId].equals("DUP")) hasDup=true;
                    if (contained_svtype[callId].equals("INV")) hasInv=true;
                    if (description.length()>0) description.append(",");
                    description.append(contained_svtype[callId]).append("_").append(contained_pos[callId]).append("_").append(contained_svlen[callId]);
                    maxLength=Math.max(maxLength,contained_svlen[callId]);
                }
                nContainedTypes=(hasDel?1:0)+(hasIns?1:0)+(hasDup?1:0)+(hasInv?1:0);
                bw.write(prefix);
                bw.write(tokens.length+"\t"+(((double)tokens.length)/containerLength)+"\t"+nContainedTypes+"\t"+maxLength+"\t"+description.toString()+"\t");
                bw.write(value.get(0));
                for (i=1; i<value.size(); i++) bw.write(","+value.get(i));
                bw.newLine();
            }

            // Next iteration
            contained2samples.clear();
            containerId++;
            str1=br1.readLine();
        }
        br1.close(); bw.close();
    }


    private static final boolean isPresent(String gt) {
        if (gt.charAt(0)=='1') return true;
        if (gt.length()<3 || (gt.charAt(1)!='|' && gt.charAt(1)!='/')) return false;
        return gt.charAt(2)=='1';
    }


    private static final byte gt2byte(String gt) {
        final int left = gt.charAt(0)=='1'?1:0;
        final int right = (gt.length()<3 || (gt.charAt(1)!='|' && gt.charAt(1)!='/')) ? 0 : (gt.charAt(2)=='1'?2:0);
        return (byte)(left|right);
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