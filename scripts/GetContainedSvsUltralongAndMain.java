import java.util.*;
import java.io.*;


/**
 * Finds ultralong interval calls that contain or overlap calls in the main SV 
 * VCF. Prints a BED file (not necessarily sorted) with format:
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
public class GetContainedSvsUltralongAndMain {
    
    /**
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
        
        boolean hasDel, hasIns, hasDup, hasInv, samplesLoaded;
        int i, j;
        int lastContained, maxLength, callId, nContained, nContainedTypes, containerId, nSamples, containerLength;
        String str1, str2, keyStr, prefix;
        StringBuilder key = new StringBuilder();
        StringBuilder description = new StringBuilder();
        BufferedReader br1, br2;
        BufferedWriter bw;
        int[] contained_svlen;
        String[] tokens, samples, container_tokens, contained_svtype;
        String[][] contained_tokens;
        ArrayList<String> value;
        HashMap<String,ArrayList<String>> contained2samples;
        Iterator<Map.Entry<String,ArrayList<String>>> iterator;
        Map.Entry<String,ArrayList<String>> entry;

        // Allocating reused data structures
        samples=null;
        contained_tokens = new String[CAPACITY][];
        contained_svtype = new String[CAPACITY];
        contained_svlen = new int[CAPACITY];
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
                if (lastContained==contained_tokens.length) {
                    String[][] contained_tokens_new = new String[contained_tokens.length*2][];
                    String[] contained_svtype_new = new String[contained_tokens.length*2];
                    int[] contained_svlen_new = new int[contained_tokens.length*2];
                    System.arraycopy(contained_tokens,0,contained_tokens_new,0,contained_tokens.length);
                    System.arraycopy(contained_svtype,0,contained_svtype_new,0,contained_svtype.length);
                    System.arraycopy(contained_svlen,0,contained_svlen_new,0,contained_svlen.length);
                    contained_tokens=contained_tokens_new;
                    contained_svtype=contained_svtype_new;
                    contained_svlen=contained_svlen_new;
                }
                contained_tokens[lastContained]=str2.split("\t");
                contained_svtype[lastContained]=getInfoField(contained_tokens[lastContained][7],"SVTYPE");
                contained_svlen[lastContained]=Math.abs(Integer.parseInt(getInfoField(contained_tokens[lastContained][7],"SVLEN")));
                str2=br2.readLine();
            }
            br2.close();

            // Early exit
            if (lastContained+1<CONTAINED_MIN_CALLS) { containerId++; str1=br1.readLine(); continue; }
            maxLength=0;
            for (i=0; i<=lastContained; i++) maxLength=Math.max(maxLength,contained_svlen[i]);
            if (maxLength<CONTAINED_MIN_SV_LENGTH) { containerId++; str1=br1.readLine(); continue; }

            // Checking containment in every sample
            contained2samples.clear();
            for (j=0; j<nSamples; j++) {
                if (!isPresent(container_tokens[3+j])) continue;
                nContained=0; maxLength=0; key.delete(0,key.length());
                for (i=0; i<=lastContained; i++) {
                    if (isPresent(contained_tokens[i][9+j])) { 
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
                    description.append(contained_svtype[callId]).append("_").append(contained_tokens[callId][1]).append("_").append(contained_svlen[callId]);
                    maxLength=Math.max(maxLength,contained_svlen[callId]);
                }
                nContainedTypes=(hasDel?1:0)+(hasIns?1:0)+(hasDup?1:0)+(hasInv?1:0);
                bw.write(prefix);
                bw.write(tokens.length+"\t"+String.format("%.3g",((double)tokens.length)/containerLength)+"\t"+nContainedTypes+"\t"+maxLength+"\t"+description.toString()+"\t");
                bw.write(value.get(0));
                for (i=1; i<value.size(); i++) bw.write(","+value.get(i));
                bw.newLine();
            }

            // Next iteration
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