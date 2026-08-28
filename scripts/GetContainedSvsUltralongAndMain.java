import java.util.*;
import java.io.*;


/**
 * Finds ultralong interval calls that contain or overlap calls in the main SV 
 * VCF. Prints a BED file (not necessarily sorted) with format:
 * 
 * container_CHROM \t container_START \t container_END \t container_SVTYPE \t N \t T \t L \t C \t sample1,...,sampleX
 * 
 * where:
 * `N` is the number of contained calls;
 * `T` is the number of distinct SVTYPEs of contained calls;
 * `L` is the max length of a contained call;
 * `C` is a concise description of the calls in the composite event;
 * `sampleI` lists every sample with that contained set.
 * 
 * Remark: this program could be easily made much faster. It should also be made
 * to work on a chunk of container+contained VCFs, rather than on just one.
 */
public class GetContainedSvsUltralongAndMain {
    
    /**
     * @param args 
     * 0: sample order is assumed to be identical between the contained and the 
     *    container VCF;
     * 3: >=1; only contained sets with >=this number of calls in the main VCF
     *    are reported;
     * 4: only contained sets with a call >=this length are reported;
     * 6: format: `CHROM,START,END,GT1,...,GTn` (comma-separated).
     */
    public static void main(String[] args) throws IOException {
        final String CONTAINED_VCF = args[0];
        final int N_CONTAINED_RECORDS = Integer.parseInt(args[1]);
        final String CONTAINER_SVTYPE = args[2];
        final int MIN_CALLS = Integer.parseInt(args[3]);
        final int MIN_SV_LENGTH = Integer.parseInt(args[4]);
        final String OUTPUT_BED = args[5];
        final String CONTAINER_DATA = args[6];
        
        boolean hasDel, hasIns, hasDup, hasInv;
        int i, j;
        int lastContained, maxLength, callId, nContained, nContainedTypes;
        String str, keyStr;
        StringBuilder key = new StringBuilder();
        StringBuilder description = new StringBuilder();
        BufferedReader br;
        BufferedWriter bw;
        int[] contained_svlen;
        String[] tokens, samples, container_tokens, contained_svtype;
        String[][] contained_tokens;
        ArrayList<String> value;
        HashMap<String,ArrayList<String>> contained2samples;
        Iterator<Map.Entry<String,ArrayList<String>>> iterator;
        Map.Entry<String,ArrayList<String>> entry;

        if (N_CONTAINED_RECORDS<MIN_CALLS) return;

        // Loading container and contained tokens
        container_tokens = CONTAINER_DATA.split(",");
        final String PREFIX = container_tokens[0]+"\t"+container_tokens[1]+"\t"+container_tokens[2]+"\t"+CONTAINER_SVTYPE+"\t";
        final int N_SAMPLES = container_tokens.length-3;
        samples = new String[N_SAMPLES];
        contained_tokens = new String[N_CONTAINED_RECORDS][];
        contained_svtype = new String[N_CONTAINED_RECORDS];
        contained_svlen = new int[N_CONTAINED_RECORDS];
        br = new BufferedReader(new InputStreamReader(new FileInputStream(CONTAINED_VCF)));
        str=br.readLine(); lastContained=-1;
        while (str!=null) { 
            if (str.charAt(0)=='#') { 
                if (str.startsWith("#CHROM")) {
                    tokens=str.split("\t");
                    for (i=0; i<N_SAMPLES; i++) samples[i]=tokens[9+i];
                }
                str=br.readLine(); continue;
            }
            lastContained++;
            contained_tokens[lastContained]=str.split("\t");
            contained_svtype[lastContained]=getInfoField(contained_tokens[lastContained][7],"SVTYPE");
            contained_svlen[lastContained]=Integer.parseInt(getInfoField(contained_tokens[lastContained][7],"SVLEN"));
            str=br.readLine();
        }
        br.close();

        // Checking containment in every sample
        contained2samples = new HashMap<String,ArrayList<String>>();
        for (j=0; j<N_SAMPLES; j++) {
            if (!isPresent(container_tokens[3+j])) continue;
            nContained=0; key.delete(0,key.length());
            for (i=0; i<=lastContained; i++) {
                if (isPresent(contained_tokens[i][9+j])) { 
                    nContained++; 
                    if (key.length()==0) key.append(i+"");
                    else key.append("-"+i);
                }
            }
            if (nContained==0 || nContained<MIN_CALLS) continue;
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
        bw = new BufferedWriter(new FileWriter(OUTPUT_BED));
        iterator=contained2samples.entrySet().iterator();
        while (iterator.hasNext()) {
            entry=iterator.next();
            keyStr=entry.getKey(); value=entry.getValue();
            tokens=keyStr.split("-");
            maxLength=0; description.delete(0,description.length());
            hasDel=false; hasIns=false; hasDup=false; hasInv=false;
            for (i=0; i<tokens.length; i++) {
                callId=Integer.parseInt(tokens[i]);
                maxLength=Math.max(maxLength,contained_svlen[callId]);
                if (contained_svtype[callId].equals("DEL")) hasDel=true;
                if (contained_svtype[callId].equals("INS")) hasIns=true;
                if (contained_svtype[callId].equals("DUP")) hasDup=true;
                if (contained_svtype[callId].equals("INV")) hasInv=true;
                if (description.length()>0) description.append(",");
                description.append(contained_svtype[callId]).append("_").append(contained_tokens[callId][1]).append("_").append(contained_svlen[callId]);
            }
            if (maxLength>=MIN_SV_LENGTH) {
                nContainedTypes=(hasDel?1:0)+(hasIns?1:0)+(hasDup?1:0)+(hasInv?1:0);
                bw.write(PREFIX);
                bw.write(tokens.length+"\t"+nContainedTypes+"\t"+maxLength+"\t"+description.toString()+"\t");
                bw.write(value.get(0));
                for (i=1; i<value.size(); i++) bw.write(","+value.get(i));
                bw.newLine();
            }
        }
        bw.close();
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