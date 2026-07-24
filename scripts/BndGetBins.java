import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a VCF that contains only BND records, the program prints an output BED 
 * (i.e. zero-based, left-inclusive, right-exclusive) with 4 bins per record,
 * one bin on each side of each endpoint of the BND.
 */
public class BndGetBins {
    
    private static HashMap<String,Integer> fai;
    
    /**
     * Output format: CHROM,START,END,RECORDID,BINID
     *
     * @param args
     * 2: fixed length of a breakpoint bin;
     * 3: additional length added to the right of the left-side bin (resp. to
     * the left of the right-side bin).
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String INPUT_FAI = args[1];
        final int BREAKPOINT_BIN_LENGTH = Integer.parseInt(args[2]);
        final int SLACK_BP = Integer.parseInt(args[3]);
        
        char separator;
        int p, q, first;
        int chromLength, pos;
        String str, chrom, id, alt, info;
        BufferedReader br;
        String[] tokens;
        
        loadFai(INPUT_FAI);
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            id=tokens[2];

            // First breakpoint
            chrom=tokens[0];
            chromLength=fai.get(chrom).intValue();
            pos=Integer.parseInt(tokens[1]);  // 1-based, exclusive.
            p=pos-BREAKPOINT_BIN_LENGTH;
            System.out.println(chrom+"\t"+(p>=0?p:0)+"\t"+(pos+SLACK_BP<=chromLength?pos+SLACK_BP:chromLength)+"\t"+id+"\t0");
            if (BREAKPOINT_BIN_LENGTH>0) System.out.println(chrom+"\t"+(pos-SLACK_BP>=0?pos-SLACK_BP:0)+"\t"+(pos+BREAKPOINT_BIN_LENGTH<=chromLength?pos+BREAKPOINT_BIN_LENGTH:chromLength)+"\t"+id+"\t1");

            // Second breakpoint
            alt=tokens[4];
            p=alt.indexOf('['); q=alt.indexOf(']'); first=-1; separator='_';
            if (p>=0) { separator='['; first=p; }
            else if (q>=0) { separator=']'; first=q; }
            else {
                System.err.println("ERROR: unrecognized BND ALT: "+alt);
                System.exit(1);
            }
            p=alt.indexOf(':',first+1);
            chrom=alt.substring(first+1,p);
            q=alt.indexOf(separator,p+1);
            pos=Integer.parseInt(alt.substring(p+1,q));
            chromLength=fai.get(chrom).intValue();
            p=pos-BREAKPOINT_BIN_LENGTH;
            System.out.println(chrom+"\t"+(p>=0?p:0)+"\t"+(pos+SLACK_BP<=chromLength?pos+SLACK_BP:chromLength)+"\t"+id+"\t2");
            if (BREAKPOINT_BIN_LENGTH>0) System.out.println(chrom+"\t"+(pos-SLACK_BP>=0?pos-SLACK_BP:0)+"\t"+(pos+BREAKPOINT_BIN_LENGTH<=chromLength?pos+BREAKPOINT_BIN_LENGTH:chromLength)+"\t"+id+"\t3");
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
    }
    
    
    private static final void loadFai(String path) throws IOException {
        String str;
        BufferedReader br;
        String[] tokens;
        
        fai = new HashMap<String,Integer>();
        br = new BufferedReader(new FileReader(path));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            fai.put(tokens[0],Integer.valueOf(tokens[1]));
            str=br.readLine();
        }
        br.close();
    }

}