import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * A faster, smaller-memory implementation of `truvari divide`, which assumes
 * the input VCF contains only one chromosome, is sorted, and has sequence-
 * resolved REF/ALT.
 *
 * Remark: the outputs are GZIP-compressed (not BGZIP-compressed) VCFs, so they 
 * should be converted to BGZIP by the user.
 *
 * Remark: this could be made much faster, but it is fast enough for now.
 */
public class TruvariDivide {
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String OUTPUT_DIR = args[1];
        final int BUFFER = Integer.parseInt(args[2]);
        final int MIN_RECORDS_PER_VCF = Integer.parseInt(args[3]);
        
        int i, p, q, r, s, t;
        int chunkID, pos, refLength, altLength, nRecords;
        int first, last, maxLast;  // One-based, inclusive.
        String str, headerStr;
        StringBuilder header;
        BufferedReader br;
        BufferedWriter bw;
        
        // Loading header
        header = new StringBuilder();
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)!='#') break;
            header.append(str); header.append('\n');
            str=br.readLine();
        }
        br.close();
        headerStr=header.toString();
        
        // Splitting
        chunkID=0;
        bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(OUTPUT_DIR+"/chunk_"+chunkID+".vcf.gz"))));
        bw.write(headerStr);
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine(); nRecords=0; first=0; last=0; maxLast=0;
        while (str!=null) {
            if (str.charAt(0)=='#') { str=br.readLine(); continue; }
            p=str.indexOf('\t');
            q=str.indexOf('\t',p+1);
            r=str.indexOf('\t',q+1);
            s=str.indexOf('\t',r+1);
            t=str.indexOf('\t',s+1);
            pos=Integer.parseInt(str.substring(p+1,q));
            refLength=s-r-1; altLength=t-s-1;  // Length of each VCF field
            if (refLength>1 && altLength==1) {
                // DEL
                first=pos+1;
                last=pos+refLength-1;
            }
            else if (refLength==1 && altLength>1) {
                // INS
                first=pos;
                last=pos+1;
            }
            else if (refLength>1 && altLength>1) {
                // Substitution
                first=pos+1;
                last=pos+refLength-1;
            }
            else {
                System.err.println("ERROR: Record type not recognized:");
                System.err.println(str.substring(0,t));
                System.exit(1);
            }
            if (first>maxLast+BUFFER && nRecords>=MIN_RECORDS_PER_VCF) {
                bw.close();
                System.err.println("chunk="+chunkID+" nRecords="+nRecords);
                chunkID++;
                bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(OUTPUT_DIR+"/chunk_"+chunkID+".vcf.gz"))));
                bw.write(headerStr);
                bw.write(str); bw.newLine();
                nRecords=1;
                maxLast=last;
            }
            else {
                bw.write(str); bw.newLine();
                nRecords++;
                if (last>maxLast) maxLast=last;
            }
            str=br.readLine();
        }
        br.close(); bw.close();
        System.err.println("Created "+(chunkID+1)+" truvari chunks");
    }
    
}