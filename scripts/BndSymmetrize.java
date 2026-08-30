import java.io.*;
import java.util.zip.GZIPInputStream;


/**
 * Given a BND-only VCF, the program adds a symmetrized duplicate of every 
 * record (or just inter-chromosomal records). This is typically run after 
 * `BndCanonize.java` to make it adhere to VCF conventions.
 * 
 * Remarks: 
 * 1. the output VCF is not necessarily sorted;
 * 2. for simplicity, a symmetrized record uses N in REF and ALT;
 * 3. the ID of a symmetrized record is the original ID with "_sym" appended;
 * 4. for simplicity, BNDs are assumed to follow the simple form without 
 *    inserted sequence.
 */
public class BndSymmetrize {

    /**
     * @param args 1: symmetrize also intra-chromosomal BNDs (1=yes).
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final boolean SYMMETRIZE_INTRA_CHROMOSOMAL = Integer.parseInt(args[1])==1;
        
        int nRecordsInput, nRecordsSymmetrized, nRecordsOutput;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        nRecordsInput=0; nRecordsSymmetrized=0; nRecordsOutput=0;
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.startsWith("#CHROM")) System.out.println("##INFO=<ID=DUPLICATE,Number=0,Type=Flag,Description=\"The BND record is the symmetrized duplicate of an original record from a canonized VCF.\">");
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecordsInput++;
            System.out.println(str); nRecordsOutput++;
            tokens=str.split("\t");
            if (symmetrize(tokens,SYMMETRIZE_INTRA_CHROMOSOMAL)) {
                nRecordsSymmetrized++;
                System.out.println(String.join("\t",tokens));
                nRecordsOutput++;
            }
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecordsInput = "+nRecordsInput);
        System.err.println("nRecordsSymmetrized = "+nRecordsSymmetrized);
        System.err.println("nRecordsOutput = "+nRecordsOutput);
    }


    /**
     * Symmetrizes a BND record stored in `tokens` by changing only
     * CHROM,POS,ID,REF,ALT and adding a new flag to INFO.
     * 
     * @return true iff the record was symmetrized.
     */
    private static final boolean symmetrize(String[] tokens, boolean symmetrizeIntraChromosomal) {
        boolean refDirection, altDirection;
        char c, separator;
        int p, q;
        int first;
        String refChrom, refPos, alt, altChrom, altPos;

        refChrom=tokens[0]; refPos=tokens[1]; alt=tokens[4];

        // Extracting key quantities
        c=alt.charAt(0);
        refDirection=(c!='[')&&(c!=']');   // True = Left
        altDirection=alt.indexOf(']')>=0;  // True = Left
        p=alt.indexOf('['); q=alt.indexOf(']'); first=-1; separator='_';
        if (p>=0) { separator='['; first=p; }
        else if (q>=0) { separator=']'; first=q; }
        else {
            System.err.println("ERROR: unrecognized ALT = "+alt);
            System.exit(1);
        }
        p=alt.indexOf(':',first+1);
        altChrom=alt.substring(first+1,p);
        if (!symmetrizeIntraChromosomal && altChrom.equals(refChrom)) return false;
        q=alt.indexOf(separator,p+1);
        altPos=alt.substring(p+1,q);

        // Symmetrizing
        tokens[0]=altChrom; tokens[1]=altPos; tokens[2]+="_sym"; tokens[3]="N";
        separator=refDirection?']':'[';
        if (altDirection) tokens[4]="N"+separator+refChrom+":"+refPos+separator;
        else tokens[4]=separator+refChrom+":"+refPos+separator+"N";
        tokens[7]=tokens[7].equals(".")?"DUPLICATE":tokens[7]+";DUPLICATE";
        return true;
    }

}
