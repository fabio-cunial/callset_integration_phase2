import java.io.*;
import java.util.zip.GZIPInputStream;


/**
 * Given a BND-only VCF, the program symmetrizes every record.
 * This is typically run after `BndCanonize.java`.
 * 
 * Remarks: 
 * 1. the output VCF is not necessarily sorted;
 * 2. for simplicity, a symmetrized record uses N in REF and ALT;
 * 3. for simplicity, BNDs are assumed to follow the simple form without 
 *    inserted sequence).
 */
public class BndSymmetrize {

    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            symmetrize(tokens);
            System.out.println(String.join("\t", tokens));
            str=br.readLine();
        }
        br.close();
    }


    /**
     * Symmetrizes a BND record stored in `tokens` by changing only
     * CHROM,POS,REF,ALT and adding the SYMMETRIZED flag to INFO.
     */
    private static final void symmetrize(String[] tokens) {
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
        q=alt.indexOf(separator,p+1);
        altPos=alt.substring(p+1,q);

        // Symmetrizing
        tokens[0]=altChrom; tokens[1]=altPos; tokens[3]="N";
        separator=refDirection?']':'[';
        if (altDirection) tokens[4]="N"+separator+refChrom+":"+refPos+separator;
        else tokens[4]=separator+refChrom+":"+refPos+separator+"N";
        tokens[7]=tokens[7].equals(".")?"SYMMETRIZED":tokens[7]+";SYMMETRIZED";
    }

}