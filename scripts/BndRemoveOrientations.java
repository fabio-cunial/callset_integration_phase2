import java.util.*;
import java.io.*;
import java.util.zip.GZIPInputStream;


/**
 * Given a BND-only VCF, the program forces every ALT to have a fixed form 
 * `A[B:C[`, i.e. a fixed orientation on both sides. This is useful e.g. if one 
 * wants to filter BNDs without taking their orientations into account.
 */
public class BndRemoveOrientations {
    
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
            removeOrientations(tokens);
            System.out.println(String.join("\t",tokens));
            str=br.readLine();
        }
        br.close();
    }


    private static final void removeOrientations(String[] tokens) {
        int p, q;
        int first;
        String alt;
        
        alt=tokens[4];
        alt=alt.replace(']','[');
        p=alt.indexOf('[');
        if (p==0) {
            q=alt.indexOf('[',1);
            alt=alt.substring(q+1)+alt.substring(0,q+1);
        }
        tokens[4]=alt;
    }

}