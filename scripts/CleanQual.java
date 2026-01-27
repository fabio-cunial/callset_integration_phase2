import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * Overwrites QUAL with a given constant, and overwrites FILTER with PASS.
 */
public class CleanQual {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String FORCE_QUAL = args[1];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        int i;
        int nRecords;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine(); nRecords=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");
            tokens[5]=FORCE_QUAL;
            tokens[6]="PASS";
            
            // Outputting
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) { System.out.print('\t'); System.out.print(tokens[i]); }
            System.out.println();
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
    }

}