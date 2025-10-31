import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * Makes sure that REF and ALT are uppercase and contain only characters in
 * {A,C,G,T,N}. Overwrites QUAL with a given constant and FILTER with PASS.
 */
public class CleanRefAltQual {
    
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
        StringBuilder buffer;
        BufferedReader br;
        String[] tokens;
        
        buffer = new StringBuilder();
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
            replaceNonstandardChars(tokens[3].toUpperCase(),buffer);
            tokens[3]=buffer.toString();
            replaceNonstandardChars(tokens[4].toUpperCase(),buffer);
            tokens[4]=buffer.toString();
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
    
    
	private static final void replaceNonstandardChars(String str, StringBuilder out) {
        final int LENGTH = str.length();
        
        char c;
        int i;
        
        out.delete(0,out.length());
        for (i=0; i<LENGTH; i++) {
            c=str.charAt(i);
            if (c=='A' || c=='C' || c=='G' || c=='T') out.append(c);
            else out.append('N');
        }
	}

}