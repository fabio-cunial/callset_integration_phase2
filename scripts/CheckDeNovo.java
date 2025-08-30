import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class CheckDeNovo {
    
    /**
     * @param args
     * 
     */
    public static void main(String[] args) throws IOException {
        final String TRIO_MATRIX_CSV = args[0];
        
        final int QUANTUM = 10000;  // Arbitrary
        
        boolean isMissing, isAltChild, isDeNovo;
        int i, j;
        int nTrios, nRecords;
        String str;
        BufferedReader br;
        String[] tokens;
        int[] numerator, denominator;
        
        
        // Computing the number of trios
        br = new BufferedReader(new FileReader(TRIO_MATRIX_CSV));
        str=br.readLine();
        tokens=str.split(",");
        br.close();
        nTrios=tokens.length/3;
        System.err.println(nTrios+" trios");
        
        // Computing rates
        numerator = new int[nTrios];
        Arrays.fill(numerator,0);
        denominator = new int[nTrios];
        Arrays.fill(denominator,0);
        br = new BufferedReader(new FileReader(TRIO_MATRIX_CSV));
        nRecords=0;
        str=br.readLine();
        while (str!=null) { 
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records...");
            tokens=str.split(",");
            for (i=0; i<nTrios; i++) {
                if (tokens[3*i].indexOf("1")<0) {
                    // Not ALT in the child
                    continue;
                }
                /*isMissing=false;
                for (j=3*i; j<3*i+3; j++) {
                    try { 
                        if (tokens[j].indexOf(".")>=0) { isMissing=true; break; }
                    } catch(Exception e) { System.err.println("ERROR in trio "+i+"-th: "+str); }
                }
                if (isMissing) {
                    // At least one GT in the triplet is missing
                    continue;
                }*/
                denominator[i]++;
                if (tokens[3*i+1].indexOf("1")<0 && tokens[3*i+2].indexOf("1")<0) numerator[i]++;
            }
            str=br.readLine(); 
        }
        br.close();
        
        // Outputting
        for (i=0; i<nTrios; i++) System.out.println(((double)numerator[i])/denominator[i]);
    }
    
}