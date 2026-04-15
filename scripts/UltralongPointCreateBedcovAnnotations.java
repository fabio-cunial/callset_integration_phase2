import java.util.*;
import java.io.*;


/**
 * Given the BED file created by `samtools bedcov` (with format CHROM,START,END,
 * ID,BEDCOV), the program reformats it as VCF_ID,BIN. Output values are 
 * normalized by BIN_LENGTH.
 */
public class UltralongPointCreateBedcovAnnotations {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_BED = args[0];
        final int BIN_LENGTH = Integer.parseInt(args[1]);
        
        int i;
        String str, id;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_BED)));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            id=tokens[3].substring(0,tokens[3].lastIndexOf("_"));
            System.out.printf("%s\t%.3f\n",id,(Double.parseDouble(tokens[4])/BIN_LENGTH));
            str=br.readLine();
        }
        br.close();
    }

}