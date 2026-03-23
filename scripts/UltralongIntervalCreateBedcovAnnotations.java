import java.util.*;
import java.io.*;


/**
 * Given the BED file created by `samtools bedcov` (with format CHROM,START,END,
 * ID,BEDCOV), the program writes all the bins of the same VCF record to a
 * single row with format VCF_ID,BIN_1,BIN_2,...,BIN_N. Output values are 
 * normalized by SVLEN.
 */
public class UltralongIntervalCreateBedcovAnnotations {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_BED = args[0];
        final int N_BINS = Integer.parseInt(args[1]);
        
        int i;
        int svlen, currentBin, pos;
        String str, id, currentID, chrom;
        BufferedReader br;
        double[] bins;
        String[] tokens;
        
        bins = new double[N_BINS];
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_BED)));
        str=br.readLine(); currentID=""; currentBin=-1;
        while (str!=null) {
            tokens=str.split("\t");
            chrom=tokens[0];
            pos=Integer.parseInt(tokens[1]);
            id=tokens[3].substring(0,tokens[3].lastIndexOf("_"));
            if (!id.equals(currentID)) {
                if (currentID.length()!=0) {
                    System.out.print(currentID);
                    for (i=0; i<N_BINS; i++) System.out.printf("\t%.3f",bins[i]);
                    System.out.println();
                }
                currentID=id; currentBin=-1;
                Arrays.fill(bins,0.0);
            }
            svlen=Integer.parseInt(tokens[2])-Integer.parseInt(tokens[1]);
            bins[++currentBin]=Double.parseDouble(tokens[4])/svlen;

            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.out.print(currentID);
        for (i=0; i<N_BINS; i++) System.out.printf("\t%.3f",bins[i]);
        System.out.println();
    }

}