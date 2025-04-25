import java.io.*;

/**
 * Splits every chromosome into consecutive intervals of fixed size.
 *
 * Remark: the output is a 0-based, half-open BED. To split a single-sample VCF
 * in preparation of a parallel bcftools merge, one should use the BED as 
 * follows:
 *
 * bcftools view --regions-file one_line_of_bed.bed --regions-overlap pos
 */
public class SplitForBcftoolsMerge {
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FAI = args[0];
        final int CHUNK_LENGH_BP = Integer.parseInt(args[1]);
        
        int i;
        int chrLength;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(INPUT_FAI));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            chrLength=Integer.parseInt(tokens[1]);
            i=0;
            while (i<chrLength) {
                System.out.println(tokens[0]+"\t"+i+"\t"+Math.min(i+CHUNK_LENGH_BP,chrLength));
                i+=CHUNK_LENGH_BP;
            }
            str=br.readLine();
        }
        br.close();
    }
    
}