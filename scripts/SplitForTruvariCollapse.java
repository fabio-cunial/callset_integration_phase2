import java.io.*;

/**
 * Partitions a chromosome into consecutive chunks of a given minimum length.
 *
 * Remark: the output is a 0-based, half-open CSV that covers the entire chr.  
 * To split the bcftools merge VCF in preparation of a parallel truvari collapse
 * one should use the BED as follows:
 *
 * bcftools view --regions-file one_line_of_bed.bed --regions-overlap pos
 * 
 * Remark: a better approach would partition the chromosome by keeping the
 * number of calls in each element of the partition balanced.
 *
 * @param args
 * 0: chromosome ID;
 * 1: BED file of truvari chunks created by `density_counter.py`; all chunks are
 *    used, including those with too many calls;
 * 3: min length of a truvari collapse chunk (chunks in output might be longer);
 * 4: the program extends an input chunk [start,end) to [start,end+X), just for 
 *    boundary safety.
 */
public class SplitForTruvariCollapse {
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String CHROMOSOME = args[0];
        final String INPUT_BED = args[1];
        final String INPUT_FAI = args[2];
        final int CHUNK_LENGTH_BP = Integer.parseInt(args[3]);
        final int SLACK = Integer.parseInt(args[4]);
        
        int i;
        int from, to, currentFrom, currentTo, nInputChunks, chrLength;
        String str;
        BufferedReader br;
        String[] tokens;
        
        // Loading chromosome length
        br = new BufferedReader(new FileReader(INPUT_FAI));
        str=br.readLine();
        chrLength=0;
        while (str!=null) {
            tokens=str.split("\t");
            if (tokens[0].equals(CHROMOSOME)) {
                chrLength=Integer.parseInt(tokens[1]);
                break;
            }
            str=br.readLine();
        }
        br.close();
        if (chrLength==0) {
            System.err.println("ERROR: "+CHROMOSOME+" nt found in "+INPUT_FAI);
            System.exit(1);
        }
        
        // Loading the number of input chunks
        br = new BufferedReader(new FileReader(INPUT_BED));
        str=br.readLine();
        nInputChunks=0;
        while (str!=null) {
            nInputChunks++;
            str=br.readLine();
        }
        br.close();
        
        // Computing output chunks
        br = new BufferedReader(new FileReader(INPUT_BED));
        currentFrom=0; currentTo=0; tokens=null;
        for (i=0; i<nInputChunks; i++) {
            str=br.readLine();
            tokens=str.split("\t");
            from=Integer.parseInt(tokens[1]); to=Integer.parseInt(tokens[2]);
            if (to>currentTo) currentTo=to;
            if (currentTo-currentFrom>=CHUNK_LENGTH_BP) {
                if (i==nInputChunks-1) currentTo=chrLength;
                System.out.println(tokens[0]+","+currentFrom+","+(currentTo+SLACK));
                currentFrom=currentTo+SLACK; currentTo=currentFrom;
            }
        }
        br.close();
        if (currentTo<chrLength) System.out.println(tokens[0]+","+currentFrom+","+chrLength);
    }
    
}