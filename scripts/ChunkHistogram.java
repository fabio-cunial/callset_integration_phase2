import java.util.Arrays;
import java.io.*;

/**
 * 
 */
public class ChunkHistogram {
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String CHUNKS_BED = args[0];
        final int[] BINS = new int[] {1,2,3,4,5,6,7,8,9, 10,20,30,40,50,60,70,80,90, 100,200,300,400,500,600,700,800,900, 1000,2000,3000,4000,5000,6000,7000,8000,9000, 10000};
        
        int i;
        String str;
        BufferedReader br;
        long[] histogram;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(CHUNKS_BED));
        histogram = new long[BINS.length];
        Arrays.fill(histogram,0);
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            i=Arrays.binarySearch(BINS,Integer.parseInt(tokens[3]));
            if (i<0) i=-i-1;
            if (i==BINS.length) i=BINS.length-1;
            histogram[i]++;
            str=br.readLine();
        }
        br.close();
        for (i=0; i<BINS.length; i++) System.out.println(BINS[i]+","+histogram[i]);
    }
    
}