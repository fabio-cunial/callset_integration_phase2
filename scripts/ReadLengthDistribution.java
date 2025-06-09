import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class ReadLengthDistribution {
    
    /**
     * @param args
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TXT = args[0];
        final int MIN_LENGTH = Integer.parseInt(args[1]);
        final int MAX_LENGTH = Integer.parseInt(args[2]);
        final int QUANTUM = Integer.parseInt(args[3]);
        
        final int N_BINS = (MAX_LENGTH-MIN_LENGTH+QUANTUM-1)/QUANTUM;
        
        int i;
        int length;
        BufferedReader br;
        String str;
        long[] histogram;
        
        histogram = new long[N_BINS];
        br = new BufferedReader(new FileReader(INPUT_TXT));
        str=br.readLine();
        while (str!=null) {
            length=Integer.parseInt(str);
            i=(length-MIN_LENGTH)/QUANTUM;
            if (i<0) i=0;
            else if (i>=N_BINS) i=N_BINS-1;
            histogram[i]++;
            str=br.readLine();
        }
        br.close();
        for (i=0; i<N_BINS; i++) System.out.println((MIN_LENGTH+i*QUANTUM)+","+histogram[i]);
    }
    
}