import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class SumFileSizes {
    
    /**
     * @param args 0 the output of `gsutil ls -l`.
     */
    public static void main(String[] args) throws IOException {
        final String GSUTIL_LS_FILE = args[0];
        
        final long denominator = 1000000000;  // GB
        
        int p;
        long sum;
        String str;
        BufferedReader br;
        
        br = new BufferedReader(new FileReader(GSUTIL_LS_FILE));
        str=br.readLine(); sum=0L;
        while (str!=null) {
            p=str.indexOf(" ");
            sum+=Long.parseLong(str.substring(0,p));
            str=br.readLine();
        }
        br.close();
        System.out.println(""+(sum+denominator-1)/denominator);  // Ceil
    }
    
}