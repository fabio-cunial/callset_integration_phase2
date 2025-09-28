import java.io.*;
import java.text.NumberFormat;


/**
 * 
 */
public class SumCoverages {
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV = args[0];
        
        int i;
        double sum;
        String str;
        BufferedReader br;
        NumberFormat formatter;
        String[] tokens, tokensPrime;
        
        formatter=NumberFormat.getInstance();
        formatter.setMaximumFractionDigits(1);
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_TSV)));
        str=br.readLine(); System.out.println(str);
        str=br.readLine();
        while (str!=null) { 
            tokens=str.split("\t");
            tokensPrime=tokens[1].substring(1,tokens[1].length()-1).split(",");
            sum=0;
            for (i=0; i<tokensPrime.length; i++) sum+=Double.parseDouble(tokensPrime[i]);
            System.out.println(tokens[0]+"\t"+formatter.format(sum));
            str=br.readLine();
        }
        br.close();
    }
    
}
