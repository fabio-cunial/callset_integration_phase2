import java.io.*;

/**
 * 
 */
public class KanpigScore2 {
    
    /**
     * @param 
     */
    public static void main(String[] args) throws IOException {
        final String KS_CSV = args[0];
        
        final int QUANTUM = 10000;  // Arbitrary
        
        int i;
        int ks1, ks2;
        String str;
        BufferedReader br;
        int[] histogram_true_ks, histogram_true_ks1, histogram_true_ks2, histogram_false_ks, histogram_false_ks1, histogram_false_ks2;
        String[] tokens;
        
        histogram_true_ks = new int[101];
        histogram_true_ks1 = new int[101];
        histogram_true_ks2 = new int[101];
        histogram_false_ks = new int[101];
        histogram_false_ks1 = new int[101];
        histogram_false_ks2 = new int[101];
        br = new BufferedReader(new FileReader(KS_CSV));
        str=br.readLine();
        while (str!=null) { 
            tokens=str.split(",");
            if (tokens.length==3 && tokens[2].equals(".")) {
                str=br.readLine();
                continue;
            }
            if (tokens.length==4) {
                ks1=Integer.parseInt(tokens[2]);
                ks2=Integer.parseInt(tokens[3]);
            }
            else {
                ks1=Integer.parseInt(tokens[2]);
                ks2=-1;
            }
            if (tokens[0].equals("1")) {
                if (ks2!=-1) { histogram_true_ks1[ks1]++; histogram_true_ks2[ks2]++; }
                else histogram_true_ks[ks1]++;
            }
            else {
                if (ks2!=-1) { histogram_false_ks1[ks1]++; histogram_false_ks2[ks2]++; }
                else histogram_false_ks[ks1]++;
            }
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        for (i=0; i<histogram_true_ks.length; i++) System.out.println(histogram_true_ks[i]+","+histogram_true_ks1[i]+","+histogram_true_ks2[i]+","+histogram_false_ks[i]+","+histogram_false_ks1[i]+","+histogram_false_ks2[i]);
    }
    
}