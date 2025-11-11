import java.io.*;

/**
 * 
 */
public class KanpigScore {
    
    /**
     * @param 
     */
    public static void main(String[] args) throws IOException {
        final String KS_CSV = args[0];
        
        final int QUANTUM = 10000;  // Arbitrary
        
        int i, p, q;
        int ks1, ks2;
        String str, genotype;
        BufferedReader br;
        int[] histogram_present_ks, histogram_present_ks1, histogram_present_ks2, histogram_absent_ks, histogram_absent_ks1, histogram_absent_ks2;
        String[] tokens;
        
        histogram_present_ks = new int[101];
        histogram_present_ks1 = new int[101];
        histogram_present_ks2 = new int[101];
        histogram_absent_ks = new int[101];
        histogram_absent_ks1 = new int[101];
        histogram_absent_ks2 = new int[101];
        br = new BufferedReader(new FileReader(KS_CSV));
        str=br.readLine();
        while (str!=null) { 
            p=str.indexOf(",");
            if (str.substring(p+1).equals(".")) {
                str=br.readLine();
                continue;
            }
            genotype=str.substring(0,p);
            q=str.indexOf(",",p+1);
            if (q>=0) {
                ks1=Integer.parseInt(str.substring(p+1,q));
                ks2=Integer.parseInt(str.substring(q+1));
            }
            else {
                ks1=Integer.parseInt(str.substring(p+1));
                ks2=-1;
            }
            if (genotype.indexOf("1")>=0) {
                if (ks2!=-1) { histogram_present_ks1[ks1]++; histogram_present_ks2[ks2]++; }
                else histogram_present_ks[ks1]++;
            }
            else {
                if (ks2!=-1) { histogram_absent_ks1[ks1]++; histogram_absent_ks2[ks2]++; }
                else histogram_absent_ks[ks1]++;
            }
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        for (i=0; i<histogram_present_ks.length; i++) System.out.println(histogram_present_ks[i]+","+histogram_present_ks1[i]+","+histogram_present_ks2[i]+","+histogram_absent_ks[i]+","+histogram_absent_ks1[i]+","+histogram_absent_ks2[i]);
    }
    
}