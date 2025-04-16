import java.io.*;
import java.util.zip.*;


/**
 * 
 */
public class CopyFormat {
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String VCF_GZ_1 = args[0];
        final String VCF_GZ_2 = args[1];
        
        int i, j, p;
        int count;
        String str1, str2;
        BufferedReader br1, br2;
        String[] tokens1, tokens2;
        
        br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(VCF_GZ_1))));
        br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(VCF_GZ_2))));
        str1=br1.readLine();
        while (str1.charAt(0)=='#') str1=br1.readLine();
        str2=br2.readLine();
        while (str2.charAt(0)=='#') {
            System.out.println(str2);
            str2=br2.readLine();
        }
        while (str1!=null) {
            tokens1=str1.split("\t"); tokens2=str2.split("\t");
            if (tokens1.length!=tokens1.length) {
                System.err.println("ERROR: different number of columns in the two files:");
                System.err.println(str1);
                System.err.println(str2);
                System.exit(1);
            }
            for (i=9; i<tokens1.length; i++) {
                count=0; p=-1;
                for (j=0; j<tokens1[i].length(); j++) {
                    if (tokens1[i].charAt(j)!=':') continue;
                    count++;
                    if (count==9) { p=j; break; }
                }
                tokens2[i]+=tokens1[i].substring(p);
            }
            System.out.print(tokens2[0]);
            for (i=1; i<tokens2.length; i++) System.out.print("\t"+tokens2[i]);
            System.out.println();
            str1=br1.readLine(); str2=br2.readLine();
            if ((str1==null)!=(str2==null)) {
                System.err.println("ERROR: the two files do not contain the same number of lines.");
                System.exit(1);
            }
        }
        br1.close(); br2.close();
    }
    
}