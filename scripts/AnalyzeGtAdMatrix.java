import java.io.*;


/** 
 * Remark: the procedure considers only diploid GTs without any `.`.
 */
public class AnalyzeGtAdMatrix {
    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV = args[0];
        
        final int AD_MAX = 100;
        
        boolean isMissing;
        int i, j, p, q;
        int nAlts, adRef, adAlt, nRows;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        long[][] histogram_ref, histogram_alt;
        String[] tokens;
        
        // Collecting stats
        histogram_ref = new long[4][AD_MAX+1];
        histogram_alt = new long[4][AD_MAX+1];
        br = new BufferedReader(new FileReader(INPUT_TSV));
        str=br.readLine(); nRows=0;
        while (str!=null) {
            tokens=str.split("\t");
            for (i=0; i<tokens.length; i++) {
                p=tokens[i].indexOf(",");
                if (p!=3 || tokens[i].charAt(0)=='.' || tokens[i].charAt(2)=='.') continue;
                nAlts=(tokens[i].charAt(0)=='1'?1:0)+(tokens[i].charAt(2)=='1'?1:0);
                q=tokens[i].indexOf(",",p+1);
                adRef=Integer.parseInt(tokens[i].substring(p+1,q));
                if (adRef>AD_MAX) adRef=AD_MAX;
                adAlt=Integer.parseInt(tokens[i].substring(q+1));
                if (adAlt>AD_MAX) adAlt=AD_MAX;
                if (nAlts>0) { histogram_ref[0][adRef]++; histogram_alt[0][adAlt]++; }
                if (nAlts==0) { histogram_ref[1][adRef]++; histogram_alt[1][adAlt]++; }
                else if (nAlts==1) { histogram_ref[2][adRef]++; histogram_alt[2][adAlt]++; }
                else if (nAlts==2) { histogram_ref[3][adRef]++; histogram_alt[3][adAlt]++; }
            }
            nRows++;
            if (nRows%1000==0) System.err.println("Processed "+nRows+" rows...");
            if (nRows%100000==0) {
                System.err.println("histogram_ref:");
                for (i=0; i<histogram_ref.length; i++) {
                    for (j=0; j<histogram_ref[i].length; j++) System.err.print(histogram_ref[i][j]+",");
                    System.err.println();
                }
                System.err.println("histogram_alt:");
                for (i=0; i<histogram_alt.length; i++) {
                    for (j=0; j<histogram_alt[i].length; j++) System.err.print(histogram_alt[i][j]+",");
                    System.err.println();
                }
            }
            str=br.readLine();
        }
        br.close();
        
        // Writing the full histograms to files
        bw = new BufferedWriter(new FileWriter("histogram_ref.txt"));
        for (i=0; i<histogram_ref.length; i++) {
            for (j=0; j<histogram_ref[i].length; j++) bw.write(histogram_ref[i][j]+",");
            bw.newLine();
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter("histogram_alt.txt"));
        for (i=0; i<histogram_alt.length; i++) {
            for (j=0; j<histogram_alt[i].length; j++) bw.write(histogram_alt[i][j]+",");
            bw.newLine();
        }
        bw.close();
    }
    
}