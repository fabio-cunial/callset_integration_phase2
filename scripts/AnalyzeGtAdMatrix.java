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
        final double FRACTIONAL_QUANTUM = 0.01;
        
        boolean isMissing;
        int i, j, p, q;
        int nAlts, adRef, adAlt, nRows;
        double adRef_fractional, adAlt_fractional;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        long[][] histogram_ref, histogram_alt;
        long[][] histogram_ref_fractional, histogram_alt_fractional;
        String[] tokens;
        
        // Collecting stats
        histogram_ref = new long[4][AD_MAX+1];
        histogram_alt = new long[4][AD_MAX+1];
        histogram_ref_fractional = new long[4][(int)Math.ceil(1.0/FRACTIONAL_QUANTUM)+1];
        histogram_alt_fractional = new long[4][(int)Math.ceil(1.0/FRACTIONAL_QUANTUM)+1];
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
                adAlt=Integer.parseInt(tokens[i].substring(q+1));
                adRef_fractional=((double)adRef)/(adRef+adAlt);
                adAlt_fractional=((double)adAlt)/(adRef+adAlt);
                if (adRef>AD_MAX) adRef=AD_MAX;
                if (adAlt>AD_MAX) adAlt=AD_MAX;
                if (nAlts>0) {
                    histogram_ref[0][adRef]++; 
                    histogram_alt[0][adAlt]++; 
                    histogram_ref_fractional[0][(int)(adRef_fractional/FRACTIONAL_QUANTUM)]++;
                    histogram_alt_fractional[0][(int)(adAlt_fractional/FRACTIONAL_QUANTUM)]++;
                }
                if (nAlts==0) {
                    histogram_ref[1][adRef]++;
                    histogram_alt[1][adAlt]++;
                    histogram_ref_fractional[1][(int)(adRef_fractional/FRACTIONAL_QUANTUM)]++;
                    histogram_alt_fractional[1][(int)(adAlt_fractional/FRACTIONAL_QUANTUM)]++;
                }
                else if (nAlts==1) { 
                    histogram_ref[2][adRef]++; 
                    histogram_alt[2][adAlt]++; 
                    histogram_ref_fractional[2][(int)(adRef_fractional/FRACTIONAL_QUANTUM)]++;
                    histogram_alt_fractional[2][(int)(adAlt_fractional/FRACTIONAL_QUANTUM)]++;
                }
                else if (nAlts==2) { 
                    histogram_ref[3][adRef]++; 
                    histogram_alt[3][adAlt]++; 
                    histogram_ref_fractional[3][(int)(adRef_fractional/FRACTIONAL_QUANTUM)]++;
                    histogram_alt_fractional[3][(int)(adAlt_fractional/FRACTIONAL_QUANTUM)]++;
                }
            }
            nRows++;
            if (nRows%1000==0) System.err.println("Processed "+nRows+" rows...");
            if (nRows%100000==0) printHistograms(histogram_ref,histogram_alt,histogram_ref_fractional,histogram_alt_fractional);
            str=br.readLine();
        }
        br.close();
        printHistograms(histogram_ref,histogram_alt,histogram_ref_fractional,histogram_alt_fractional);
        
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
        bw = new BufferedWriter(new FileWriter("histogram_ref_fractional.txt"));
        for (i=0; i<histogram_ref_fractional.length; i++) {
            for (j=0; j<histogram_ref_fractional[i].length; j++) bw.write(histogram_ref_fractional[i][j]+",");
            bw.newLine();
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter("histogram_alt_fractional.txt"));
        for (i=0; i<histogram_alt_fractional.length; i++) {
            for (j=0; j<histogram_alt_fractional[i].length; j++) bw.write(histogram_alt_fractional[i][j]+",");
            bw.newLine();
        }
        bw.close();
    }
    
    
    /**
     * To STDERR
     */
    private static final void printHistograms(long[][] histogram_ref, long[][] histogram_alt, long[][] histogram_ref_fractional, long[][] histogram_alt_fractional) {
        int i, j;
        
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
        System.err.println("histogram_ref_fractional:");
        for (i=0; i<histogram_ref_fractional.length; i++) {
            for (j=0; j<histogram_ref_fractional[i].length; j++) System.err.print(histogram_ref_fractional[i][j]+",");
            System.err.println();
        }
        System.err.println("histogram_alt_fractional:");
        for (i=0; i<histogram_alt_fractional.length; i++) {
            for (j=0; j<histogram_alt_fractional[i].length; j++) System.err.print(histogram_alt_fractional[i][j]+",");
            System.err.println();
        }
    }
    
}