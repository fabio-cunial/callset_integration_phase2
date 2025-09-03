import java.io.*;


/** 
 * 
 */
public class AnalyzeGtAdMatrix {
    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV = args[0];
        
        final int AD_MAX = 50;
        
        boolean isMissing;
        int i, j, p, q;
        int nAlts, adRef, adAlt;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        long[][] histogram_ref, histogram_alt;
        String[] tokens;
        
        // Collecting stats
        histogram_ref = new int[4][AD_MAX+1];
        histogram_alt = new int[4][AD_MAX+1];
        br = new BufferedReader(new FileReader(INPUT_TSV));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            for (i=0; i<tokens.length; i++) {
                p=tokens[i].indexOf(",");
                isMissing=false; nAlts=0;
                if (p==3) {
                    isMissing=tokens[i].charAt(0)=='.' && tokens[i].charAt(2)=='.';
                    nAlts=(tokens[i].charAt(0)=='1'?1:0)+(tokens[i].charAt(2)=='1'?1:0);
                }
                else if (p==1) {
                    isMissing=tokens[i].charAt(0)=='.';
                    nAlts=tokens[i].charAt(0)=='1'?1:0;
                }
                if (isMissing) continue;
                q=tokens[i].indexOf(",",p+1);
                adRef=Integer.parseInt(tokens[i].substring(p+1,q));
                if (adRef>AD_MAX) adRef=AD_MAX;
                adAlt=Integer.parseInt(tokens[i].substring(q+1));
                if (adAlt>AD_MAX) adAlt=AD_MAX;
                histogram_ref[0][adRef]++; histogram_alt[0][adAlt]++;
                if (nAlts==0) { histogram_ref[1][adRef]++; histogram_alt[1][adAlt]++; }
                else if (nAlts==1) { histogram_ref[2][adRef]++; histogram_alt[2][adAlt]++; }
                else if (nAlts==2) { histogram_ref[3][adRef]++; histogram_alt[3][adAlt]++; }
            }
            str=br.readLine();
        }
        br.close();
        
        // Outputting
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