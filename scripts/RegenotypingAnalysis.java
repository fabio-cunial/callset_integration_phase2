import java.util.zip.GZIPInputStream;
import java.io.*;


/** 
 * Measures the effects of re-genotyping an inter-sample VCF. Downstream
 * analyses can extract just the fields they need by doing e.g.: 
 * 
 * bcftools query -f '%COUNT1,%COUNT2' annotated.vcf
 */
public class RegenotypingAnalysis {
    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String BEFORE_REGENOTYPING_VCF_GZ = args[0];
        final String AFTER_REGENOTYPING_VCF_GZ = args[1];
        
        int i, j, k;
        int record;
        String str1, str2;
        BufferedReader br1, br2;
        int[] counts1 = new int[3];
        int[] counts2 = new int[3];
        String[] tokens1, tokens2;
        int[][] deltas = new int[7][7];
        
        br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(BEFORE_REGENOTYPING_VCF_GZ))));
        str1=br1.readLine();
        while (str1!=null) {
            if (str1.charAt(0)!='#') break;
            if (str1.substring(0,6).equalsIgnoreCase("#CHROM")) {
                System.out.println("##INFO=<ID=COUNT0,Number=1,Type=Integer,Description=\"after-before for GT=1\">");
                System.out.println("##INFO=<ID=COUNT1,Number=1,Type=Integer,Description=\"after-before for GT=0/1\">");
                System.out.println("##INFO=<ID=COUNT2,Number=1,Type=Integer,Description=\"after-before for GT=1/1\">");
                for (j=0; j<deltas.length; j++) {
                    for (k=0; k<deltas[j].length; k++) System.out.println("##INFO=<ID=DELTA_"+j+"_"+k+",Number=1,Type=Integer,Description=\"after-before for every possible pair of GTs\">");
                }
            }
            System.out.println(str1);
            str1=br1.readLine();
        }
        br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(AFTER_REGENOTYPING_VCF_GZ))));
        str2=br2.readLine();
        while (str2!=null) {
            if (str2.charAt(0)!='#') break;
            str2=br2.readLine();
        }
        record=0;
        while (str1!=null && str2!=null) {
            tokens1=str1.split("\t"); tokens2=str2.split("\t");
            if (tokens1.length!=tokens2.length) {
                System.err.println("ERROR: different number of samples at record "+record);
                System.exit(1);
            }
            getCounts(tokens1,tokens2,counts1,counts2,deltas);
            for (j=0; j<counts1.length; j++) tokens1[7]+=";COUNT"+j+"="+(counts2[j]-counts1[j]);
            for (j=0; j<deltas.length; j++) {
                for (k=0; k<deltas[j].length; k++) tokens1[7]+=";DELTA_"+j+"_"+k+"="+deltas[j][k];
            }
            // Outputting a VCF with only the new annotations and no samples,
            // to speed up downstream analyses.
            System.out.print(tokens1[0]);
            for (i=1; i<=7; i++) System.out.print("\t"+tokens1[i]);
            System.out.println("\tGT\t0/1");
            record++;
            if (record%1000==0) System.err.println("Processed "+record+" records...");
            str1=br1.readLine(); str2=br2.readLine();
        }
        br1.close(); br2.close();
    }
    
    
    /**
     * @param counts* (output array) number of GTs equal to: 
     * 0    1
     * 1    0/1
     * 2    1/1
     * in `tokens*`;
     * @param deltas (output array) every row and column corresponds to a GT, as
     * follows:
     * 0    .
     * 1    0
     * 2    1
     * 3    ./.
     * 4    0/0
     * 5    0/1
     * 6    1/1
     */
    private static final void getCounts(String[] tokens1, String[] tokens2, int[] counts1, int[] counts2, int[][] deltas) {
        final int nColumns = tokens1.length;
        
        int i, j, p;
        int type1, type2;
        
        for (i=0; i<counts1.length; i++) counts1[i]=0;
        for (i=0; i<counts2.length; i++) counts2[i]=0;
        for (i=0; i<deltas.length; i++) {
            for (j=0; j<deltas[i].length; j++) deltas[i][j]=0;
        }
        for (i=9; i<nColumns; i++) {
            // First array
            type1=-1;
            p=-1;
            for (j=0; j<tokens1[i].length(); j++) {
                if (tokens1[i].charAt(j)==':') { p=j; break; }
            }
            if (p==1 || (p==-1 && tokens1[i].length()==1)) {
                if (tokens1[i].charAt(0)=='.') type1=0;
                else if (tokens1[i].charAt(0)=='0') type1=1;
                else if (tokens1[i].charAt(0)=='1') type1=2;
            }
            else if (p==3 || (p==-1 && tokens1[i].length()>=3)) {
                if (tokens1[i].charAt(0)=='.' && tokens1[i].charAt(2)=='.') type1=3;
                else if (tokens1[i].charAt(0)=='0' && tokens1[i].charAt(2)=='0') type1=4;
                else if ((tokens1[i].charAt(0)=='0' && tokens1[i].charAt(2)=='1') || (tokens1[i].charAt(0)=='1' && tokens1[i].charAt(2)=='0')) type1=5;
                else if (tokens1[i].charAt(0)=='1' && tokens1[i].charAt(2)=='1') type1=6;
            }
            else {
                System.err.println("ERROR: wrong ':' position at record "+i);
                System.exit(1);
            }
            // Second array
            type2=-1;
            p=-1;
            for (j=0; j<tokens2[i].length(); j++) {
                if (tokens2[i].charAt(j)==':') { p=j; break; }
            }
            if (p==1 || (p==-1 && tokens2[i].length()==1)) {
                if (tokens2[i].charAt(0)=='.') type2=0;
                else if (tokens2[i].charAt(0)=='0') type2=1;
                else if (tokens2[i].charAt(0)=='1') type2=2;
            }
            else if (p==3 || (p==-1 && tokens2[i].length()>=3)) {
                if (tokens2[i].charAt(0)=='.' && tokens2[i].charAt(2)=='.') type2=3;
                else if (tokens2[i].charAt(0)=='0' && tokens2[i].charAt(2)=='0') type2=4;
                else if ((tokens2[i].charAt(0)=='0' && tokens2[i].charAt(2)=='1') || (tokens2[i].charAt(0)=='1' && tokens2[i].charAt(2)=='0')) type2=5;
                else if (tokens2[i].charAt(0)=='1' && tokens2[i].charAt(2)=='1') type2=6;
            }
            else {
                System.err.println("ERROR: wrong ':' position at record "+i);
                System.exit(1);
            }
            // Stats
            if (type1==2) counts1[0]++;
            else if (type1==5) counts1[1]++;
            else if (type1==6) counts1[2]++;
            if (type2==2) counts2[0]++;
            else if (type2==5) counts2[1]++;
            else if (type2==6) counts2[2]++;
            if (type1!=-1 && type2!=-1) deltas[type1][type2]++;
        }
    }
    
}