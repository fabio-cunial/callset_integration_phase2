import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * 
 */
public class GetGenotypingPriors {
    
    /**
     * See https://github.com/ACEnglish/truvari/wiki/MatchIds
     *
     * @param args 
     * 
     */
    public static void main(String[] args) throws IOException {
        final String TPCOMP_VCF_GZ = args[0];
        final String TPBASE_VCF_GZ = args[1];
        final String FP_VCF_GZ = args[2];
        final String OUTPUT_PREFIX = args[3];
        
        final int AD_MAX = 50;  // Arbitrary
        final int AD_INDEX = 7;
        final String MATCH_ID_STR = "MatchId=";
        final int MATCH_ID_STR_LENGTH = MATCH_ID_STR.length();
        
        int i, j, p, q;
        int nOnes, adRef, adAlt, record;
        String str, gt, matchIds;
        HashMap<String,String> baseGts;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens, tokensPrime;
        int[][] distribution_00, distribution_01, distribution_11;
        
        
        // Loading base GTs
        baseGts = new HashMap<String,String>();
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(TPBASE_VCF_GZ))));
        str=br.readLine(); record=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            p=tokens[9].indexOf(":");
            if (p!=3) {
                // Not an autosome
                record++;
                str=br.readLine();
                continue;
            }
            gt=tokens[9].substring(0,p);
            p=tokens[7].indexOf(MATCH_ID_STR);
            q=tokens[7].indexOf(";",p+MATCH_ID_STR_LENGTH);
            matchIds=tokens[7].substring(p+MATCH_ID_STR_LENGTH,q==-1?tokens[7].length():q);
            p=matchIds.indexOf(",");
            baseGts.put(matchIds.substring(0,p),gt);
            record++;
            if (record%1000==0) System.err.println("Processed "+record+" records...");
            str=br.readLine();
        }
        br.close();
        
        // Computing distributions for 0/1 and 1/1.
        distribution_01 = new int[AD_MAX+1][AD_MAX+1];
        distribution_11 = new int[AD_MAX+1][AD_MAX+1];
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(TPCOMP_VCF_GZ))));
        str=br.readLine(); record=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            p=tokens[9].indexOf(":");
            if (p!=3) {
                // Not an autosome
                record++;
                str=br.readLine();
                continue;
            }
            tokensPrime=tokens[9].split(":");
            if (tokensPrime[AD_INDEX].indexOf(".")>=0) {
                // No AD values
                record++;
                str=br.readLine();
                continue;
            }
            p=tokensPrime[AD_INDEX].indexOf(",");
            adRef=Integer.parseInt(tokensPrime[AD_INDEX].substring(0,p));
            adAlt=Integer.parseInt(tokensPrime[AD_INDEX].substring(p+1));
            p=tokens[7].indexOf(MATCH_ID_STR);
            q=tokens[7].indexOf(";",p+MATCH_ID_STR_LENGTH);
            matchIds=tokens[7].substring(p+MATCH_ID_STR_LENGTH,q==-1?tokens[7].length():q);
            p=matchIds.indexOf(",");
            gt=baseGts.get(matchIds.substring(0,p));
            if (gt==null) {
                System.err.println("ERROR: base MatchId not present in the base VCF: "+matchIds.substring(0,p));
                System.exit(1);
            }
            nOnes=(gt.charAt(0)=='1'?1:0)+(gt.charAt(2)=='1'?1:0);
            if (nOnes==1) distribution_01[adRef+adAlt>AD_MAX?AD_MAX:adRef+adAlt][adAlt>AD_MAX?AD_MAX:adAlt]++;
            else if (nOnes==2) distribution_11[adRef+adAlt>AD_MAX?AD_MAX:adRef+adAlt][adAlt>AD_MAX?AD_MAX:adAlt]++;
            record++;
            if (record%1000==0) System.err.println("Processed "+record+" records...");
            str=br.readLine();
        }
        br.close();
        
        // Computing distributions for 0/0.
        distribution_00 = new int[AD_MAX+1][AD_MAX+1];
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(FP_VCF_GZ))));
        str=br.readLine(); record=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            p=tokens[9].indexOf(":");
            if (p!=3) {
                // Not an autosome
                record++;
                str=br.readLine();
                continue;
            }
            tokensPrime=tokens[9].split(":");
            if (tokensPrime[AD_INDEX].indexOf(".")>=0) {
                // No AD values
                record++;
                str=br.readLine();
                continue;
            }
            p=tokensPrime[AD_INDEX].indexOf(",");
            adRef=Integer.parseInt(tokensPrime[AD_INDEX].substring(0,p));
            adAlt=Integer.parseInt(tokensPrime[AD_INDEX].substring(p+1));            
            distribution_00[adRef+adAlt>AD_MAX?AD_MAX:adRef+adAlt][adAlt>AD_MAX?AD_MAX:adAlt]++;
            record++;
            if (record%1000==0) System.err.println("Processed "+record+" records...");
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_01.csv"));
        for (i=0; i<distribution_01.length; i++) {
            for (j=0; j<distribution_01[i].length; j++) bw.write(distribution_01[i][j]+",");
            bw.newLine();
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_11.csv"));
        for (i=0; i<distribution_11.length; i++) {
            for (j=0; j<distribution_11[i].length; j++) bw.write(distribution_11[i][j]+",");
            bw.newLine();
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_00.csv"));
        for (i=0; i<distribution_00.length; i++) {
            for (j=0; j<distribution_00[i].length; j++) bw.write(distribution_00[i][j]+",");
            bw.newLine();
        }
        bw.close();
    }
    
}