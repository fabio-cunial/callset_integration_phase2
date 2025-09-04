import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class CheckDeNovo {
    
    /**
     * @param args
     * 0: a tab-separated matrix where each call has format `GT,AD_REF,AD_ALT`.
     */
    public static void main(String[] args) throws IOException {
        final String TRIO_MATRIX_TSV = args[0];
        final int MAX_AD = Integer.parseInt(args[1]);
        
        
        boolean isMissing, isAltChild, isDeNovo;
        int i, j;
        int nTrios, nRecords;
        String str;
        BufferedReader br;
        String[] tokens, tokensPrime;
        int[] numerator, denominator;  // De novo rates for each trio
        // Number of de novos with a given GT, AD_ALT_CHILD=X, and
        // AD_ALT_FATHER+AD_ALT_MOTHER=Y.
        long[][] denovo_01_adAlt, denovo_11_adAlt;
        long[][] nondenovo_01_adAlt, denovo_11_adAlt;
        // Number of de novos with a given GT, AD_REF_CHILD=X, and
        // AD_REF_FATHER+AD_REF_MOTHER=Y.
        long[][] denovo_01_adRef, denovo_11_adRef;
        long[][] denovo_01_adRef, denovo_11_adRef;
        
        
        // Computing the number of trios
        br = new BufferedReader(new FileReader(TRIO_MATRIX_TSV));
        str=br.readLine();
        tokens=str.split("\t");
        br.close();
        nTrios=tokens.length/3;
        System.err.println(nTrios+" trios");
        
        // Computing rates
        numerator = new int[nTrios];
        denominator = new int[nTrios];
        denovo_01_adAlt = new long[MAX_AD+1][MAX_AD+1];
        denovo_11_adAlt = new long[MAX_AD+1][MAX_AD+1];
        denovo_01_adRef = new long[MAX_AD+1][MAX_AD+1];
        denovo_11_adRef = new long[MAX_AD+1][MAX_AD+1];
        
        nondenovo_01_adAlt = new long[MAX_AD+1][MAX_AD+1];
        nondenovo_11_adAlt = new long[MAX_AD+1][MAX_AD+1];
        nondenovo_01_adRef = new long[MAX_AD+1][MAX_AD+1];
        nondenovo_11_adRef = new long[MAX_AD+1][MAX_AD+1];
        
        br = new BufferedReader(new FileReader(TRIO_MATRIX_TSV));
        nRecords=0;
        str=br.readLine();
        while (str!=null) { 
            tokens=str.split("\t");
            for (i=0; i<nTrios; i++) {
                tokensPrime=tokens[3*i+0].split(",");
                gtChild=tokensPrime[0]; adRefChild=Integer.parseInt(tokensPrime[1]); adAltChild=Integer.parseInt(tokensPrime[2]);
                tokensPrime=tokens[3*i+1].split(",");
                gtFather=tokensPrime[0]; adRefFather=Integer.parseInt(tokensPrime[1]); adAltFather=Integer.parseInt(tokensPrime[2]);
                tokensPrime=tokens[3*i+2].split(",");
                gtMother=tokensPrime[0]; adRefMother=Integer.parseInt(tokensPrime[1]); adAltMother=Integer.parseInt(tokensPrime[2]);
                if (gtChild.indexOf("1")<0 || gtChild.length()!=3 || gtFather.length()!=3 || gtMother.length()!=3) {
                    // Must occur in the child and be on an autosome.
                    continue;
                }
                isMissing=gtFather.indexOf(".")>=0||gtMother.indexOf(".")>=0;
                if (isMissing) {
                    // At least one GT in the triplet is missing
                    continue;
                }
                denominator[i]++;
                if (gtFather.indexOf("1")<0 && gtMother.indexOf("1")<0) {
                    // De novo
                    numerator[i]++;
                    nAltChild=(gtChild.charAt(0)=='1'?1:0)+(gtChild.charAt(2)=='1'?1:0);
                    x=adAltChild>MAX_AD?MAX_AD:adAltChild;
                    y=adAltFather+adAltMother;
                    if (y>MAX_AD) y=MAX_AD;
                    if (nAltChild==1) denovo_01_adAlt[x][y]++;
                    else if (nAltChild==2) denovo_11_adAlt[x][y]++;
                    x=adRefChild>MAX_AD?MAX_AD:adRefChild;
                    y=adRefFather+adRefMother;
                    if (y>MAX_AD) y=MAX_AD;
                    if (nAltChild==1) denovo_01_adRef[x][y]++;
                    else if (nAltChild==2) denovo_11_adRef[x][y]++;
                }
                else {
                    // Not a de novo
                    nAltChild=(gtChild.charAt(0)=='1'?1:0)+(gtChild.charAt(2)=='1'?1:0);
                    x=adAltChild>MAX_AD?MAX_AD:adAltChild;
                    y=adAltFather+adAltMother;
                    if (y>MAX_AD) y=MAX_AD;
                    if (nAltChild==1) nondenovo_01_adAlt[x][y]++;
                    else if (nAltChild==2) nondenovo_11_adAlt[x][y]++;
                    x=adRefChild>MAX_AD?MAX_AD:adRefChild;
                    y=adRefFather+adRefMother;
                    if (y>MAX_AD) y=MAX_AD;
                    if (nAltChild==1) nondenovo_01_adRef[x][y]++;
                    else if (nAltChild==2) nondenovo_11_adRef[x][y]++;
                }
            }
            nRecords++;
            if (nRecords%10000==0) System.err.println("Processed "+nRecords+" records...");
            str=br.readLine(); 
        }
        br.close();
        
        // Outputting
        --------------->
        for (i=0; i<nTrios; i++) System.out.println(((double)numerator[i])/denominator[i]);
    }
    
}