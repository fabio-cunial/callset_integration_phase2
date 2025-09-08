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
        
        boolean isMissing;
        int i, j, x, y;
        int nTrios, nRecords, nAltChild, adAltChild, adAltFather, adAltMother, adRefChild, adRefFather, adRefMother;
        int minDepth, maxDepth, avgDepth;
        String str, gtChild, gtFather, gtMother;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens, tokensPrime1, tokensPrime2, tokensPrime3;
        int[] numerator, denominator;  // De novo rates for each trio
        long[][] depth_to_numerator, depth_to_denominator;
        // Number of de novos with a given GT, AD_ALT_CHILD=X, and
        // AD_ALT_FATHER+AD_ALT_MOTHER=Y.
        long[][] denovo_01_adAlt, denovo_11_adAlt;
        long[][] nondenovo_01_adAlt, nondenovo_11_adAlt;
        // Number of de novos with a given GT, AD_REF_CHILD=X, and
        // AD_REF_FATHER+AD_REF_MOTHER=Y.
        long[][] denovo_01_adRef, denovo_11_adRef;
        long[][] nondenovo_01_adRef, nondenovo_11_adRef;
        
        // Computing the number of trios
        br = new BufferedReader(new FileReader(TRIO_MATRIX_TSV));
        str=br.readLine();
        tokens=str.split("\t");
        br.close();
        nTrios=tokens.length/3;
        System.err.println(nTrios+" trios");
        
        // Allocating space
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
        depth_to_numerator = new long[MAX_AD+1][nTrios];
        depth_to_denominator = new long[MAX_AD+1][nTrios];
        
        // Computing counts
        br = new BufferedReader(new FileReader(TRIO_MATRIX_TSV));
        nRecords=0;
        str=br.readLine();
        while (str!=null) { 
            tokens=str.split("\t");
            for (i=0; i<nTrios; i++) {
                tokensPrime1=tokens[3*i+0].split(",");
                gtChild=tokensPrime1[0];
                tokensPrime2=tokens[3*i+1].split(",");
                gtFather=tokensPrime2[0];
                tokensPrime3=tokens[3*i+2].split(",");
                gtMother=tokensPrime3[0];
                isMissing=gtFather.indexOf(".")>=0||gtMother.indexOf(".")>=0;
                if (isMissing) {
                    // At least one GT in the triplet is missing
                    continue;
                }
                if (gtChild.indexOf("1")<0 || gtChild.length()!=3 || gtFather.length()!=3 || gtMother.length()!=3) {
                    // Must occur in the child and be on an autosome.
                    continue;
                }
                adRefChild=Integer.parseInt(tokensPrime1[1]); adAltChild=Integer.parseInt(tokensPrime1[2]);
                adRefFather=Integer.parseInt(tokensPrime2[1]); adAltFather=Integer.parseInt(tokensPrime2[2]);
                adRefMother=Integer.parseInt(tokensPrime3[1]); adAltMother=Integer.parseInt(tokensPrime3[2]);
                denominator[i]++;
                minDepth=adRefChild+adAltChild; maxDepth=minDepth; avgDepth=minDepth;
                if (adRefFather+adAltFather<minDepth) minDepth=adRefFather+adAltFather;
                if (adRefMother+adAltMother<minDepth) minDepth=adRefMother+adAltMother;
                if (adRefFather+adAltFather>maxDepth) maxDepth=adRefFather+adAltFather;
                if (adRefMother+adAltMother>maxDepth) maxDepth=adRefMother+adAltMother;
                avgDepth+=adRefFather+adAltFather+adRefMother+adAltMother;
                avgDepth/=3;
                depth_to_denominator[(adRefChild+adAltChild)>MAX_AD?MAX_AD:(adRefChild+adAltChild)][i]++;
                if (gtFather.indexOf("1")<0 && gtMother.indexOf("1")<0) {
                    numerator[i]++;
                    depth_to_numerator[(adRefChild+adAltChild)>MAX_AD?MAX_AD:(adRefChild+adAltChild)][i]++;
                    // De novo
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
        for (i=0; i<nTrios; i++) System.out.println(((double)numerator[i])/denominator[i]);
        printMatrix(denovo_01_adAlt,TRIO_MATRIX_TSV+"_denovo_01_adAlt.txt");
        printMatrix(denovo_11_adAlt,TRIO_MATRIX_TSV+"_denovo_11_adAlt.txt");
        printMatrix(denovo_01_adRef,TRIO_MATRIX_TSV+"_denovo_01_adRef.txt");
        printMatrix(denovo_11_adRef,TRIO_MATRIX_TSV+"_denovo_11_adRef.txt");
        printMatrix(nondenovo_01_adAlt,TRIO_MATRIX_TSV+"_nondenovo_01_adAlt.txt");
        printMatrix(nondenovo_11_adAlt,TRIO_MATRIX_TSV+"_nondenovo_11_adAlt.txt");
        printMatrix(nondenovo_01_adRef,TRIO_MATRIX_TSV+"_nondenovo_01_adRef.txt");
        printMatrix(nondenovo_11_adRef,TRIO_MATRIX_TSV+"_nondenovo_11_adRef.txt");
        bw = new BufferedWriter(new FileWriter(TRIO_MATRIX_TSV+"_depth_to_denovo_rate.txt"));
        for (i=0; i<depth_to_numerator.length; i++) {
            for (j=0; j<depth_to_numerator[i].length; j++) bw.write((depth_to_denominator[i][j]==0?"0":(((double)depth_to_numerator[i][j])/depth_to_denominator[i][j]))+",");
            bw.newLine();
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter(TRIO_MATRIX_TSV+"_depth_to_ncalls.txt"));
        for (i=0; i<depth_to_numerator.length; i++) {
            for (j=0; j<depth_to_numerator[i].length; j++) bw.write((depth_to_numerator[i][j]+depth_to_denominator[i][j])+",");
            bw.newLine();
        }
        bw.close();
    }
    
    
    private static final void printMatrix(long[][] matrix, String fileName) throws IOException {
        int i, j;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(fileName));
        for (i=0; i<matrix.length; i++) {
            for (j=0; j<matrix[i].length; j++) bw.write(matrix[i][j]+",");
            bw.newLine();
        }
        bw.close();
    }
    
}