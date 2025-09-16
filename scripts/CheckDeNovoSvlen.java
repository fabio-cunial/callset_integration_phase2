import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class CheckDeNovoSvlen {
    
    /**
     * @param args
     * 0: a tab-separated matrix where each row has format:
     * `SVTYPE,SVLEN,GT,AD_REF,AD_ALT \t GT,AD_REF,AD_ALT \t GT,AD_REF,AD_ALT`.
     */
    public static void main(String[] args) throws IOException {
        final String TRIO_MATRIX_TSV = args[0];
        final int MAX_AD = Integer.parseInt(args[1]);
        
        final int MAX_LENGTH = 10000;
        final int[] LENGTHS = new int[] {20,50,100,200,300,400,500,600,700,800,1000,2000,3000,4000,5000,6000,7000,8000,9000,MAX_LENGTH};
        final int N_LENGTHS = LENGTHS.length;
        
        
        boolean isMissing;
        int i, j, x, y, p, q;
        int nTrios, nRecords, nAltChild, adAltChild, adAltFather, adAltMother, adRefChild, adRefFather, adRefMother;
        int minDepth, maxDepth, avgDepth, svlen, depthIndex, lengthIndex;
        String str, gtChild, gtFather, gtMother, svtype;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens, tokensPrime1, tokensPrime2, tokensPrime3;
        long[][] ins_depth_to_numerator, ins_depth_to_denominator;
        long[][] del_depth_to_numerator, del_depth_to_denominator;
        long[][] ins_len_to_numerator, ins_len_to_denominator;
        long[][] del_len_to_numerator, del_len_to_denominator;
        
        // Computing the number of trios
        br = new BufferedReader(new FileReader(TRIO_MATRIX_TSV));
        str=br.readLine();
        tokens=str.split("\t");
        br.close();
        nTrios=tokens.length/3;
        System.err.println(nTrios+" trios");
        
        // Allocating space
        ins_depth_to_numerator = new long[MAX_AD+1][nTrios];
        ins_depth_to_denominator = new long[MAX_AD+1][nTrios];
        del_depth_to_numerator = new long[MAX_AD+1][nTrios];
        del_depth_to_denominator = new long[MAX_AD+1][nTrios];
        ins_len_to_numerator = new long[N_LENGTHS][nTrios];
        ins_len_to_denominator = new long[N_LENGTHS][nTrios];
        del_len_to_numerator = new long[N_LENGTHS][nTrios];
        del_len_to_denominator = new long[N_LENGTHS][nTrios];
        
        // Computing counts
        br = new BufferedReader(new FileReader(TRIO_MATRIX_TSV));
        nRecords=0;
        str=br.readLine();
        while (str!=null) { 
            tokens=str.split("\t");
            p=tokens[0].indexOf(",");
            svtype=tokens[0].substring(0,p);
            q=tokens[0].indexOf(",",p+1);
            svlen=Integer.parseInt(tokens[0].substring(p+1,q));
            tokens[0]=tokens[0].substring(q+1);
            for (i=0; i<nTrios; i++) {
                tokensPrime1=tokens[3*i+0].split(",");
                gtChild=tokensPrime1[0];
                tokensPrime2=tokens[3*i+1].split(",");
                gtFather=tokensPrime2[0];
                tokensPrime3=tokens[3*i+2].split(",");
                gtMother=tokensPrime3[0];
                isMissing=gtChild.indexOf(".")>=0||gtFather.indexOf(".")>=0||gtMother.indexOf(".")>=0;
                if (isMissing) {
                    // At least one GT in the triplet is missing
                    continue;
                }
                if (gtChild.indexOf("1")<0 || gtChild.length()!=3 || gtFather.length()!=3 || gtMother.length()!=3) {
                    // Must occur in the child and be on an autosome.
                    continue;
                }
                adRefChild=tokensPrime1[1].indexOf(".")>=0?0:Integer.parseInt(tokensPrime1[1]); 
                adAltChild=tokensPrime1[2].indexOf(".")>=0?0:Integer.parseInt(tokensPrime1[2]);
                adRefFather=tokensPrime2[1].indexOf(".")>=0?0:Integer.parseInt(tokensPrime2[1]); 
                adAltFather=tokensPrime2[2].indexOf(".")>=0?0:Integer.parseInt(tokensPrime2[2]);
                adRefMother=tokensPrime3[1].indexOf(".")>=0?0:Integer.parseInt(tokensPrime3[1]); 
                adAltMother=tokensPrime3[2].indexOf(".")>=0?0:Integer.parseInt(tokensPrime3[2]);
                minDepth=adRefChild+adAltChild; maxDepth=minDepth; avgDepth=minDepth;
                if (adRefFather+adAltFather<minDepth) minDepth=adRefFather+adAltFather;
                if (adRefMother+adAltMother<minDepth) minDepth=adRefMother+adAltMother;
                if (adRefFather+adAltFather>maxDepth) maxDepth=adRefFather+adAltFather;
                if (adRefMother+adAltMother>maxDepth) maxDepth=adRefMother+adAltMother;
                avgDepth+=adRefFather+adAltFather+adRefMother+adAltMother;
                avgDepth/=3;
                depthIndex=avgDepth>MAX_AD?MAX_AD:avgDepth;
                lengthIndex=Arrays.binarySearch(LENGTHS,0,N_LENGTHS,svlen);
                if (lengthIndex<0) lengthIndex=0;
                if (lengthIndex>N_LENGTHS-1) lengthIndex=N_LENGTHS-1;
                if (svtype.equalsIgnoreCase("INS")) {
                    ins_depth_to_denominator[depthIndex][i]++;
                    ins_len_to_denominator[lengthIndex][i]++;
                }
                else if (svtype.equalsIgnoreCase("DEL")) {
                    del_depth_to_denominator[depthIndex][i]++;
                    del_len_to_denominator[lengthIndex][i]++;
                }
                else {
                    // Not an INS or DEL
                    continue;
                }
                if (gtFather.indexOf("1")<0 && gtMother.indexOf("1")<0) {
                    // De novo
                    if (svtype.equalsIgnoreCase("INS")) {
                        ins_depth_to_numerator[depthIndex][i]++;
                        ins_len_to_numerator[lengthIndex][i]++;
                    }
                    else {
                        del_depth_to_numerator[depthIndex][i]++;
                        del_len_to_numerator[lengthIndex][i]++;
                    }
                }
                else {
                    // Not a de novo
                }
            }
            nRecords++;
            if (nRecords%10000==0) System.err.println("Processed "+nRecords+" records...");
            str=br.readLine(); 
        }
        br.close();
        
        // Outputting
        printMatrix(ins_depth_to_numerator,ins_depth_to_denominator,TRIO_MATRIX_TSV+"_ins_depth_to_denovo_rate.txt");
        printMatrix(del_depth_to_numerator,del_depth_to_denominator,TRIO_MATRIX_TSV+"_del_depth_to_denovo_rate.txt");
        printMatrix(ins_len_to_numerator,ins_len_to_denominator,TRIO_MATRIX_TSV+"_ins_len_to_denovo_rate.txt");
        printMatrix(del_len_to_numerator,del_len_to_denominator,TRIO_MATRIX_TSV+"_del_len_to_denovo_rate.txt");
    }
    
    
    private static final void printMatrix(long[][] numerator_matrix, long[][] denominator_matrix, String fileName) throws IOException {
        int i, j;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(fileName));
        for (i=0; i<numerator_matrix.length; i++) {
            for (j=0; j<numerator_matrix[i].length; j++) bw.write((denominator_matrix[i][j]==0?"0":(((double)numerator_matrix[i][j])/denominator_matrix[i][j]))+",");
            bw.newLine();
        }
        bw.close();
    }
    
}