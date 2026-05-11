import java.util.Arrays;
import java.util.Vector;
import java.io.*;


/**
 * Given the output of `samtools depth` over a region, the program computes a
 * simple estimate of the leftmost and rightmost positions that delimit an 
 * increase in coverage.
 */
public class UltralongDepthGetBreakpoints {
    
    /**
     * @args
     */
    public static void main(String[] args) throws IOException {
        final String SAMTOOLS_DEPTH_FILE = args[0];
        final int SAMTOOLS_DEPTH_N_POSITIONS = Integer.parseInt(args[1]);
        final int BIN_LENGTH = Integer.parseInt(args[2]);
        final double BIN_COVERAGE_RATIO = Double.parseDouble(args[3]);
        
        int i, j, p;
        int length, sum1, sum2, leftBreakpoint, rightBreakpoint, posStart, posEnd, validPairs, maxLength, bestLeft, bestRight, lowDepth1, lowDepth2;
        long firstPos;
        String str;
        BufferedReader br;
        short[] depth;
        Vector<Integer> leftBreakpoints, rightBreakpoints, lowDepths1, lowDepths2;

        // Loading all depths in memory
        depth = new short[SAMTOOLS_DEPTH_N_POSITIONS]; firstPos=-1;
        br = new BufferedReader(new FileReader(SAMTOOLS_DEPTH_FILE));
        str=br.readLine(); i=-1;
        while (str!=null) {
            length=str.length();
            p=0;
            while (p<length) {
                if (str.charAt(p)=='\t') break;
                else p++;
            }
            p++; posStart=p;
            while (p<length) {
                if (str.charAt(p)=='\t') break;
                else p++;
            }
            posEnd=p;
            if (firstPos==-1) firstPos=Integer.parseInt(str.substring(posStart,posEnd));
            p++;
            depth[++i]=Short.parseShort(str.substring(p));
            str=br.readLine();
        }
        br.close();

        // Collecting all left breakpoints, if any.
        leftBreakpoints = new Vector<Integer>();
        lowDepths1 = new Vector<Integer>();
        sum1=0;
        for (i=0; i<BIN_LENGTH; i++) sum1+=depth[i];
        sum2=0;
        for (i=BIN_LENGTH; i<2*BIN_LENGTH; i++) sum2+=depth[i];
        if (sum2>=sum1*BIN_COVERAGE_RATIO) {
            leftBreakpoints.add(i-BIN_LENGTH);
            lowDepths1.add(sum1);
        }
        for (i=2*BIN_LENGTH; i<SAMTOOLS_DEPTH_N_POSITIONS; i++) {
            sum1+=depth[i-BIN_LENGTH]-depth[i-2*BIN_LENGTH];
            sum2+=depth[i]-depth[i-BIN_LENGTH];
            if (sum2>=sum1*BIN_COVERAGE_RATIO) {
                leftBreakpoints.add(i-BIN_LENGTH);
                lowDepths1.add(sum1);
            }
        }
        if (leftBreakpoints.isEmpty()) return;
        System.err.println(leftBreakpoints.size()+" left breakpoints");

        // Collecting all right breakpoints, if any.
        rightBreakpoints = new Vector<Integer>();
        lowDepths2 = new Vector<Integer>();
        sum2=0;
        for (i=SAMTOOLS_DEPTH_N_POSITIONS-1; i>=SAMTOOLS_DEPTH_N_POSITIONS-BIN_LENGTH; i--) sum2+=depth[i];
        sum1=0;
        for (i=SAMTOOLS_DEPTH_N_POSITIONS-BIN_LENGTH-1; i>=SAMTOOLS_DEPTH_N_POSITIONS-2*BIN_LENGTH; i--) sum1+=depth[i];
        if (sum1>=sum2*BIN_COVERAGE_RATIO) {
            rightBreakpoints.add(i+BIN_LENGTH);
            lowDepths2.add(sum2);
        }
        for (i=SAMTOOLS_DEPTH_N_POSITIONS-2*BIN_LENGTH-1; i>=0; i--) {
            sum1+=depth[i]-depth[i+BIN_LENGTH];
            sum2+=depth[i+BIN_LENGTH]-depth[i+2*BIN_LENGTH];
            if (sum1>=sum2*BIN_COVERAGE_RATIO) {
                rightBreakpoints.add(i+BIN_LENGTH);
                lowDepths2.add(sum2);
            }
        }
        if (rightBreakpoints.isEmpty()) return;
        System.err.println(rightBreakpoints.size()+" right breakpoints");
        
        // Finding a longest valid pair, if any.
        // Remark: medians are computed in a very naive way and should be made 
        // faster.
        validPairs=0; maxLength=0; bestLeft=-1; bestRight=-1;
        for (i=0; i<leftBreakpoints.size(); i++) {
            leftBreakpoint=leftBreakpoints.get(i);
            lowDepth1=lowDepths1.get(i);
            for (j=0; j<rightBreakpoints.size(); j++) {
                rightBreakpoint=rightBreakpoints.get(j);
                if (rightBreakpoint<=leftBreakpoint) continue;
                lowDepth2=lowDepths2.get(j);
                short[] newArray = Arrays.copyOfRange(depth,leftBreakpoint,rightBreakpoint+1);
                Arrays.sort(newArray);
                if (newArray[newArray.length/2]<(lowDepth1<lowDepth2?lowDepth1:lowDepth2)*BIN_COVERAGE_RATIO) return;
                validPairs++;
                if (rightBreakpoint-leftBreakpoint>maxLength) {
                    maxLength=rightBreakpoint-leftBreakpoint;
                    bestLeft=leftBreakpoint; bestRight=rightBreakpoint;
                }
            }
        }
        System.err.println("validPairs="+validPairs+" maxLength="+maxLength);
        if (validPairs==0) return;

        // Outputting
        System.out.println((firstPos+bestLeft)+"\t"+(firstPos+bestRight));
    }

}