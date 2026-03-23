import java.util.*;
import java.io.*;


/**
 * 
 */
public class UltralongIntervalIntersectClips {
    
    /**
     * @param args 
     * 1: 0-based, inclusive;
     * 2: 0-based, exclusive.
     */
    public static void main(String[] args) throws IOException {
        final String CLIPS1_FILE = args[0];
        final int N_CLIPS1 = Integer.parseInt(args[1]);
        final boolean CLIPS1_IS_LEFT_MAXIMAL = Integer.parseInt(args[2])==1;
        final String CLIPS2_FILE = args[3];
        final int N_CLIPS2 = Integer.parseInt(args[4]);
        final boolean CLIPS2_IS_LEFT_MAXIMAL = Integer.parseInt(args[5])==1;
        
        boolean isRc, matchFound;
        char c;
        int i, j;
        int refPos, readPos, cigarLength, readLength, length;
        String str, cigar, output;
        BufferedReader br;
        BufferedWriter bwLeft, bwRight;
        boolean[] marked1, marked2;
        String[] clips1, clips2;
        
        // Loading both files
        if (N_CLIPS1==0 || N_CLIPS2==0) {
            System.out.println("0");
            System.exit(0);
        }
        clips1 = new String[N_CLIPS1][4]; marked1 = new boolean[N_CLIPS1];
        br = new BufferedReader(new InputStreamReader(new FileInputStream(CLIPS1_FILE)));
        str=br.readLine(); i=-1;
        while (str!=null) {
            clips1[++i]=str.split("\t");
            str=br.readLine();
        }
        br.close();
        clips2 = new String[N_RIGHT_CLIPS][4]; marked2 = new boolean[N_CLIPS2];
        br = new BufferedReader(new InputStreamReader(new FileInputStream(CLIPS2_FILE)));
        str=br.readLine(); i=-1;
        while (str!=null) {
            clips2[++i]=str.split("\t");
            str=br.readLine();
        }
        br.close();
        
        // Naive implementation that performs all pairwise comparisons
        Arrays.fill(marked1,false); Arrays.fill(marked2,false);
        for (i=0; i<N_CLIPS1; i++) {
            readId1=clips1[i][0];
            isRc1=clips1[i][1].charAt(0)==1;
            readPos1=Integer.parseInt(clips1[i][2]);
            readLength1=Integer.parseInt(clips1[i][3]);
            for (j=0; j<N_CLIPS2; j++) {
                readId2=clips2[i][0];
                k=readId2.compareTo(readId1);
                if (k<0) continue;
                else if (k>0) break;
                isRc2=clips2[i][1].charAt(0)==1;
                readPos2=Integer.parseInt(clips2[i][2]);
                readLength2=Integer.parseInt(clips2[i][3]);
                if (canBeChained(isRc1,radPos1,readLength1,CLIPS1_IS_LEFT_MAXIMAL,isRc2,readPos2,readLength2,CLIPS2_IS_LEFT_MAXIMAL)) { marked1[i]=true; marked2[j]=true; }
            }
        }
        
        // Outputting
        
        
        
    }
    
    
    /**
     * @return TRUE iff two maximal alignments of the same read can be connected
     * with one another by traversing the read in some direction.
     */
    private static final boolean canBeChained(boolean isRc1, int pos1, int length1, boolean isLeftMaximal1, boolean isRc2, int pos2, int length2, boolean isLeftMaximal2) {
        if (!isRc1) {
            if (isLeftMaximal1) {
                if (isLeftMaximal2) return isRc2 && (length2-pos2<=pos1+slack && length2-pos2>=pos1-slack);
                else return !isRc2 && pos2<=pos1+slack && pos2>=pos1-slack;
            }
            else {
                if (isLeftMaximal2) !isRc2 && pos2<=pos1+slack && pos2>=pos1-slack;
                else return isRc2 && (length2-pos2<=pos1+slack && length2-pos2>=pos1-slack);
            }
        }
        else {
            if (isLeftMaximal1) {
                if (isLeftMaximal2) return !isRc2 && (pos2<=length1-pos1+slack && pos2>=length1-pos1-slack);
                else return isRc2 && pos2<=pos1+slack && pos2>=pos1-slack;
            }
            else {
                if (isLeftMaximal2) isRc2 && pos2<=pos1+slack && pos2>=pos1-slack;
                else return !isRc2 && (pos2<=length1-pos1+slack && pos2>=length1-pos1-slack);
            }
        }
    }

}