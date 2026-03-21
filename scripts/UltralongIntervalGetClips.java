import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Creates two output files with the read IDs of all reads with a left (resp.
 * right) soft- or hard-clip.
 */
public class UltralongIntervalGetClips {
    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_SAM = args[0];
        final int START = Integer.parseInt(args[1]);  // 0-based, inclusive.
        final int END = Integer.parseInt(args[2]);  // 0-based, exclusive.
        final String OUTPUT_PREFIX = args[3];
        
        boolean matchFound;
        char c;
        int i, j;
        int pos, cigarLength;
        String str, cigar;
        BufferedReader br;
        BufferedWriter bwLeft, bwRight;
        String[] tokens;
        
        bwLeft = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_left.clips"));
        bwRight = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_right.clips"));
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_SAM)));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            pos=Integer.parseInt(tokens[3]);  // 1-based, first match.
            cigar=tokens[5];
            cigarLength=cigar.length();
            i=0; matchFound=false;
            for (j=1; j<cigarLength; j++) {
                c=cigar.charAt(j);
                if (c=='S' || c=='H') {
                    if (pos-1>=START && pos-1<END) {
                        if (!matchFound) { bwLeft.write(tokens[0]); bwLeft.newLine(); }
                        else  { bwRight.write(tokens[0]); bwRight.newLine(); }
                    }
                    i=j+1;
                }
                else if (c=='M' || c=='D' || c=='N' || c=='=' || c=='X') {
                    matchFound=true;
                    pos+=Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                }
                else if (c=='I' || c=='P') i=j+1;
            }

            // Next iteration
            str=br.readLine();
        }
        br.close(); bwLeft.close(); bwRight.close();
    }

}