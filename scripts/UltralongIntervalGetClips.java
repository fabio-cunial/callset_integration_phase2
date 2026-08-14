import java.util.*;
import java.io.*;


/**
 * Given a SAM file and a window, the program outputs two files containing a
 * record for every alignment that has a left (resp. right) soft- or hard-clip
 * at a position inside the window. We call such alignments left- (resp. right-)
 * maximal.
 *
 * Output format: READ_ID,IS_RC,READ_POS,READ_LENGTH
 * 
 * The program also outputs a file with the number of CIGAR INS inside the 
 * window, and with the number of CIGAR DEL that start or end inside the window:
 * 
 * Output format: N_INS_i,N_DEL_START_i,N_DEL_END_i for every i
 */
public class UltralongIntervalGetClips {
    
    /**
     * @param args 
     * 1: 0-based, inclusive;
     * 2: 0-based, exclusive;
     * 4: comma-separated, not necessarily sorted sequence; only CIGAR INDELs
     *    >=this are counted.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_SAM = args[0];
        final int START = Integer.parseInt(args[1]);
        final int END = Integer.parseInt(args[2]);
        final int MIN_CLIP_LENGTH = Integer.parseInt(args[3]);
        final String[] MIN_INDEL_LENGTH = args[4].split(",");
        final String OUTPUT_PREFIX = args[5];
        
        final int RC_MASK = 16;
        
        boolean isRc, matchFound;
        char c;
        int i, j, k;
        int refPos, readPos, newReadPos, cigarLength, readLength, length;
        String str, cigar, output;
        BufferedReader br;
        BufferedWriter bwLeft, bwRight, bwIndel;
        int[] minIndelLengths, nIns, nDelStart, nDelEnd;
        String[] tokens;

        minIndelLengths = new int[MIN_INDEL_LENGTH.length];
        for (i=0; i<MIN_INDEL_LENGTH.length; i++) minIndelLengths[i]=Integer.parseInt(MIN_INDEL_LENGTH[i]);
        nIns = new int[minIndelLengths.length]; 
        nDelStart = new int[minIndelLengths.length]; 
        nDelEnd = new int[minIndelLengths.length];
        bwLeft = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_leftmaximal.txt"));
        bwRight = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_rightmaximal.txt"));
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_SAM)));
        str=br.readLine(); 
        while (str!=null) {
            tokens=str.split("\t");
            refPos=Integer.parseInt(tokens[3]);  // 1-based, first mapping pos.
            cigar=tokens[5];
            cigarLength=cigar.length();
            isRc=(Integer.parseInt(tokens[1])&RC_MASK)!=0;
            readLength=tokens[9].length();

            // Updating readLength using hardclips
            i=0; c=cigar.charAt(i);
            while (i<cigarLength && c>='0' && c<='9') { i++; c=cigar.charAt(i); }
            if (c=='H') readLength+=Integer.parseInt(cigar.substring(0,i));
            if (cigar.charAt(cigarLength-1)=='H') {
                i=cigar.length()-2; c=cigar.charAt(i);
                while (i>=0 && c>='0' && c<='9') { i--; c=i>=0?cigar.charAt(i):'0'; }
                readLength+=Integer.parseInt(cigar.substring(i+1,cigarLength-1));
            }

            // Parsing the whole CIGAR
            readPos=0;  // 0-based, first mapping pos.
            i=0; matchFound=false;
            for (j=1; j<cigarLength; j++) {
                c=cigar.charAt(j);
                if (c=='M' || c=='=' || c=='X') {
                    matchFound=true;
                    length=Integer.parseInt(cigar.substring(i,j));
                    refPos+=length; readPos+=length;
                    i=j+1;
                }
                else if (c=='I') {
                    length=Integer.parseInt(cigar.substring(i,j));
                    if (refPos-1>=START && refPos-1<END) {
                        for (k=0; k<minIndelLengths.length; k++) {
                            if (length>=minIndelLengths[k]) nIns[k]++;
                        }
                    }
                    readPos+=length;
                    i=j+1;
                }
                else if (c=='D') {
                    matchFound=true;
                    length=Integer.parseInt(cigar.substring(i,j));
                    if (refPos-1>=START && refPos-1<END) {
                        for (k=0; k<minIndelLengths.length; k++) {
                            if (length>=minIndelLengths[k]) nDelStart[k]++;
                        }
                    }
                    refPos+=length;
                    if (refPos-1>=START && refPos-1<END) {
                        for (k=0; k<minIndelLengths.length; k++) {
                            if (length>=minIndelLengths[k]) nDelEnd[k]++;
                        }
                    }
                    i=j+1;
                }
                else if (c=='N') {
                    refPos+=Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                }
                else if (c=='S' || c=='H') {
                    newReadPos=readPos+Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                    if ( refPos-1>=START && refPos-1<END &&
                         ( (matchFound && readLength-readPos>=MIN_CLIP_LENGTH) || (!matchFound && newReadPos>=MIN_CLIP_LENGTH) )
                       ) {
                        output=tokens[0]+"\t"+(isRc?"1\t":"0\t")+(matchFound?readPos:newReadPos)+"\t"+(readLength+"\n");
                        if (!matchFound) bwLeft.write(output);
                        else bwRight.write(output);
                    }
                    readPos=newReadPos;
                }
                else if (c=='P') i=j+1;
            }

            // Next iteration
            str=br.readLine();
        }
        br.close(); bwLeft.close(); bwRight.close();

        // Outputting
        bwIndel = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_indel.txt"));
        bwIndel.write(nIns[0]+"\t"+nDelStart[0]+"\t"+nDelEnd[0]);
        for (k=1; k<minIndelLengths.length; k++) bwIndel.write("\t"+nIns[k]+"\t"+nDelStart[k]+"\t"+nDelEnd[k]);
        bwIndel.newLine();
        bwIndel.close();
    }

}