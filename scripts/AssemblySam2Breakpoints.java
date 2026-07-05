import java.util.Vector;
import java.io.*;


/**
 * Given an RNAME-sorted assembly-to-ref SAM, the program prints a CHR,POS
 * breakpoint for every clipped first/last position of an alignment that is far
 * enough from the previous/next alignment.
 * 
 * Remark: only standard chromosomes are printed in output.
 */
public class AssemblySam2Breakpoints {
    
    /**
     * @param args
     * 0: all alignments with the same RNAME are assumed to form a contiguous
     *    block.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_SAM = args[0];
        final int MIN_ALIGNMENT_DISTANCE = Integer.parseInt(args[1]);
        final int MIN_CLIP_LENGTH = Integer.parseInt(args[2]);
        final String OUTPUT_FILE = args[3];
        
        char c;
        int i, j, p, q;
        int refPos, refPosOriginal, readPos, cigarLength, leftClipLength, rightClipLength, readLength, length, nRecords;
        String str, readId, chrId, lastChrId, cigar, mapq;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        Vector<Alignment> alignments;
        
        alignments = new Vector<Alignment>();
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_SAM)));
        bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(OUTPUT_FILE)));
        str=br.readLine(); nRecords=0; lastChrId="";
        while (str!=null) {
            if (str.charAt(0)=='@') { str=br.readLine(); continue; }
            p=str.indexOf('\t');
            readId=str.substring(0,p);
            q=str.indexOf('\t',p+1);
            p=q; q=str.indexOf('\t',p+1);
            chrId=str.substring(p+1,q);
            if (lastChrId.length()>0 && !chrId.equals(lastChrId)) {
                if (isStandardChromosome(lastChrId)) printBreakpoints(lastChrId,alignments,MIN_ALIGNMENT_DISTANCE,MIN_CLIP_LENGTH,bw);
                alignments.clear();
            }
            lastChrId=chrId;
            p=q; q=str.indexOf('\t',p+1);
            refPos=Integer.parseInt(str.substring(p+1,q))-1;  // 0-based, first matching pos.
            refPosOriginal=refPos;
            refPos--;  // Moving the pointer right before the first matching pos.
            p=q; q=str.indexOf('\t',p+1);
            mapq=str.substring(p+1,q);
            p=q; q=str.indexOf('\t',p+1);
            cigar=str.substring(p+1,q); cigarLength=cigar.length();
            readPos=-1;  // 0-based
            i=0; leftClipLength=0; rightClipLength=0;
            for (j=1; j<cigarLength; j++) {
                c=cigar.charAt(j);
                if (c=='M' || c=='=' || c=='X') {
                    length=Integer.parseInt(cigar.substring(i,j));
                    refPos+=length; readPos+=length;
                    i=j+1;
                }
                else if (c=='I') {
                    readPos+=Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                }
                else if (c=='D') {
                    refPos+=Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                }
                else if (c=='N') {
                    refPos+=Integer.parseInt(cigar.substring(i,j));
                    i=j+1;
                }
                else if (c=='S' || c=='H') {
                    length=Integer.parseInt(cigar.substring(i,j));
                    if (i==0) leftClipLength=length;
                    else if (j==cigarLength-1) rightClipLength=length;
                    readPos+=length;
                    i=j+1;
                }
                else if (c=='P') i=j+1;
            }
            readLength=readPos+1;
            alignments.add(new Alignment(refPosOriginal,refPos,leftClipLength,rightClipLength));
            nRecords++;
            if (nRecords%1000==0) System.err.println("Processed "+nRecords+" records...");
            str=br.readLine();
        }
        br.close();
        if (lastChrId.length()>0) {
            if (isStandardChromosome(lastChrId)) printBreakpoints(lastChrId,alignments,MIN_ALIGNMENT_DISTANCE,MIN_CLIP_LENGTH,bw);
            alignments.clear();
        }
        bw.close();
    }


    /**
     * For every alignment, prints its first (resp. last) position iff it is 
     * clipped and it is far enough from the previous (resp. next) alignment.
     * 
     * Remark: this is a simple quadratic scan just for simplicity. Could be 
     * made much faster.
     * 
     * @param alignments assumed to be all and only the alignments in a chr;
     * @param out format: CHR,POS.
     */
    private static final void printBreakpoints(String chr, Vector<Alignment> alignments, int minDistance, int minClipLength, BufferedWriter out) throws IOException {
        int i, j;
        int min, max;
        Alignment a, b;

        alignments.sort(null);
        for (i=0; i<alignments.size(); i++) {
            a=alignments.get(i);
            if (a.leftClipLength>=minClipLength) {
                max=-1;
                for (j=i-1; j>=0; j--) {
                    b=alignments.get(j);
                    if (b.chrLast>max) max=b.chrLast;
                    if (max>=a.chrFirst) break;
                }
                if (max==-1 || a.chrFirst-max>=minDistance) { out.write(chr+","+a.chrFirst); out.newLine(); }
            }
            if (a.rightClipLength>=minClipLength) {
                min=Integer.MAX_VALUE;
                for (j=i+1; j<alignments.size(); j++) {
                    b=alignments.get(j);
                    if (b.chrFirst<min) min=b.chrFirst;
                    if (min<=a.chrLast || b.chrFirst>a.chrLast) break;
                }
                if (min==Integer.MAX_VALUE || min-a.chrLast>=minDistance) { out.write(chr+","+a.chrLast); out.newLine(); }
            }
        }
    }


    /**
     * Designed to be fast albeit not completely correct.
     */
    private static final boolean isStandardChromosome(String str) {
        if (!str.startsWith("chr")) return false;
        final int length = str.length();
        char c;

        if (length==4) {
            c=str.charAt(3);
            return c=='X' || c=='Y' || c=='M' || c=='x' || c=='y' || c=='m' || (c>='1' && c<='9');
        }
        else if (length==5) {
            c=str.charAt(3);
            if (c!='1' && c!='2') return false;
            c=str.charAt(4);
            return c>='0' && c<='9';
        }
        else return false;
    }


    private static class Alignment implements Comparable {
        int chrFirst, chrLast;
        int leftClipLength, rightClipLength;

        public Alignment(int chrFirst, int chrLast, int leftClipLength, int rightClipLength) {
            this.chrFirst=chrFirst; this.chrLast=chrLast;
            this.leftClipLength=leftClipLength; this.rightClipLength=rightClipLength;
        }

        public int compareTo(Object o) {
            Alignment a = (Alignment)o;
            if (this.chrFirst<a.chrFirst) return -1;
            else if (this.chrFirst>a.chrFirst) return 1;
            if (this.chrLast<a.chrLast) return -1;
            else if (this.chrLast>a.chrLast) return 1;
            return 0;
        }
    }

}