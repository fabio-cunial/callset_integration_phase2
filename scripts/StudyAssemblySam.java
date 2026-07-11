import java.util.Vector;
import java.util.HashMap;
import java.io.*;


/**
 * Given a QNAME-sorted assembly-to-ref SAM, the program prints a record for
 * every alignment with just enough information to build a permutation/
 * rearrangement diagram.
 */
public class StudyAssemblySam {
    
    /**
     * @param args
     * 0: all alignments with the same QNAME are assumed to form a contiguous
     *    block (e.g. `samtools sort -@ 8 -m 3G -n -O BAM hap1.bam`);
     * 1: opened for append;
     * 2: opened for append.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_SAM = args[0];
        final String OUTPUT_LENGTHS_FILE = args[1];
        final String OUTPUT_MATRIX_FILE = args[2];
        final int DISTANCE_THRESHOLD = Integer.parseInt(args[3]);

        final int RC_MASK = 16;
        
        boolean isRc;
        char c;
        int i, j, p, q;
        int refPos, refPosOriginal, readPos, cigarLength, leftClipLength, rightClipLength, readLength, length, nRecords;
        String str, readId, chrId, lastReadId, cigar, mapq;
        BufferedReader br;
        BufferedWriter bwMatrix, bwLengths;
        String[] tokens;
        Vector<Alignment> alignments;
        HashMap<String,Integer> tmpMap;
        
        alignments = new Vector<Alignment>();
        tmpMap = new HashMap<String,Integer>();
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_SAM)));
        bwMatrix = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(OUTPUT_MATRIX_FILE,true)));
        bwLengths = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(OUTPUT_LENGTHS_FILE,true)));
        str=br.readLine(); nRecords=0; lastReadId="";
        while (str!=null) {
            if (str.charAt(0)=='@') { str=br.readLine(); continue; }
            p=str.indexOf('\t');
            readId=str.substring(0,p);
            if (lastReadId.length()>0 && !readId.equals(lastReadId)) {
                printLengths(alignments,DISTANCE_THRESHOLD,tmpMap,bwLengths);
                printAlignments(alignments,bwMatrix);
                alignments.clear();
            }
            lastReadId=readId;
            q=str.indexOf('\t',p+1);
            isRc=(Integer.parseInt(str.substring(p+1,q))&RC_MASK)!=0;
            p=q; q=str.indexOf('\t',p+1);
            chrId=str.substring(p+1,q);
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
            alignments.add(new Alignment(readId,readLength,leftClipLength,readLength-rightClipLength-1,chrId,refPosOriginal,refPos,isRc,mapq));
            nRecords++;
            if (nRecords%1000==0) System.err.println("Processed "+nRecords+" records...");
            str=br.readLine();
        }
        br.close();
        if (lastReadId.length()>0) {
            printLengths(alignments,DISTANCE_THRESHOLD,tmpMap,bwLengths);
            printAlignments(alignments,bwMatrix);
            alignments.clear();
        }
        bwLengths.close(); bwMatrix.close();
    }


    private static final void printAlignments(Vector<Alignment> alignments, BufferedWriter out) throws IOException {
        final int nAlignments = alignments.size();
        int i;
        Alignment a;

        for (i=0; i<nAlignments; i++) {
            a=alignments.get(i);
            out.write(a.readId+"\t"+a.readLength+"\t"+a.readFirst+"\t"+a.readLast+"\t"+a.chrId+"\t"+a.chrFirst+"\t"+a.chrLast+"\t"+(a.isRc?"RC":"FW")+"\t"+a.mapq+"\n");
        }
    }


    /**
     * @param alignments assumed to be all and only the alignments in a contig;
     * @param tmpMap temporary space;
     * @param out format:
     * 
     * readLength,alignmentType,mapq,mostFrequentChr,mostFrequentOrientation,alignmentChr,alignmentOrientation
     * 
     * where `alignmentType` is:
     * 0: same `(chr,orientation)` combination as the most frequent alignment in 
     *    the contig, and position concordant (within `distanceThreshold`) with
     *    the previous and next alignment with the same `(chr,orientation)` 
     *    combination (if any);
     * 1: same as 0, but the position is discordant;
     * 2: same `chr` as the most frequent alignment in the contig, but different 
     *    `orientation`;
     * 3: different `chr` from the most frequent alignment in the contig.
     * 
     * Remark: only contigs whose most frequent chromosome is a standard one are
     * printed in output.
     */
    private static final void printLengths(Vector<Alignment> alignments, int distanceThreshold, HashMap<String,Integer> tmpMap, BufferedWriter out) throws IOException {
        boolean mostFrequentIsRc, alignmentIsRc;
        int i;
        int maxCount, alignmentType;
        String key, alignmentChr;
        Integer value;
        Alignment a;
        String mostFrequentChr;

        // Computing the most frequent (chr,orientation) pair.
        tmpMap.clear();
        for (i=0; i<alignments.size(); i++) {
            a=alignments.get(i);
            key=a.chrId+"@"+(a.isRc?"r":"f");
            value=tmpMap.get(key);
            if (value==null) tmpMap.put(key,1);
            else tmpMap.put(key,value+1);
        }
        maxCount=0; mostFrequentChr=""; mostFrequentIsRc=false;
        for (String tmpKey : tmpMap.keySet()) {
            value=tmpMap.get(tmpKey);
            if (value>maxCount) { maxCount=value; mostFrequentChr=tmpKey.substring(0,tmpKey.indexOf('@')); mostFrequentIsRc=tmpKey.charAt(tmpKey.length()-1)=='r'; }
        }
        if (!isStandardChromosome(mostFrequentChr)) return;

        // Printing the output
        alignments.sort(null);
        for (i=0; i<alignments.size(); i++) {
            a=alignments.get(i);
            alignmentChr=a.chrId; alignmentIsRc=a.isRc;
            if (alignmentChr.equals(mostFrequentChr)) {
                if (alignmentIsRc==mostFrequentIsRc) {
                    if ( (i==0||!alignments.get(i-1).chrId.equals(alignmentChr)||!alignments.get(i-1).isRc==alignmentIsRc?true:(Math.abs(a.readFirst-alignments.get(i-1).readLast)<=distanceThreshold)) &&
                         (i==alignments.size()-1||!alignments.get(i+1).chrId.equals(alignmentChr)||!alignments.get(i+1).isRc==alignmentIsRc?true:(Math.abs(alignments.get(i+1).readFirst-a.readLast)<=distanceThreshold))
                       ) alignmentType=0;
                    else alignmentType=1;
                }
                else alignmentType=2;
            }
            else alignmentType=3;
            out.write((alignments.get(i).readLast-alignments.get(i).readFirst)+"\t"+alignmentType+"\t"+alignments.get(i).mapq+"\t"+mostFrequentChr+"\t"+(mostFrequentIsRc?"RC":"FW")+"\t"+alignmentChr+"\t"+(alignmentIsRc?"RC":"FW")+"\n");
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
        String readId;
        String chrId;
        int readLength, readFirst, readLast, chrFirst, chrLast;
        boolean isRc;
        String mapq;

        public Alignment(String readId, int readLength, int readFirst, int readLast, String chrId, int chrFirst, int chrLast, boolean isRc, String mapq) {
            this.readId=readId; this.readLength=readLength;
            this.readFirst=isRc?readLength-readLast-1:readFirst;
            this.readLast=isRc?readLength-readFirst-1:readLast;
            this.chrId=chrId; this.chrFirst=chrFirst; this.chrLast=chrLast; this.isRc=isRc; this.mapq=mapq;
        }

        public int compareTo(Object o) {
            Alignment a = (Alignment)o;
            if (this.readFirst<a.readFirst) return -1;
            else if (this.readFirst>a.readFirst) return 1;
            return 0;
        }
    }

}