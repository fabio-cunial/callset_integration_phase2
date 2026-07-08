import java.util.Arrays;
import java.util.ArrayList;
import java.io.*;


/**
 * Given a QNAME-sorted assembly-to-ref SAM, the program prints points in the
 * reference where colinearity is violated by the assembly.
 * 
 * Remark: only violations involving standard chromosomes are printed.
 */
public class AssemblySam2Breakpoints2 {
    
    /**
     * @param args
     * 0: all alignments with the same QNAME (i.e. assembly contig ID) are 
     *    assumed to form a contiguous block (e.g. `samtools sort -@ 8 -m 3G -n 
     *    -O SAM hap1.bam`);
     * 1: max distance (on an assembled contig) between two alignments for them
     *    to be considered adjacent;
     * 2: min distance (on the reference) between two alignments (that are 
     *    adjacent on some contig) for them to be considered a colinearity 
     *    violation;
     * 3: 0=a CSV of points, with format `CHR,POS,ID`, not necessarily sorted;
     *    1=a BND VCF, not necessarily sorted or canonized.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_SAM = args[0];
        final int MAX_ADJACENCY_DISTANCE = Integer.parseInt(args[1]);
        final int MIN_VIOLATION_DISTANCE = Integer.parseInt(args[2]);
        final int OUTPUT_MODE = Integer.parseInt(args[3]);
        
        boolean isRc;
        char c;
        int i, j, p, q;
        int refPos, refPosOriginal, readPos, cigarLength, leftClipLength, rightClipLength, readLength, length, nRecords, flag, readPosLast, mapq;
        String str, readId, chrId, lastReadId, cigar;
        BufferedReader br;
        String[] tokens;
        ArrayList<Alignment> alignments;
        
        alignments = new ArrayList<Alignment>();
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_SAM)));
        str=br.readLine(); nRecords=0; lastReadId="";
        while (str!=null) {
            if (str.charAt(0)=='@') { str=br.readLine(); continue; }
            p=str.indexOf('\t');
            readId=str.substring(0,p);
            if (lastReadId.length()!=0 && !readId.equals(lastReadId)) {
                printBreakpoints(alignments,MAX_ADJACENCY_DISTANCE,MIN_VIOLATION_DISTANCE,OUTPUT_MODE);
                alignments.clear();
            }
            lastReadId=readId;
            q=str.indexOf('\t',p+1);
            flag=Integer.parseInt(str.substring(p+1,q));
            isRc=(flag&16)!=0;
            p=q; q=str.indexOf('\t',p+1);
            chrId=str.substring(p+1,q);
            p=q; q=str.indexOf('\t',p+1);
            refPos=Integer.parseInt(str.substring(p+1,q))-1;  // 0-based, first matching pos.
            refPosOriginal=refPos;
            refPos--;  // Moving the pointer right before the first matching pos.
            p=q; q=str.indexOf('\t',p+1);
            mapq=Integer.parseInt(str.substring(p+1,q));
            p=q; q=str.indexOf('\t',p+1);
            cigar=str.substring(p+1,q); cigarLength=cigar.length();
            readPos=-1;  // 0-based
            i=0; leftClipLength=0; rightClipLength=0; readPosLast=-1;
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
                    else if (j==cigarLength-1) {
                        readPosLast=readPos;
                        rightClipLength=length;
                    }
                    readPos+=length;
                    i=j+1;
                }
                else if (c=='P') i=j+1;
            }
            if (readPosLast==-1) readPosLast=readPos;
            readLength=readPos+1;
            alignments.add(new Alignment(chrId,refPosOriginal,refPos,readId,leftClipLength,readPosLast,isRc,readLength,mapq));
            nRecords++;
            if (nRecords%1000==0) System.err.println("Processed "+nRecords+" records...");
            str=br.readLine();
        }
        br.close();
        if (lastReadId.length()!=0) {
            printBreakpoints(alignments,MAX_ADJACENCY_DISTANCE,MIN_VIOLATION_DISTANCE,OUTPUT_MODE);
            alignments.clear();
        }
    }


    /**
     * Remark: only violations involving standard chromosomes are printed.
     * 
     * Remark: gaps in an assembled contig induce connected components in the 
     * alignments DAG. One could mark the boundaries of such gaps as breakpoints
     * in the reference, but they are not new adjacencies so we don't do it
     * here.
     * 
     * Remark: we compute a longest chain in every connected component to try to
     * explain as much of the assembled contig as possible. Among all longest 
     * chains we pick the one with smallest number of violations, to try to be
     * as compatible as possible with the reference.
     * 
     * @param alignments assumed to be all and only the alignments in an 
     * assembled contig.
     */
    private static final void printBreakpoints(ArrayList<Alignment> alignments, int adjacencyThreshold, int violationThreshold, int outputMode) throws IOException {
        final int N_ALIGNMENTS = alignments.size();

        int i, j;
        int idGenerator, nComponents;
        long aLength, bLength, nViolations;
        Alignment a, b;
        long[] minViolations;
        long[] maxLength;
        Alignment[] minAlignment;

        // Computing a longest chain for every connected component
        nComponents=0;
        Alignment.ORDER=1; alignments.sort(null);  // Read order
        for (i=0; i<N_ALIGNMENTS; i++) {
            a=alignments.get(i);
            if (a.connectedComponent==-1) a.connectedComponent=nComponents++;
            aLength=a.readLast-a.readFirst+1;
            for (j=i+1; j<N_ALIGNMENTS; j++) {
                b=alignments.get(j);
                if (b.readFirst>a.readLast+adjacencyThreshold) break;
                if (b.readFirst<a.readLast-adjacencyThreshold) continue;
                b.connectedComponent=a.connectedComponent;
                a.hasChildren=true;
                bLength=a.maxLength+aLength;
                if (bLength>b.maxLength) { b.maxLength=bLength; b.parent=a; }
            }
        }
        maxLength = new long[nComponents];
        for (i=0; i<N_ALIGNMENTS; i++) {
            a=alignments.get(i);
            aLength=a.maxLength+a.readLast-a.readFirst+1;
            if (aLength>maxLength[a.connectedComponent]) maxLength[a.connectedComponent]=aLength;
        }
        for (i=0; i<nComponents; i++) System.err.println("Connected component "+i+" has longest chain of length "+maxLength[i]+"bp");

        // Finding a longest chain with smallest number of violations, for every 
        // connected component.
        minAlignment = new Alignment[nComponents];
        minViolations = new long[nComponents]; 
        Arrays.fill(minViolations,Long.MAX_VALUE);
        for (i=0; i<N_ALIGNMENTS; i++) {
            a=alignments.get(i);
            aLength=a.maxLength+a.readLast-a.readFirst+1;
            if (aLength!=maxLength[a.connectedComponent]) continue;
            b=a; nViolations=0;
            while (b.parent!=null) {
                if (isViolation(b,b.parent,violationThreshold)) nViolations++;
                b=b.parent;
            }
            if (nViolations<minViolations[a.connectedComponent]) { minViolations[a.connectedComponent]=nViolations; minAlignment[a.connectedComponent]=a; }
        }
        for (i=0; i<nComponents; i++) System.err.println("Connected component "+i+" has longest chain with "+minViolations[i]+" violations");

        // Printing breakpoints for every connected component
        idGenerator=0;
        for (i=0; i<nComponents; i++) {
            if (minViolations[i]==0) continue;
            a=minAlignment[i];
            while (a.parent!=null) {
                if (isStandardChromosome(a.chrId) && isStandardChromosome(a.parent.chrId) && isViolation(a,a.parent,violationThreshold)) printViolation(a,a.parent,idGenerator++,outputMode);
                a=a.parent;
            }
        }
        System.err.println("Created "+idGenerator+" breakpoints");
    }


    /**
     * @return TRUE iff the adjacency `parent->alignment` violates colinearity
     * WRT the reference.
     */
    private static final boolean isViolation(Alignment alignment, Alignment parent, int distanceThreshold) {
        return !parent.chrId.equals(alignment.chrId) || parent.isRc!=alignment.isRc || (!alignment.isRc && alignment.chrFirst>parent.chrLast+distanceThreshold) || (alignment.isRc && parent.chrFirst>alignment.chrLast+distanceThreshold);
    }


    /**
     * @param printMode
     * 0: a CSV of points, with format `CHR,POS,ID`; not necessarily sorted;
     * 1: a BND VCF; not necessarily sorted or canonized.
     */
    private static final void printViolation(Alignment alignment, Alignment parent, int id, int printMode) {
        if (printMode==0) {
            System.out.println(parent.chrId+","+(parent.isRc?parent.chrFirst:parent.chrLast)+","+id);
            System.out.println(alignment.chrId+","+(alignment.isRc?alignment.chrLast:alignment.chrFirst)+","+id);
        }
        else if (printMode==1) {
            final String refChrom = parent.chrId;
            final int refPos = (parent.isRc?parent.chrFirst:parent.chrLast)+1;
            final String altChrom = alignment.chrId;
            final int altPos = (alignment.isRc?alignment.chrLast:alignment.chrFirst)+1;
            String alt;
            if (!alignment.isRc) alt="["+altChrom+":"+altPos+"[";
            else alt="]"+altChrom+":"+altPos+"]";
            alt=parent.isRc?alt+"N":"N"+alt;
            System.out.println(refChrom+"\t"+refPos+"\t"+id+"\tN\t"+alt+"\t"+Math.min(alignment.mapq,parent.mapq)+"\tPASS\tSVTYPE=BND\tGT\t0/1");
        }
    }


    /**
     * Designed to be fast albeit not fully correct
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
        public static int ORDER = 0;  // 0=chr, 1=read.

        /**
         * Properties of the alignment
         */
        public boolean isRc;
        public String chrId, readId;
        public int chrFirst, chrLast, readFirst, readLast;  // 0-based, inclusive.
        public int readLength, mapq;

        /**
         * Max length of a sequence of chainable alignments that leads up to
         * the beginning of this alignment.
         */
        public long maxLength;
        public boolean hasChildren;
        public Alignment parent;
        public int connectedComponent;


        public Alignment(String chrId, int chrFirst, int chrLast, String readId, int readFirst, int readLast, boolean isRc, int readLength, int mapq) {
            this.chrId = chrId; this.chrFirst=chrFirst; this.chrLast=chrLast;
            this.readId = readId; this.readFirst=readFirst; this.readLast=readLast;
            this.readLength=readLength; this.mapq=mapq;
            this.isRc=isRc;
            maxLength=0; parent=null; hasChildren=false; connectedComponent=-1;
        }

        public int compareTo(Object o) {
            Alignment a = (Alignment)o;
            if (ORDER==0) {
                if (this.chrFirst<a.chrFirst) return -1;
                else if (this.chrFirst>a.chrFirst) return 1;
                if (this.chrLast<a.chrLast) return -1;
                else if (this.chrLast>a.chrLast) return 1;
            } else {
                if (this.readFirst<a.readFirst) return -1;
                else if (this.readFirst>a.readFirst) return 1;
                if (this.readLast<a.readLast) return -1;
                else if (this.readLast>a.readLast) return 1;
            }
            return 0;
        }
    }

}