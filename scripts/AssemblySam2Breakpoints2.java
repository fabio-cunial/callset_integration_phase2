import java.util.Arrays;
import java.util.ArrayList;
import java.io.*;


/**
 * Given a QNAME-sorted assembly-to-ref SAM, the program prints points in the
 * reference where colinearity is violated by the assembly.
 */
public class AssemblySam2Breakpoints2 {
    /**
     * Max distance (on an assembled contig) between two alignments for them to
     * be considered adjacent.
     */
    private static int MAX_ADJACENCY_DISTANCE;

    /**
     * Min distance (on the reference) between two alignments (that are adjacent
     * on some contig) for them to be considered a colinearity violation.
     */
    private static int MIN_VIOLATION_DISTANCE;

    /**
     * 0=no constraint; 
     * 1=print only BND sides on standard chromosomes, but the two sides of a 
     *   BND can occur on any chromosome;
     * 2=print only BNDs where both sides belong to a standard chromosome.
     * 
     * One may want to print only violations involving standard chromosomes, 
     * since non-standard chrs are highly repetitive and BNDs between them and 
     * standard chrs might not represent new adjacencies (since such chrs are
     * unplaced).
     */
    private static int CHROMOSOME_MODE;

    /**
     * Consider as breakpoints also the projections onto ref. of the first/last
     * positions on the contig of the longest chain in a connected component.
     */
    private static boolean PRINT_CHAIN_START_END;

    /**
     * Slack for deciding whether a component is contained in another one in an
     * assembled contig.
     */
    private static int CONTAINMENT_SLACK_BP;

    /**
     * Min length of a CIGAR INS/DEL for it to create breakpoints
     */
    private static int MIN_INTERNAL_SV_LENGTH;

    /**
     * Temporary, reused space.
     */
    private static long[] minViolations, maxLength;
    private static Alignment[] minAlignment;
    private static ArrayList<Interval> componentIntervals = new ArrayList<Interval>();  // Forward orientation in the contig
    private static int[] nBreakpoints = new int[6];  // total, violation1, violation2, violation3, first/last, internal

    
    /**
     * Remark: the program prints a CSV of breakpoints, where every row 
     * represents a single side of an adjacency, with format:
     * 
     * CHR,POS, MATE_CHR,MATE_POS
     * 
     * Rows with MATE_CHR=-1, MATE_POS=-1 represent breakpoints in the reference
     * that connect to new sequence in the assembly that is not covered by any
     * alignment with the reference. Every row with a mate has a symmetrical
     * row. The direction of the adjacencies is not reported.
     * 
     * Remark: the file is not necessarily sorted and it may contain multiple 
     * rows with identical or similar CHR,POS values. This is fine, since users
     * downstream will match such true points to a VCF.
     * 
     * @param args
     * 0: all alignments with the same QNAME (i.e. assembly contig ID) are 
     *    assumed to form a contiguous block (e.g. `samtools sort -@ 4 -n -O SAM
     *    hap1.bam`).
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_SAM = args[0];
        MAX_ADJACENCY_DISTANCE=Integer.parseInt(args[1]);
        MIN_VIOLATION_DISTANCE=Integer.parseInt(args[2]);
        CHROMOSOME_MODE=Integer.parseInt(args[3]);
        PRINT_CHAIN_START_END=Integer.parseInt(args[4])==1;
        CONTAINMENT_SLACK_BP=Integer.parseInt(args[5]);
        MIN_INTERNAL_SV_LENGTH=Integer.parseInt(args[6]);

        final boolean DEBUG = false;
        
        boolean isRc;
        char c;
        int i, j, p, q;
        int refPos, refPosOriginal, readPos, cigarLength, leftClipLength, rightClipLength, readLength, length, nRecords, flag, readPosLast, mapq;
        String str, readId, chrId, lastReadId, cigar;
        BufferedReader br;
        String[] tokens;
        ArrayList<Alignment> alignments;
        ArrayList<Integer> internalBreakpoints;  // Forward ref orientation
        
        alignments = new ArrayList<Alignment>();
        internalBreakpoints = new ArrayList<Integer>();
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_SAM)));
        str=br.readLine(); nRecords=0; lastReadId="";
        while (str!=null) {
            if (str.charAt(0)=='@') { str=br.readLine(); continue; }
            p=str.indexOf('\t');
            readId=str.substring(0,p);
            if (lastReadId.length()!=0 && !readId.equals(lastReadId)) {
                if (DEBUG) printNearestNeighborDistances(alignments);
                printBreakpoints(alignments);
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
            i=0; leftClipLength=0; rightClipLength=0; readPosLast=-1; internalBreakpoints.clear();
            for (j=1; j<cigarLength; j++) {
                c=cigar.charAt(j);
                if (c=='M' || c=='=' || c=='X') {
                    length=Integer.parseInt(cigar.substring(i,j));
                    refPos+=length; readPos+=length;
                    i=j+1;
                }
                else if (c=='I') {
                    length=Integer.parseInt(cigar.substring(i,j));
                    if (length>=MIN_INTERNAL_SV_LENGTH) {
                        internalBreakpoints.add(refPos);
                        internalBreakpoints.add(refPos);
                    }
                    readPos+=length;
                    i=j+1;
                }
                else if (c=='D') {
                    length=Integer.parseInt(cigar.substring(i,j));
                    if (length>=MIN_INTERNAL_SV_LENGTH) {
                        internalBreakpoints.add(refPos);
                        internalBreakpoints.add(refPos+length+1);
                    }
                    refPos+=length;
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
            alignments.add(new Alignment(chrId,refPosOriginal,refPos,readId,leftClipLength,readPosLast,isRc,readLength,mapq,internalBreakpoints));
            nRecords++;
            if (nRecords%1000==0) System.err.println("Processed "+nRecords+" records...");
            str=br.readLine();
        }
        br.close();
        if (lastReadId.length()!=0) {
            if (DEBUG) printNearestNeighborDistances(alignments);
            printBreakpoints(alignments);
            alignments.clear();
        }
        System.err.println("Total breakpoints created: "+nBreakpoints[0]);
        System.err.println("Colinearity violation 1 (chr): "+nBreakpoints[1]);
        System.err.println("Colinearity violation 2 (orientation): "+nBreakpoints[2]);
        System.err.println("Colinearity violation 3 (distance): "+nBreakpoints[3]);
        System.err.println("First/last of chain: "+nBreakpoints[4]);
        System.err.println("Internal breakpoints: "+nBreakpoints[5]);
    }


    /**
     * Prints all the breakpoints of an assembled contig.
     * 
     * Remark: regions of an assembled contig with no alignments induce
     * connected components in the alignments DAG. One could mark the boundaries
     * of such regions as breakpoints in the reference (since they connect the
     * ref with new sequence), but of course they are not new adjacencies 
     * between existing regions of the ref.
     * 
     * Remark: we compute a longest chain in every connected component, to try
     * to explain as much of the assembled contig as possible using the ref. 
     * Among all longest chains, we pick one with smallest number of violations,
     * to try to be as compatible as possible with the ref.
     * 
     * @param alignments assumed to be all and only the alignments in an 
     * assembled contig.
     */
    private static final void printBreakpoints(ArrayList<Alignment> alignments) throws IOException {
        final int N_ALIGNMENTS = alignments.size();

        int i, j;
        int nComponents, nComponentsWithViolations, nContainedComponents, violationType, nBreakpointsBefore;
        long aLength, bLength, nViolations;
        Alignment a, b;
        Interval tmpInterval;

        nBreakpointsBefore=nBreakpoints[0];

        // Computing connected components and their longest chains
        nComponents=0;
        Alignment.ORDER=1; alignments.sort(null);  // Forward read order
        for (i=0; i<N_ALIGNMENTS; i++) {
            a=alignments.get(i);
            if (a.connectedComponent==-1) a.connectedComponent=nComponents++;
            aLength=a.readLast-a.readFirst+1;
            for (j=i+1; j<N_ALIGNMENTS; j++) {
                b=alignments.get(j);
                if (b.readFirst>a.readLast+MAX_ADJACENCY_DISTANCE) break;
                if (b.readFirst<a.readLast-MAX_ADJACENCY_DISTANCE) continue;
                b.connectedComponent=a.connectedComponent;
                a.hasChildren=true;
                bLength=a.maxLength+aLength;
                if (bLength>b.maxLength) { b.maxLength=bLength; b.parent=a; }
            }
        }
        if (maxLength==null || maxLength.length<nComponents) maxLength = new long[nComponents];
        Arrays.fill(maxLength,0,nComponents,0);
        for (i=0; i<N_ALIGNMENTS; i++) {
            a=alignments.get(i);
            aLength=a.maxLength+a.readLast-a.readFirst+1;
            if (aLength>maxLength[a.connectedComponent]) maxLength[a.connectedComponent]=aLength;
        }
        System.err.println("Contig "+alignments.get(0).readId+" has "+nComponents+" connected components, whose longest chains have length (bp):");
        for (i=0; i<nComponents; i++) System.err.println(maxLength[i]+"");


        if (alignments.get(0).readId.equals("HG00097#2#CM094075.1")) {
            for (i=0; i<N_ALIGNMENTS; i++) {
                a=alignments.get(i);
                System.err.println("Alignment "+i+": component="+a.connectedComponent+", "+a.chrId+", isRc="+a.isRc+", chrFirst="+a.chrFirst+" chrLast="+a.chrLast+", read "+a.readId+", readFirst="+a.readFirst+" readLast="+a.readLast);
            }
        }


        // Finding a longest chain with smallest number of violations, for every 
        // connected component.
        if (minViolations==null || minViolations.length<nComponents) minViolations = new long[nComponents]; 
        if (minAlignment==null || minAlignment.length<nComponents) minAlignment = new Alignment[nComponents];
        Arrays.fill(minViolations,0,nComponents,Long.MAX_VALUE);
        Arrays.fill(minAlignment,0,nComponents,null);
        for (i=0; i<N_ALIGNMENTS; i++) {
            a=alignments.get(i);
            aLength=a.maxLength+a.readLast-a.readFirst+1;
            if (aLength!=maxLength[a.connectedComponent]) continue;
            b=a; nViolations=0;
            while (b.parent!=null) {
                if (getViolationType(b,b.parent)>0) nViolations++;
                b=b.parent;
            }
            if (nViolations<minViolations[a.connectedComponent]) { minViolations[a.connectedComponent]=nViolations; minAlignment[a.connectedComponent]=a; }
        }
        nComponentsWithViolations=0;
        for (i=0; i<nComponents; i++) {
            if (minViolations[i]>0) {
                nComponentsWithViolations++;
                System.err.println("Contig "+alignments.get(0).readId+", component "+i+" has longest chain with "+minViolations[i]+" violations");
            }
        }
        System.err.println("Contig "+alignments.get(0).readId+" has "+nComponentsWithViolations+" components with violations");

        // Detecting chains that are contained in other chains on the contig
        componentIntervals.clear();
        for (i=0; i<nComponents; i++) {
            a=minAlignment[i];
            tmpInterval = new Interval(-1,a.readLast,i);
            while (a.parent!=null) a=a.parent;
            tmpInterval.first=a.readFirst;
            tmpInterval.isContained=false;
            componentIntervals.add(tmpInterval);
        }
        Interval.ORDER=0; componentIntervals.sort(null);
        for (i=0; i<nComponents; i++) {
            if (componentIntervals.get(i).isContained) continue;
            for (j=i+1; j<nComponents; j++) {
                if (componentIntervals.get(j).last>componentIntervals.get(i).last+CONTAINMENT_SLACK_BP) break;
                componentIntervals.get(j).isContained=true;
            }
            for (j=i-1; j>=0; j--) {
                if (componentIntervals.get(j).first<componentIntervals.get(i).first-CONTAINMENT_SLACK_BP) break;
                if (componentIntervals.get(j).last<=componentIntervals.get(i).last+CONTAINMENT_SLACK_BP) componentIntervals.get(j).isContained=true;
            }
        }
        Interval.ORDER=1; componentIntervals.sort(null);
        nContainedComponents=0;
        for (i=0; i<nComponents; i++) { if (componentIntervals.get(i).isContained) nContainedComponents++; }
        System.err.println("Contig "+alignments.get(0).readId+" has "+nContainedComponents+" components contained in other components");

        // Printing breakpoints using only connected components that are not 
        // contained in other connected components on the contig.
        for (i=0; i<nComponents; i++) {
            if (componentIntervals.get(i).isContained) continue;
            a=minAlignment[i];
            if (PRINT_CHAIN_START_END) printFirstLast(a,true);
            a.printInternalBreakpoints();
            while (a.parent!=null) {
                violationType=getViolationType(a,a.parent);
                if (violationType>0) printViolation(a,a.parent,violationType);
                a=a.parent;
                a.printInternalBreakpoints();
            }
            if (PRINT_CHAIN_START_END) printFirstLast(a,false);
        }
        System.err.println("Contig "+alignments.get(0).readId+" induced "+(nBreakpoints[0]-nBreakpointsBefore)+" breakpoints");
    }


    /**
     * @return 0 iff the adjacency `parent->alignment` does not violate 
     * colinearity WRT the reference; >0 otherwise.
     */
    private static final int getViolationType(Alignment alignment, Alignment parent) {
        if (!parent.chrId.equals(alignment.chrId)) return 1;
        if (parent.isRc!=alignment.isRc) return 2;
        if ((!alignment.isRc && alignment.chrFirst>parent.chrLast+MIN_VIOLATION_DISTANCE) || (alignment.isRc && parent.chrFirst>alignment.chrLast+MIN_VIOLATION_DISTANCE)) return 3;
        return 0;
    }


    private static final void printViolation(Alignment alignment, Alignment parent, int violationType) {
        final boolean isStandardChrAlignment = isStandardChromosome(alignment.chrId);
        final boolean isStandardChrParent = isStandardChromosome(parent.chrId);
        final boolean printAlignment = CHROMOSOME_MODE==0 || (CHROMOSOME_MODE==1 && isStandardChrAlignment) || (CHROMOSOME_MODE==2 && isStandardChrAlignment && isStandardChrParent);
        final boolean printParent = CHROMOSOME_MODE==0 || (CHROMOSOME_MODE==1 && isStandardChrParent) || (CHROMOSOME_MODE==2 && isStandardChrAlignment && isStandardChrParent);

        if (printParent) {
            nBreakpoints[0]++; nBreakpoints[violationType]++;
            System.out.println(parent.chrId+","+(parent.isRc?parent.chrFirst:parent.chrLast)+","+alignment.chrId+","+(alignment.isRc?alignment.chrLast:alignment.chrFirst));
        }
        if (printAlignment) {
            nBreakpoints[0]++; nBreakpoints[violationType]++;
            System.out.println(alignment.chrId+","+(alignment.isRc?alignment.chrLast:alignment.chrFirst)+","+parent.chrId+","+(parent.isRc?parent.chrFirst:parent.chrLast));
        }

        /* Old code for printing a VCF:
        else if (OUTPUT_MODE==1 && printParent && printAlignment) {
            nBreakpoints[0]++; nBreakpoints[violationType]++;
            final String refChrom = parent.chrId;
            final int refPos = (parent.isRc?parent.chrFirst:parent.chrLast)+1;
            final String altChrom = alignment.chrId;
            final int altPos = (alignment.isRc?alignment.chrLast:alignment.chrFirst)+1;
            String alt;
            if (!alignment.isRc) alt="["+altChrom+":"+altPos+"[";
            else alt="]"+altChrom+":"+altPos+"]";
            alt=parent.isRc?alt+"N":"N"+alt;
            System.out.println(refChrom+"\t"+refPos+"\t"+nBreakpoints[0]+"\tN\t"+alt+"\t"+Math.min(alignment.mapq,parent.mapq)+"\tPASS\tSVTYPE=BND\tGT\t0/1");
        }
        */
    }


    /**
     * @param last prints the projection on the ref. of the last (TRUE) of first
     * (FALSE) position of `alignment` in contig order.
     */
    private static void printFirstLast(Alignment alignment, boolean last) {
        if (CHROMOSOME_MODE==0 || (CHROMOSOME_MODE==1 && isStandardChromosome(alignment.chrId))) {
            nBreakpoints[0]++; nBreakpoints[4]++;
            if (last) System.out.println(alignment.chrId+","+(alignment.isRc?alignment.chrFirst:alignment.chrLast)+",-1,-1");
            else System.out.println(alignment.chrId+","+(alignment.isRc?alignment.chrLast:alignment.chrFirst)+",-1,-1");
        }
    }


    /**
     * Prints the distance, along an assembled contig, between every alignment
     * and the closest neighbor to its right that starts after its last pos.
     */
    private static final void printNearestNeighborDistances(ArrayList<Alignment> alignments) {
        final int N_ALIGNMENTS = alignments.size();

        int i, j;
        Alignment a, b;

        Alignment.ORDER=1; alignments.sort(null);  // Forward read order
        System.err.println("Contig "+alignments.get(0).readId+", right nearest neighbor distances:");
        for (i=0; i<N_ALIGNMENTS; i++) {
            a=alignments.get(i);
            for (j=i+1; j<N_ALIGNMENTS; j++) {
                b=alignments.get(j);
                if (b.readFirst>a.readLast) { System.err.println((b.readFirst-a.readLast)+""); break; }
            }
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
        /**
         * 0: position on the chromosome; every alignment is assumed to use the
         *    same chromosome;
         * 1: position on the read; every alignment is assumed to use the same
         *    read.
         */
        public static int ORDER = 1;

        /**
         * Properties of the alignment
         */
        public boolean isRc;
        public String chrId, readId;
        public int chrFirst, chrLast;  // 0-based, inclusive, forward orientation in the reference.
        public int readFirst, readLast;  // 0-based, inclusive, forward orientation in the read.
        public int readLength, mapq;

        /**
         * Max length of a sequence of chainable alignments that leads up to
         * the beginning of this alignment (and excludes this alignment).
         */
        public long maxLength;
        
        public Alignment parent;
        public boolean hasChildren;
        public int connectedComponent;

        /**
         * Format: `refPosFrom1,refPosTo1, ..., refPosFromN,refPosToN`.
         */
        public ArrayList<Integer> internalBreakpoints;


        /**
         * @param chrFirst,chrLast assumed to be in forward orientation in the
         * reference;
         * @param readFirst,readLast in forward orientation in the read if 
         * `isRc=false`, in reverse orientation if `isRc=true`.
         */
        public Alignment(String chrId, int chrFirst, int chrLast, String readId, int readFirst, int readLast, boolean isRc, int readLength, int mapq, ArrayList<Integer> ib) {
            this.chrId = chrId; this.chrFirst=chrFirst; this.chrLast=chrLast;
            this.readId = readId; this.isRc=isRc; this.readLength=readLength;
            this.readFirst=isRc?readLength-readLast-1:readFirst;
            this.readLast=isRc?readLength-readFirst-1:readLast;
            this.mapq=mapq;

            maxLength=0; parent=null; hasChildren=false; connectedComponent=-1;
            internalBreakpoints = new ArrayList<Integer>(ib);
        }


        /**
         * Remark: INS breakpoints are printed even though we don't know
         * whether they connect to new sequence, or to another region of the 
         * reference.
         */
        public void printInternalBreakpoints() {
            if (CHROMOSOME_MODE!=0 && !isStandardChromosome(chrId)) return;

            int i;
            int from, to;
            final int length = internalBreakpoints.size();

            for (i=0; i<length; i+=2) {
                from=internalBreakpoints.get(i);
                to=internalBreakpoints.get(i+1);
                if (to==from) {
                    System.out.println(chrId+","+from+",-1,-1");
                    nBreakpoints[0]++; nBreakpoints[5]++;
                }
                else {
                    System.out.println(chrId+","+from+","+chrId+","+to);
                    System.out.println(chrId+","+to+","+chrId+","+from);
                    nBreakpoints[0]+=2; nBreakpoints[5]+=2;
                }
            }
        }


        public int compareTo(Object o) {
            Alignment a = (Alignment)o;
            if (ORDER==0) {
                // chrId is not used since assumed equal
                if (this.chrFirst<a.chrFirst) return -1;
                else if (this.chrFirst>a.chrFirst) return 1;
                if (this.chrLast<a.chrLast) return -1;
                else if (this.chrLast>a.chrLast) return 1;
            } else {
                // readId is not used since assumed equal
                if (this.readFirst<a.readFirst) return -1;
                else if (this.readFirst>a.readFirst) return 1;
                if (this.readLast<a.readLast) return -1;
                else if (this.readLast>a.readLast) return 1;
            }
            return 0;
        }
    }


    private static class Interval implements Comparable {
        private static int ORDER = 0;  // 0=first, 1=id

        public boolean isContained;
        public int id;
        public int first, last;  // Forward orientation in the contig

        public Interval(int first, int last, int id) { 
            this.first=first; this.last=last; this.id=id;
            isContained=false;
        }

        public int compareTo(Object o) {
            Interval other = (Interval)o;
            if (ORDER==0) {
                if (this.first<other.first) return -1;
                else if (this.first>other.first) return 1;
                if (this.last<other.last) return -1;
                else if (this.last>other.last) return 1;
            } else {
                if (this.id<other.id) return -1;
                else if (this.id>other.id) return 1;
            }
            return 0;
        }
    }

}