import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * For each long call, computes the max R^2 coefficient with the short calls 
 * located within a max distance from the long's endpoints and not overlapping 
 * the long.
 * 
 * Remark: both the short and the long categories may include INS, DEL and 
 * replacements. Every short that affects bases within the max distance and that
 * does not overlap the long contributes to the long's R^2.
 * 
 * Remark: for speed, the program works on a small TSV projection of the
 * original VCF, where every record comes from the same chromosome, samples are 
 * indexes in the list of all samples in the cohort (rather than sample IDs),
 * every record contains only the samples with GT=ALT, in increasing order, and
 * records are sorted by POS.
 */
public class Rsquare {
    /**
     * Variant types
     */
    private static final int TYPE_SNP = 0;
    private static final int TYPE_DEL = 1;
    private static final int TYPE_INS = 2;
    private static final int TYPE_SUB = 3;

    /**
     * Array size parameters
     */
    private static final int CAPACITY = 100;  // Arbitrary
    private static final double RESIZE_FACTOR = 1.5;  // Arbitrary

    /**
     * Format of each row:
     * 
     * pos,first,last,type,length,  sample1,count1, sample2,counts2, ..., sampleK,countK
     * 
     * where:
     * - `pos` is the one-based POS from the VCF (stored just for printing
     *   the output);
     * - `[first..last]` is the interval of zero-based, inclusive reference 
     *   positions that are affected by the variant;
     * - `countX` is in {0,1,2}.
     * 
     * Remark: rows may have different lengths.
     */
    private static final int FIRST_SAMPLE_INDEX = 5;  // First index of a sample in the `active*` tables.
    private static String chrId;
    private static int[][] activeLong, activeShort;
    private static int[] activeLongLast, activeShortLast;
    private static int firstActiveLong, firstActiveShort, lastActiveLong, lastActiveShort;
    private static int totalNSamples, minAC;
    private static String[] activeShortIds;  // The ID of every active short
    private static String[][] activeLongIds;  // The ID of every active long, and the ID of its short with max R^2.

    /**
     * Format of each row:
     * 
     * tmpAvg, tmpD, R^2, AC, N, ACs, Ns
     * 
     * where:
     * - `tmp*` are just temporary variables that are cached for speeding up
     *   computation;
     * - `AC` is the total number of ALT alleles;
     * - `N` is the number of samples with GT=ALT;
     * - `ACs` is the AC of the short call that maximizes R^2;
     * - `Ns` is the N of the short call that maximizes R^2.
     */
    private static double[][] rSquaredLong;

    /**
     * Format of each row: 
     * 
     * tmpAvg, tmpD, AC, N
     */
    private static double[][] rSquaredShort;

    /**
     * Scratch space
     */
    private static int[] tmpArray1;


    /**
     * @param args 0 format: 
     * 
     * POS, refLength, altLength, ID,  sample=gtCount, ...
     * 
     * where `POS` is sorted and comes from the VCF, and `gtCount` is an integer
     * in {0,1,2}; all the records are expected to come from the same chrom;
     * @param args 1 the chromosome all records come from; must be one of the
     * autosomes (no chrX, chrY, chrM);
     * @param args 2 total number of samples in the entire cohort;
     * @param args 3 only variants of length >=this are kept as long;
     * @param args 4 only variants of length <=this are kept as short;
     * @param args 5 max bp distance to compare two variants;
     * @param args 6 min AF for a long or short variant to be considered; long
     * variants below this AF are not printed in output;
     * @param args 7 output BED with format: 
     * 
     * CHROM, P, P+1, ID, type, length, AC, N, R^2, ACs, Ns, IDs
     * 
     * where `P` is the one-based POS from the VCF minus one, `N` is the number 
     * of samples with GT=ALT and R^2=-1 iff the long cannot be compared to any
     * short.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV_GZ = args[0];
        chrId=args[1];
        totalNSamples=Integer.parseInt(args[2]);
        final int MIN_LONG_LENGTH = Integer.parseInt(args[3]);
        final int MAX_SHORT_LENGTH = Integer.parseInt(args[4]);
        final int MAX_DISTANCE_BP = Integer.parseInt(args[5]);
        final double MIN_AF = Double.parseDouble(args[6]);
        final String OUTPUT_BED = args[7];

        minAC=(int)Math.ceil(2*totalNSamples*MIN_AF);
        final int QUANTUM = 1000;  // Arbitrary

        int p1, p2, p3, p4;
        int pos, length, refLength, altLength, type;
        long nRecords, nLong, nShort, nLongWithAF, nShortWithAF;
        String str, id;
        BufferedReader br;
        BufferedWriter bw;
        
        // Allocating reused space
        activeLong = new int[CAPACITY][FIRST_SAMPLE_INDEX+2*CAPACITY];
        activeLongLast = new int[CAPACITY];
        firstActiveLong=-1; lastActiveLong=-1;
        activeShort = new int[CAPACITY][FIRST_SAMPLE_INDEX+2*CAPACITY];
        activeShortLast = new int[CAPACITY];
        activeLongIds = new String[CAPACITY][2];
        activeShortIds = new String[CAPACITY];
        firstActiveShort=-1; lastActiveShort=-1;
        tmpArray1 = new int[FIRST_SAMPLE_INDEX+2*totalNSamples];
        rSquaredLong = new double[CAPACITY][7];
        rSquaredShort = new double[CAPACITY][4];

        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_TSV_GZ))));
        bw = new BufferedWriter(new FileWriter(OUTPUT_BED));
        str=br.readLine(); nRecords=0; nLong=0; nShort=0; nLongWithAF=0; nShortWithAF=0; 
        while (str!=null) {
            p1=str.indexOf('\t');
            pos=Integer.parseInt(str.substring(0,p1));
            p2=str.indexOf('\t',p1+1);
            refLength=Integer.parseInt(str.substring(p1+1,p2));
            p3=str.indexOf('\t',p2+1);
            altLength=Integer.parseInt(str.substring(p2+1,p3));
            p4=str.indexOf('\t',p3+1);
            id=str.substring(p3+1,p4);
            type=getType(refLength,altLength);
            length=getLength(refLength,altLength,type);
            if (length>=MIN_LONG_LENGTH) {
                nLong++;
                removeInactiveShorts(pos-1-MAX_DISTANCE_BP-1);
                removeInactiveLongs(pos-1-MAX_DISTANCE_BP-1,bw);
                if (addActiveLong(id,pos,refLength,altLength,type,length,str,p4+1)) {
                    nLongWithAF++;
                    initializeMaxRsquare(MAX_DISTANCE_BP);
                }
            }
            else if (length<=MAX_SHORT_LENGTH) {
                nShort++;
                removeInactiveShorts(pos-1-MAX_DISTANCE_BP-1);
                removeInactiveLongs(pos-1-MAX_DISTANCE_BP-1,bw);
                if (addActiveShort(id,pos,refLength,altLength,type,length,str,p4+1)) {
                    nShortWithAF++;
                    updateMaxRsquare(MAX_DISTANCE_BP);
                }
            }
            else {
                // Other lengths are completely discarded
            }

            // Next iteration
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            str=br.readLine();
        }
        br.close();
        printRemainingLongs(bw);
        System.err.println("Processed "+nRecords+" total records, "+nLong+" long ("+nLongWithAF+" with the desired AF), "+nShort+" short ("+nShortWithAF+" with the desired AF).");
    }


    /**
     * @param pos one-based, from the VCF;
     * @param refLength,altLength from the VCF;
     * @param p the starting position of `str` from which to begin parsing
     * (inclusive);
     * @return FALSE iff the number of occurrences of the record is <minAC.
     */
    private static final boolean addActiveLong(String id, int pos, int refLength, int altLength, int type, int length, String str, int p) {
        int i;
        int last, nSamplesLong;
        double avg, denom;

        // Loading the record and filtering by AC
        tmpArray1[0]=pos;
        getInterval(pos,refLength,altLength,type,tmpArray1,1);
        tmpArray1[3]=type;
        tmpArray1[4]=length;
        last=loadSamples(str,p,tmpArray1,FIRST_SAMPLE_INDEX);
        nSamplesLong=(last+1-FIRST_SAMPLE_INDEX)/2;
        if (2*nSamplesLong<minAC) return false;
        avg=0.0; denom=0.0;
        for (i=FIRST_SAMPLE_INDEX+1; i<=last; i+=2) {
            avg+=tmpArray1[i];
            denom+=tmpArray1[i]*tmpArray1[i];
        }
        if (avg<minAC) return false;

        // Allocating space
        if (lastActiveLong==-1) { firstActiveLong=0; lastActiveLong=0; }
        else {
            if ((lastActiveLong+1)%activeLong.length==firstActiveLong) resizeActiveLongs();
            lastActiveLong=(lastActiveLong+1)%activeLong.length;
        }
        if (activeLong[lastActiveLong]==null || activeLong[lastActiveLong].length<last+1) activeLong[lastActiveLong] = new int[last+1];
        System.arraycopy(tmpArray1,0,activeLong[lastActiveLong],0,last+1);
        activeLongLast[lastActiveLong]=last;
        if (activeLongIds[lastActiveLong]==null) activeLongIds[lastActiveLong] = new String[2];
        activeLongIds[lastActiveLong][0]=id; activeLongIds[lastActiveLong][1]=".";

        // Initializing output counts and cached values for computing R^2
        if (rSquaredLong[lastActiveLong]==null) rSquaredLong[lastActiveLong] = new double[7];
        rSquaredLong[lastActiveLong][3]=avg;
        avg/=totalNSamples;
        denom-=totalNSamples*avg*avg;
        rSquaredLong[lastActiveLong][0]=avg;
        rSquaredLong[lastActiveLong][1]=denom;
        rSquaredLong[lastActiveLong][2]=-1.0;
        rSquaredLong[lastActiveLong][4]=nSamplesLong;
        rSquaredLong[lastActiveLong][5]=0;
        rSquaredLong[lastActiveLong][6]=0;

        return true;
    }


    /**
     * @param pos one-based, from the VCF;
     * @param refLength,altLength from the VCF;
     * @param p the starting position of `str` from which to begin parsing
     * (inclusive);
     * @return FALSE iff the number of occurrences of the record is <minAC.
     */
    private static final boolean addActiveShort(String id, int pos, int refLength, int altLength, int type, int length, String str, int p) {
        int i;
        int last, nSamplesShort;
        double avg, denom;

        // Loading the record and filtering by AC
        tmpArray1[0]=pos;
        getInterval(pos,refLength,altLength,type,tmpArray1,1);
        tmpArray1[3]=type;
        tmpArray1[4]=length;
        last=loadSamples(str,p,tmpArray1,FIRST_SAMPLE_INDEX);
        nSamplesShort=(last+1-FIRST_SAMPLE_INDEX)/2;
        if (2*nSamplesShort<minAC) return false;
        avg=0.0; denom=0.0;
        for (i=FIRST_SAMPLE_INDEX+1; i<=last; i+=2) {
            avg+=tmpArray1[i];
            denom+=tmpArray1[i]*tmpArray1[i];
        }
        if (avg<minAC) return false;

        // Allocating space
        if (lastActiveShort==-1) { firstActiveShort=0; lastActiveShort=0; }
        else {
            if ((lastActiveShort+1)%activeShort.length==firstActiveShort) resizeActiveShorts();
            lastActiveShort=(lastActiveShort+1)%activeShort.length;
        }
        if (activeShort[lastActiveShort]==null || activeShort[lastActiveShort].length<last+1) activeShort[lastActiveShort] = new int[last+1];
        System.arraycopy(tmpArray1,0,activeShort[lastActiveShort],0,last+1);
        activeShortLast[lastActiveShort]=last;
        activeShortIds[lastActiveShort]=id;

        // Initializing cached values for computing R^2
        if (rSquaredShort[lastActiveShort]==null) rSquaredShort[lastActiveShort] = new double[4];
        rSquaredShort[lastActiveShort][2]=avg;
        avg/=totalNSamples;
        denom-=totalNSamples*avg*avg;
        rSquaredShort[lastActiveShort][0]=avg;
        rSquaredShort[lastActiveShort][1]=denom;
        rSquaredShort[lastActiveShort][3]=nSamplesShort;

        return true;
    }


    /**
     * @param p the starting position of `str` from which to begin parsing;
     * @param out output array, loaded from cell `outFrom` (inclusive) onwards;
     * @return the last cell of `out` that was loaded.
     */
    private static final int loadSamples(String str, int p, int[] out, int outFrom) {
        final int STR_LENGTH = str.length();
        char c;
        int i, q;
        int last;

        last=outFrom-1; q=p;
        for (i=p; i<STR_LENGTH; i++) {
            c=str.charAt(i);
            if (c=='=' || c=='\t') { 
                out[++last]=Integer.parseInt(str.substring(q,i));
                q=i+1;
            }
        }
        out[++last]=Integer.parseInt(str.substring(q,STR_LENGTH));
        return last;
    }


    /**
     * Removes from the buffer the first few expired shorts that end before 
     * `pos` (zero-based). The procedure stops at the first short that is not 
     * expired, so it is just best-effort heuristic.
     * 
     * Remark: at the end of the procedure, `activeShort` may contain some
     * expired shorts.
     * 
     * @param pos zero-based, inclusive.
     */
    private static final void removeInactiveShorts(int pos) {
        boolean lastReached;

        if (firstActiveShort==-1) return;
        lastReached=false;
        while (true) {
            if (firstActiveShort==lastActiveShort) lastReached=true;
            if (activeShort[firstActiveShort][2]>=pos) break;
            activeShortIds[firstActiveShort]=null;
            if (lastReached) { firstActiveShort=-1; lastActiveShort=-1; break; }
            else firstActiveShort=(firstActiveShort+1)%activeShort.length;
        }
    }


    /**
     * Removes from the buffer the first few expired longs that end before `pos`
     * (zero-based). The procedure stops at the first long that is not expired, 
     * so it is just a best-effort heuristic.
     * 
     * Remark: at the end of the procedure, `activeLong` may contain some
     * expired longs.
     * 
     * @param pos zero-based, inclusive.
     */
    private static final void removeInactiveLongs(int pos, BufferedWriter bw) throws IOException {
        boolean lastReached;

        if (firstActiveLong==-1) return;
        lastReached=false;
        while (true) {
            if (firstActiveLong==lastActiveLong) lastReached=true;
            if (activeLong[firstActiveLong][2]>=pos) break;
            writeToBed(firstActiveLong,bw);
            activeLongIds[firstActiveLong][0]=null; activeLongIds[firstActiveLong][1]=null;
            if (lastReached) { firstActiveLong=-1; lastActiveLong=-1; break; }
            else firstActiveLong=(firstActiveLong+1)%activeLong.length;   
        }
    }


    /**
     * Initializes `rSquaredLong[lastActiveLong]` using all the active shorts
     * that are at distance `<=threshold` from `lastActiveLong` and do not 
     * overlap with it, if any.
     */
    private static final void initializeMaxRsquare(int threshold) {
        boolean lastReached;
        int i;
        int distance, maxIndex;
        double r, max;
        
        if (firstActiveShort==-1) return;
        i=firstActiveShort; lastReached=false; max=rSquaredLong[lastActiveLong][2]; maxIndex=-1;
        while (true) {
            if (i==lastActiveShort) lastReached=true;
            distance=getDistance(lastActiveLong,i);
            if (distance<=threshold && distance>=0) {
                r=rSquare(lastActiveLong,i);
                if (r>=max) { max=r; maxIndex=i; }
            }
            if (lastReached) break;
            i=(i+1)%activeShort.length;
        }
        rSquaredLong[lastActiveLong][2]=max;
        if (maxIndex>=0) {
            rSquaredLong[lastActiveLong][5]=rSquaredShort[maxIndex][2];
            rSquaredLong[lastActiveLong][6]=rSquaredShort[maxIndex][3];
            activeLongIds[lastActiveLong][1]=activeShortIds[maxIndex];
        }
    }


    /**
     * Uses the last active short to update the `rSquaredLong` row of every
     * active long that is at distance `<=threshold` from it and does not
     * overlap with it, if any.
     */
    private static final void updateMaxRsquare(int threshold) {
        boolean lastReached;
        int i;
        int distance;
        double r;
        
        if (firstActiveLong==-1) return;
        i=firstActiveLong; lastReached=false;
        while (true) {
            if (i==lastActiveLong) lastReached=true;
            distance=getDistance(i,lastActiveShort);
            if (distance<=threshold && distance>=0) {
                r=rSquare(i,lastActiveShort);
                if (r>=rSquaredLong[i][2]) {
                    rSquaredLong[i][2]=r;
                    rSquaredLong[i][5]=rSquaredShort[lastActiveShort][2];
                    rSquaredLong[i][6]=rSquaredShort[lastActiveShort][3];
                    activeLongIds[i][1]=activeShortIds[lastActiveShort];
                }
            }
            if (lastReached) break;
            i=(i+1)%activeLong.length;
        }
    }


    /**
     * @return the number of basepairs between two intervals, in any respective
     * order (<0 means that the intervals overlap).
     */
    private static final int getDistance(int longIndex, int shortIndex) {
        // Special case: two INS at the same POS are considered overlapping.
        if (activeLong[longIndex][2]<activeLong[longIndex][1] && activeShort[shortIndex][2]<activeShort[shortIndex][1] && activeLong[longIndex][1]==activeShort[shortIndex][1]) return -1;

        // Every other pair of calls.
        // Remark: an INS overlaps with any interval call that contains or is 
        // identical to its two adjacent positions.
        return Math.max(activeShort[shortIndex][1]-activeLong[longIndex][2]-1,activeLong[longIndex][1]-activeShort[shortIndex][2]-1);
    }


    /**
     * Empties the long buffer to disk and closes `bw`.
     */
    private static final void printRemainingLongs(BufferedWriter bw) throws IOException {
        int i;

        if (firstActiveLong!=-1) { 
            i=firstActiveLong;
            while (true) {
                writeToBed(i,bw);
                if (i==lastActiveLong) break;
                i=(i+1)%activeLong.length;
            }
        }
        bw.close();
    }


    private static final void writeToBed(int longIndex, BufferedWriter bw) throws IOException {
        bw.write( chrId+"\t" +                              // CHROM
                  (activeLong[longIndex][0]-1)+"\t" +       // P
                  activeLong[longIndex][0]+"\t" +           // P+1
                  activeLongIds[longIndex][0]+"\t" +        // ID
                  activeLong[longIndex][3]+"\t" +           // type
                  activeLong[longIndex][4]+"\t" +           // length
                  (int)(rSquaredLong[longIndex][3])+"\t" +  // AC
                  (int)(rSquaredLong[longIndex][4])+"\t" +  // N
                  rSquaredLong[longIndex][2]+"\t" +         // R^2
                  (int)(rSquaredLong[longIndex][5])+"\t" +  // ACs
                  (int)(rSquaredLong[longIndex][6])+"\t" +  // Ns
                  activeLongIds[longIndex][1]+"\n"          // IDs
                );
    }


    /**
     * Loads in `out[outFirst..outFirst+1]` the interval [first..last] (zero-
     * based, inclusive) that corresponds to all and only the positions affected
     * by the variant.
     * 
     * Remark: the interval is empty for INS. This is represented by setting
     * `last<first`.
     */
    private static final void getInterval(int pos, int refLength, int altLength, int type, int[] out, int outFirst) {
        if (type==TYPE_SNP) { out[outFirst]=pos-1; out[outFirst+1]=pos-1; }
        else if (type==TYPE_DEL) { out[outFirst]=(pos-1)+1; out[outFirst+1]=(pos-1)+(refLength-1); }
        else if (type==TYPE_INS) { out[outFirst]=pos; out[outFirst+1]=pos-1; }
        else if (type==TYPE_SUB) { out[outFirst]=pos-1; out[outFirst+1]=(pos-1)+refLength-1; }
        else {
            System.err.println("ERROR: wrong type="+type+" for lengths "+refLength+","+altLength);
            System.exit(1);
        }
    }


    /**
     * @return a non-negative integer.
     */
    private static final int getLength(int refLength, int altLength, int type) {
        if (type==TYPE_SNP) return 1;
        else if (type==TYPE_DEL) return refLength-1;
        else if (type==TYPE_INS) return altLength-1;
        else if (type==TYPE_SUB) return Math.max(refLength,altLength);
        else {
            System.err.println("ERROR: wrong type="+type+" for lengths "+refLength+","+altLength);
            System.exit(1);
            return -1;
        }
    }


    private static final int getType(int refLength, int altLength) {
        if (altLength<refLength) {
            if (altLength==1) return TYPE_DEL;
            else return TYPE_SUB;
        }
        else if (altLength>refLength) {
            if (refLength==1) return TYPE_INS;
            else return TYPE_SUB;
        }
        else {
            if (refLength==1) return TYPE_SNP;
            else return TYPE_SUB;
        }
    }


    private static final void resizeActiveLongs() {
        final int NEW_SIZE = (int)(activeLong.length*RESIZE_FACTOR);

        int[][] newArray = new int[NEW_SIZE][];
        System.arraycopy(activeLong,firstActiveLong,newArray,0,activeLong.length-firstActiveLong);
        if (lastActiveLong<activeLong.length-1) System.arraycopy(activeLong,0,newArray,activeLong.length-firstActiveLong,lastActiveLong+1);
        int[] newArrayLast = new int[NEW_SIZE];
        System.arraycopy(activeLongLast,firstActiveLong,newArrayLast,0,activeLongLast.length-firstActiveLong);
        if (lastActiveLong<activeLongLast.length-1) System.arraycopy(activeLongLast,0,newArrayLast,activeLongLast.length-firstActiveLong,lastActiveLong+1);
        double[][] newArrayRsq = new double[NEW_SIZE][];
        System.arraycopy(rSquaredLong,firstActiveLong,newArrayRsq,0,rSquaredLong.length-firstActiveLong);
        if (lastActiveLong<rSquaredLong.length-1) System.arraycopy(rSquaredLong,0,newArrayRsq,rSquaredLong.length-firstActiveLong,lastActiveLong+1);
        String[][] newArrayIds = new String[NEW_SIZE][];
        System.arraycopy(activeLongIds,firstActiveLong,newArrayIds,0,activeLongIds.length-firstActiveLong);
        if (lastActiveLong<activeLongIds.length-1) System.arraycopy(activeLongIds,0,newArrayIds,activeLongIds.length-firstActiveLong,lastActiveLong+1);
        firstActiveLong=0; lastActiveLong=activeLong.length-1;
        activeLong=newArray; activeLongLast=newArrayLast; rSquaredLong=newArrayRsq; activeLongIds=newArrayIds;
    }


    private static final void resizeActiveShorts() {
        final int NEW_SIZE = (int)(activeShort.length*RESIZE_FACTOR);

        int[][] newArray = new int[NEW_SIZE][];
        System.arraycopy(activeShort,firstActiveShort,newArray,0,activeShort.length-firstActiveShort);
        if (lastActiveShort<activeShort.length-1) System.arraycopy(activeShort,0,newArray,activeShort.length-firstActiveShort,lastActiveShort+1);
        int[] newArrayLast = new int[NEW_SIZE];
        System.arraycopy(activeShortLast,firstActiveShort,newArrayLast,0,activeShortLast.length-firstActiveShort);
        if (lastActiveShort<activeShortLast.length-1) System.arraycopy(activeShortLast,0,newArrayLast,activeShortLast.length-firstActiveShort,lastActiveShort+1);
        double[][] newArrayRsq = new double[NEW_SIZE][];
        System.arraycopy(rSquaredShort,firstActiveShort,newArrayRsq,0,rSquaredShort.length-firstActiveShort);
        if (lastActiveShort<rSquaredShort.length-1) System.arraycopy(rSquaredShort,0,newArrayRsq,rSquaredShort.length-firstActiveShort,lastActiveShort+1);
        String[] newArrayIds = new String[NEW_SIZE];
        System.arraycopy(activeShortIds,firstActiveShort,newArrayIds,0,activeShortIds.length-firstActiveShort);
        if (lastActiveShort<activeShortIds.length-1) System.arraycopy(activeShortIds,0,newArrayIds,activeShortIds.length-firstActiveShort,lastActiveShort+1);
        firstActiveShort=0; lastActiveShort=activeShort.length-1;
        activeShort=newArray; activeShortLast=newArrayLast; rSquaredShort=newArrayRsq; activeShortIds=newArrayIds;
    }


    /**
     * Remark: this procedure essentially just computes the numerator of the R^2
     * coefficient, since the denominator's quantities are already cached in 
     * `rSquaredLong` and `rSquaredShort`.
     * 
     * @return the R^2 coefficient between `activeLong[longIndex]` and 
     * `activeShort[shortIndex]`.
     */
    private static final double rSquare(int longIndex, int shortIndex) {
        final int LAST_LONG = activeLongLast[longIndex];
        final int LAST_SHORT = activeShortLast[shortIndex];
        final int N_SAMPLES_LONG = (int)rSquaredLong[longIndex][4];
        final int LOG2_N_SAMPLES_LONG = 32-Integer.numberOfLeadingZeros(N_SAMPLES_LONG);
        final int N_SAMPLES_SHORT = (int)rSquaredShort[shortIndex][3];
        final int LOG2_N_SAMPLES_SHORT = 32-Integer.numberOfLeadingZeros(N_SAMPLES_SHORT);
        final double AVG_LONG = rSquaredLong[longIndex][0];
        final double AVG_SHORT = rSquaredShort[shortIndex][0];

        int i, j;
        double n, numerator;

        numerator=0.0;
        if (N_SAMPLES_LONG>=N_SAMPLES_SHORT*LOG2_N_SAMPLES_LONG) {
            for (i=FIRST_SAMPLE_INDEX; i<=LAST_SHORT; i+=2) {
                j=binarySearch(activeLong[longIndex],LAST_LONG,activeShort[shortIndex][i]);
                if (j>=0) numerator+=activeLong[longIndex][j+1]*activeShort[shortIndex][i+1];
            }
        }
        else if (N_SAMPLES_SHORT>=N_SAMPLES_LONG*LOG2_N_SAMPLES_SHORT) {
            for (i=FIRST_SAMPLE_INDEX; i<=LAST_LONG; i+=2) {
                j=binarySearch(activeShort[shortIndex],LAST_SHORT,activeLong[longIndex][i]);
                if (j>=0) numerator+=activeLong[longIndex][i+1]*activeShort[shortIndex][j+1];
            }
        }
        else {
            i=FIRST_SAMPLE_INDEX; j=FIRST_SAMPLE_INDEX;
            while (i<=LAST_LONG && j<=LAST_SHORT) {
                if (activeLong[longIndex][i]<activeShort[shortIndex][j]) i+=2;
                else if (activeLong[longIndex][i]>activeShort[shortIndex][j]) j+=2;
                else {
                    numerator+=activeLong[longIndex][i+1]*activeShort[shortIndex][j+1];
                    i+=2; j+=2;
                }
            }
        }
        numerator-=totalNSamples*AVG_LONG*AVG_SHORT;
        n=numerator/Math.sqrt(rSquaredLong[longIndex][1]*rSquaredShort[shortIndex][1]);
        return n*n;
    }


    /**
     * @param array one of the `active*` rows;
     * @param last one of the `active*Last` values.
     */
    private static final int binarySearch(int[] array, int last, int key) {
        int low, mid, high, value;
        
        low=0;
        high=(last+1-FIRST_SAMPLE_INDEX)/2-1;
        while (low <= high) {
            mid=(low+high)>>>1;
            value=array[FIRST_SAMPLE_INDEX+mid*2];
            if (value<key) low=mid+1;
            else if (value>key) high=mid-1;
            else return FIRST_SAMPLE_INDEX+mid*2;
        }
        return -1;
    }
    
}