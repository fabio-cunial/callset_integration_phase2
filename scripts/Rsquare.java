import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * For each long call, computes the max R^2 coefficient with the short calls 
 * located within a max distance from the long's endpoints and not overlapping 
 * the long.
 * 
 * Remark: both the short and the long categories include INS, DEL and 
 * replacements. Every short that affects bases within the max distance and does
 * not overlap the long contributes to the long's R^2.
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
    private static final int FIRST_SAMPLE_INDEX = 4;

    /**
     * Format of each row:
     * 
     * first,last,type,length,  sample1,count1, sample2,counts2, ..., sampleK,countK
     * 
     * where `[first..last]` is the interval of zero-based, inclusive positions
     * that are affected by the variant, and `countX` is in {0,1,2}.
     * 
     * Remark: rows may have different lengths.
     */
    private static int[][] activeLong, activeShort;
    private static int[] activeLongLast, activeShortLast;
    private static int firstActiveLong, firstActiveShort, lastActiveLong, lastActiveShort;
    private static int totalNSamples;

    /**
     * Format of each row: 
     * 
     * tmpAvg, tmpD, R^2, AC, N
     * 
     * where `tmp*` are just temporary variables that are cached for speeding up
     * computation, and `N` is the number of samples with GT=ALT.
     */
    private static double[][] rSquaredLong;

    /**
     * Format of each row: `tmpAvg, tmpD`.
     */
    private static double[][] rSquaredShort;

    /**
     * Scratch space
     */
    private static int[] tmpArray1;


    /**
     * @param args 0 format: `POS,refLength,altLength, sample=gtCount,...`,
     * where `POS` is sorted and comes from the VCF, and `gtCount` is an integer
     * in {0,1,2}; all the records are expected to come from the same chrom;
     * @param args 4 format: `type,length,R^2,AC,N`, where `N` is the number of
     * samples with GT=ALT and R^2 is -1 iff it cannot be compared to any short.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV = args[0];
        totalNSamples=Integer.parseInt(args[1]);
        final int MIN_SV_LENGTH = Integer.parseInt(args[2]);
        final int MAX_DISTANCE_BP = Integer.parseInt(args[3]);
        final String OUTPUT_TSV = args[4];

        final int QUANTUM = 1000;  // Arbitrary

        int p1, p2, p3;
        int pos, length, refLength, altLength;
        long nRecords, nLong, nShort;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        
        // Allocating reused space
        activeLong = new int[CAPACITY][4+2*CAPACITY];
        activeLongLast = new int[CAPACITY];
        firstActiveLong=-1; lastActiveLong=-1;
        activeShort = new int[CAPACITY][4+2*CAPACITY];
        activeShortLast = new int[CAPACITY];
        firstActiveShort=-1; lastActiveShort=-1;
        tmpArray1 = new int[4+2*totalNSamples];
        rSquaredLong = new double[CAPACITY][5];
        rSquaredShort = new double[CAPACITY][2];

        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_TSV))));
        bw = new BufferedWriter(new FileWriter(OUTPUT_TSV));
        str=br.readLine(); nRecords=0; nLong=0; nShort=0;
        while (str!=null) {
            p1=str.indexOf('\t');
            pos=Integer.parseInt(str.substring(0,p1));
            p2=str.indexOf('\t',p1+1);
            refLength=Integer.parseInt(str.substring(p1+1,p2));
            p3=str.indexOf('\t',p2+1);
            altLength=Integer.parseInt(str.substring(p2+1,p3));
            length=altLength-refLength;
            if (length<=-MIN_SV_LENGTH || length>=MIN_SV_LENGTH) {
                nLong++;
                addActiveLong(pos,refLength,altLength,str,p3+1);
                removeInactiveShorts(pos-1-MAX_DISTANCE_BP-1);
                removeInactiveLongs(pos-1-MAX_DISTANCE_BP-1,bw);
                initializeMaxRsquare(MAX_DISTANCE_BP);
            }
            else {
                nShort++;
                addActiveShort(pos,refLength,altLength,str,p3+1);
                removeInactiveLongs(pos-1-MAX_DISTANCE_BP-1,bw);
                removeInactiveShorts(pos-1-MAX_DISTANCE_BP-1);
                updateMaxRsquare(MAX_DISTANCE_BP);
            }

            // Next iteration
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            str=br.readLine();
        }
        br.close();
        printRemainingLong(bw);
        System.err.println("Processed "+nRecords+" total records, "+nLong+" long, "+nShort+" short.");
    }


    /**
     * @param pos one-based, from the VCF;
     * @param refLength,altLength from the VCF;
     * @param p the starting position of `str` from which to begin parsing
     * (inclusive).
     */
    private static final void addActiveLong(int pos, int refLength, int altLength, String str, int p) {
        int i;
        int last, nSamplesLong;
        double avg, denom;

        // Allocating space
        if (lastActiveLong==-1) { firstActiveLong=0; lastActiveLong=0; }
        else {
            if ((lastActiveLong+1)%activeLong.length==firstActiveLong) resizeActiveLongs();
            lastActiveLong=(lastActiveLong+1)%activeLong.length;
        }

        // Loading the record
        getInterval(pos,refLength,altLength,tmpArray1);
        tmpArray1[2]=getType(refLength,altLength);
        tmpArray1[3]=Math.abs(altLength-refLength);
        last=loadSamples(str,p,tmpArray1,4);
        if (activeLong[lastActiveLong]==null || activeLong[lastActiveLong].length<last+1) activeLong[lastActiveLong] = new int[last+1];
        System.arraycopy(tmpArray1,0,activeLong[lastActiveLong],0,last+1);
        activeLongLast[lastActiveLong]=last;
        nSamplesLong=(last+1-FIRST_SAMPLE_INDEX)/2;

        // Initializing output counts and cached values for computing R^2
        if (rSquaredLong[lastActiveLong]==null) rSquaredLong[lastActiveLong] = new double[5];
        avg=0.0; denom=0.0;
        for (i=FIRST_SAMPLE_INDEX+1; i<=last; i+=2) {
            avg+=activeLong[lastActiveLong][i];
            denom+=activeLong[lastActiveLong][i]*activeLong[lastActiveLong][i];
        }
        rSquaredLong[lastActiveLong][3]=avg;
        avg/=totalNSamples;
        denom-=totalNSamples*avg*avg;
        rSquaredLong[lastActiveLong][0]=avg;
        rSquaredLong[lastActiveLong][1]=denom;
        rSquaredLong[lastActiveLong][2]=-1.0;
        rSquaredLong[lastActiveLong][4]=nSamplesLong;
    }


    /**
     * @param pos one-based, from the VCF;
     * @param refLength,altLength from the VCF;
     * @param p the starting position of `str` from which to begin parsing
     * (inclusive).
     */
    private static final void addActiveShort(int pos, int refLength, int altLength, String str, int p) {
        int i;
        int last;
        double avg, denom;

        // Allocating space
        if (lastActiveShort==-1) { firstActiveShort=0; lastActiveShort=0; }
        else {
            if ((lastActiveShort+1)%activeShort.length==firstActiveShort) resizeActiveShorts();
            lastActiveShort=(lastActiveShort+1)%activeShort.length;
        }

        // Loading the record
        getInterval(pos,refLength,altLength,tmpArray1);
        tmpArray1[2]=getType(refLength,altLength);
        tmpArray1[3]=Math.abs(altLength-refLength);
        last=loadSamples(str,p,tmpArray1,4);
        if (activeShort[lastActiveShort]==null || activeShort[lastActiveShort].length<last+1) activeShort[lastActiveShort] = new int[last+1];
        System.arraycopy(tmpArray1,0,activeShort[lastActiveShort],0,last+1);
        activeShortLast[lastActiveShort]=last;

        // Initializing cached values for computing R^2
        if (rSquaredShort[lastActiveShort]==null) rSquaredShort[lastActiveShort] = new double[2];
        avg=0.0; denom=0.0;
        for (i=FIRST_SAMPLE_INDEX+1; i<=last; i+=2) {
            avg+=activeShort[lastActiveShort][i];
            denom+=activeShort[lastActiveShort][i]*activeShort[lastActiveShort][i];
        }
        avg/=totalNSamples;
        denom-=totalNSamples*avg*avg;
        rSquaredShort[lastActiveShort][0]=avg;
        rSquaredShort[lastActiveShort][1]=denom;
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
            if (activeShort[firstActiveShort][1]>=pos) break;
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
            if (activeLong[firstActiveLong][1]>=pos) break;
            bw.write(activeLong[firstActiveLong][2]+"\t"+activeLong[firstActiveLong][3]+"\t"+rSquaredLong[firstActiveLong][2]+"\t"+rSquaredLong[firstActiveLong][3]+"\t"+rSquaredLong[firstActiveLong][4]+"\n");
            if (lastReached) { firstActiveLong=-1; lastActiveLong=-1; break; }
            else firstActiveLong=(firstActiveLong+1)%activeLong.length;   
        }
    }


    /**
     * Initializes `rSquaredLong[lastActiveLong][2]` using all the active shorts
     * that are at distance `<=threshold` from `lastActiveLong` and do not 
     * overlap with it, if any.
     */
    private static final void initializeMaxRsquare(int threshold) {
        boolean lastReached;
        int i, distance;
        double r, max;
        
        if (firstActiveShort==-1) return;
        i=firstActiveShort; lastReached=false; max=rSquaredLong[lastActiveLong][2];
        while (true) {
            if (i==lastActiveShort) lastReached=true;
            distance=getDistance(lastActiveLong,i);
            if (distance<=threshold && distance>=0) {
                r=rSquare(lastActiveLong,i);
                if (r>=max) max=r;
            }
            if (lastReached) break;
            i=(i+1)%activeShort.length;
        }
        rSquaredLong[lastActiveLong][2]=max;
    }


    /**
     * Uses the last active short to update the `rSquared` value of every active 
     * long that is at distance `<=threshold` from it and does not overlap with 
     * it, if any.
     */
    private static final void updateMaxRsquare(int threshold) {
        boolean lastReached;
        int i, distance;
        double r;
        
        if (firstActiveLong==-1) return;
        i=firstActiveLong; lastReached=false;
        while (true) {
            if (i==lastActiveLong) lastReached=true;
            distance=getDistance(i,lastActiveShort);
            if (distance<=threshold && distance>=0) {
                r=rSquare(i,lastActiveShort);
                if (r>=rSquaredLong[i][2]) rSquaredLong[i][2]=r;
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
        return Math.max(activeShort[shortIndex][0]-activeLong[longIndex][1]-1,activeLong[longIndex][0]-activeShort[shortIndex][1]-1);
    }


    /**
     * Empties the long buffer to disk and closes `bw`.
     */
    private static final void printRemainingLong(BufferedWriter bw) throws IOException {
        int i;

        if (firstActiveLong!=-1) { 
            i=firstActiveLong;
            while (true) {
                bw.write(activeLong[i][2]+"\t"+activeLong[i][3]+"\t"+rSquaredLong[i][2]+"\t"+rSquaredLong[i][3]+"\t"+rSquaredLong[i][4]+"\n");
                if (i==lastActiveLong) break;
                i=(i+1)%activeLong.length;
            }
        }
        bw.close();
    }


    /**
     * Loads in `out` the interval [first..last] (zero-based, inclusive) that
     * corresponds to all and only the positions affected by the variant.
     * 
     * Remark: the interval is empty for INS. This is represented by setting
     * `last<first`.
     */
    private static final void getInterval(int pos, int refLength, int altLength, int[] out) {
        final int TYPE = getType(refLength,altLength);

        if (TYPE==TYPE_SNP) { out[0]=pos-1; out[1]=pos-1; }
        else if (TYPE==TYPE_DEL) { out[0]=(pos-1)+1; out[1]=(pos-1)+(refLength-1); }
        else if (TYPE==TYPE_INS) { out[0]=pos; out[1]=pos-1; }
        else if (TYPE==TYPE_SUB) { out[0]=pos-1; out[1]=(pos-1)+refLength-1; }
        else {
            System.err.println("ERROR: wrong type="+TYPE+" for lengths "+refLength+","+altLength);
            System.exit(1);
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
        firstActiveLong=0; lastActiveLong=activeLong.length-1;
        activeLong=newArray; activeLongLast=newArrayLast; rSquaredLong=newArrayRsq;
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
        firstActiveShort=0; lastActiveShort=activeShort.length-1;
        activeShort=newArray; activeShortLast=newArrayLast; rSquaredShort=newArrayRsq;
    }


    /**
     * Remark: this procedure essentially just computes the numerator of the R^2
     * coefficient, since the denominator's quantities are already cached in 
     * `rSquaredLong` and `rSquaredShort`.
     * 
     * @return the R^2 coefficient between `activeLong[svIndex]` and 
     * `activeShort[snpIndex]`.
     */
    private static final double rSquare(int svIndex, int snpIndex) {
        final int LAST_SV = activeLongLast[svIndex];
        final int LAST_SNP = activeShortLast[snpIndex];
        final double AVG_SV = rSquaredLong[svIndex][0];
        final double AVG_SNP = rSquaredShort[snpIndex][0];
        int i, j;
        double n, numerator;

        numerator=0.0;
        i=FIRST_SAMPLE_INDEX; j=FIRST_SAMPLE_INDEX;
        while (i<=LAST_SV && j<=LAST_SNP) {
            if (activeLong[svIndex][i]<activeShort[snpIndex][j]) i+=2;
            else if (activeLong[svIndex][i]>activeShort[snpIndex][j]) j+=2;
            else {
                numerator+=activeLong[svIndex][i+1]*activeShort[snpIndex][j+1];
                i+=2; j+=2;
            }
        }
        numerator-=totalNSamples*AVG_SV*AVG_SNP;
        n=numerator/Math.sqrt(rSquaredLong[svIndex][1]*rSquaredShort[snpIndex][1]);
        return n*n;
    }
    
}