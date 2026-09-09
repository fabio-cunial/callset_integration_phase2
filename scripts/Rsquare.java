import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * For each SV, computes the max R^2 coefficient with the SNPs located within a 
 * max distance from the SV's endpoints and not overlapping the SV.
 * 
 * Remark: the distinction between SVs and SNPs is only based on length. In 
 * particular, the SNP category includes INS, DEL, replacements. Every SNP that
 * affects bases within the max distance and does not overlap the SV contributes
 * to the SV's R^2.
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
     * where `[first..last]` is an interval with zero-based, inclusive 
     * positions, and `countX` is in {0,1,2}.
     * 
     * Remark: rows may have different lengths.
     */
    private static int[][] activeSvs, activeSnps;
    private static int[] activeSvsLast, activeSnpsLast;
    private static int firstActiveSv, firstActiveSnp, lastActiveSv, lastActiveSnp;
    private static int totalNSamples;

    /**
     * Format of each row: 
     * 
     * tmpAvg, tmpD, R^2, AC, N
     * 
     * where `tmp*` are just temporary variables that are cached for speeding up
     * computation, and `N` is the number of samples with GT=ALT.
     */
    private static double[][] rSquaredSv;

    /**
     * Format of each row: `tmpAvg, tmpD`.
     */
    private static double[][] rSquaredSnp;

    /**
     * Scratch space
     */
    private static int[] tmpArray1;


    /**
     * @param args 0 format: `POS,refLength,altLength, sample=gtCount,...`,
     * where `POS` is sorted and comes from the VCF, and `gtCount` is an integer
     * in {0,1,2}; all the records are expected to come from the same chrom;
     * @param args 4 format: `svtype,svlen,R^2,AC,N`, where `N` is the number of
     * samples with GT=ALT and R^2 is -1 iff it cannot be compared to any SNP.
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
        long nRecords, nSvs, nSnps;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        
        // Allocating reused space
        activeSvs = new int[CAPACITY][4+2*CAPACITY];
        activeSvsLast = new int[CAPACITY];
        firstActiveSv=-1; lastActiveSv=-1;
        activeSnps = new int[CAPACITY][4+2*CAPACITY];
        activeSnpsLast = new int[CAPACITY];
        firstActiveSnp=-1; lastActiveSnp=-1;
        tmpArray1 = new int[4+2*totalNSamples];
        rSquaredSv = new double[CAPACITY][5];
        rSquaredSnp = new double[CAPACITY][2];

        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_TSV))));
        bw = new BufferedWriter(new FileWriter(OUTPUT_TSV));
        str=br.readLine(); nRecords=0; nSvs=0; nSnps=0;
        while (str!=null) {
            p1=str.indexOf('\t');
            pos=Integer.parseInt(str.substring(0,p1));
            p2=str.indexOf('\t',p1+1);
            refLength=Integer.parseInt(str.substring(p1+1,p2));
            p3=str.indexOf('\t',p2+1);
            altLength=Integer.parseInt(str.substring(p2+1,p3));
            length=altLength-refLength;
            if (length<=-MIN_SV_LENGTH || length>=MIN_SV_LENGTH) {
                nSvs++;
                addActiveSv(pos,refLength,altLength,str,p3+1);
                removeInactiveSnps(pos-1-MAX_DISTANCE_BP-1);
                removeInactiveSvs(pos-1-MAX_DISTANCE_BP-1,bw);
                initializeMaxRsquare(MAX_DISTANCE_BP);
            }
            else {
                nSnps++;
                addActiveSnp(pos,refLength,altLength,str,p3+1);
                removeInactiveSvs(pos-1-MAX_DISTANCE_BP-1,bw);
                removeInactiveSnps(pos-1-MAX_DISTANCE_BP-1);
                updateMaxRsquare(MAX_DISTANCE_BP);
            }

            // Next iteration
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            str=br.readLine();
        }
        br.close();
        printRemainingSvs(bw);
        System.err.println("Processed "+nRecords+" total records, "+nSvs+" SVs, "+nSnps+" SNPs.");
    }


    /**
     * @param pos one-based, from the VCF;
     * @param refLength,altLength from the VCF;
     * @param p the starting position of `str` from which to begin parsing
     * (inclusive).
     */
    private static final void addActiveSv(int pos, int refLength, int altLength, String str, int p) {
        int i;
        int last, nSamplesSv;
        double avg, denom;

        // Allocating space
        if (lastActiveSv==-1) { firstActiveSv=0; lastActiveSv=0; }
        else {
            if ((lastActiveSv+1)%activeSvs.length==firstActiveSv) resizeActiveSvs();
            lastActiveSv=(lastActiveSv+1)%activeSvs.length;
        }

        // Loading the record
        getInterval(pos,refLength,altLength,tmpArray1);
        tmpArray1[2]=getType(refLength,altLength);
        tmpArray1[3]=Math.abs(altLength-refLength);
        last=loadSamples(str,p,tmpArray1,4);
        if (activeSvs[lastActiveSv]==null || activeSvs[lastActiveSv].length<last+1) activeSvs[lastActiveSv] = new int[last+1];
        System.arraycopy(tmpArray1,0,activeSvs[lastActiveSv],0,last+1);
        activeSvsLast[lastActiveSv]=last;
        nSamplesSv=(last+1-FIRST_SAMPLE_INDEX)/2;

        // Initializing output counts and cached values for computing R^2
        if (rSquaredSv[lastActiveSv]==null) rSquaredSv[lastActiveSv] = new double[5];
        avg=0.0; denom=0.0;
        for (i=FIRST_SAMPLE_INDEX+1; i<=last; i+=2) {
            avg+=activeSvs[lastActiveSv][i];
            denom+=activeSvs[lastActiveSv][i]*activeSvs[lastActiveSv][i];
        }
        rSquaredSv[lastActiveSv][3]=avg;
        avg/=totalNSamples;
        denom-=totalNSamples*avg*avg;
        rSquaredSv[lastActiveSv][0]=avg;
        rSquaredSv[lastActiveSv][1]=denom;
        rSquaredSv[lastActiveSv][2]=-1.0;
        rSquaredSv[lastActiveSv][4]=nSamplesSv;
    }


    /**
     * @param pos one-based, from the VCF;
     * @param refLength,altLength from the VCF;
     * @param p the starting position of `str` from which to begin parsing
     * (inclusive).
     */
    private static final void addActiveSnp(int pos, int refLength, int altLength, String str, int p) {
        int i;
        int last;
        double avg, denom;

        // Allocating space
        if (lastActiveSnp==-1) { firstActiveSnp=0; lastActiveSnp=0; }
        else {
            if ((lastActiveSnp+1)%activeSnps.length==firstActiveSnp) resizeActiveSnps();
            lastActiveSnp=(lastActiveSnp+1)%activeSnps.length;
        }

        // Loading the record
        getInterval(pos,refLength,altLength,tmpArray1);
        tmpArray1[2]=getType(refLength,altLength);
        tmpArray1[3]=Math.abs(altLength-refLength);
        last=loadSamples(str,p,tmpArray1,4);
        if (activeSnps[lastActiveSnp]==null || activeSnps[lastActiveSnp].length<last+1) activeSnps[lastActiveSnp] = new int[last+1];
        System.arraycopy(tmpArray1,0,activeSnps[lastActiveSnp],0,last+1);
        activeSnpsLast[lastActiveSnp]=last;

        // Initializing cached values for computing R^2
        if (rSquaredSnp[lastActiveSnp]==null) rSquaredSnp[lastActiveSnp] = new double[2];
        avg=0.0; denom=0.0;
        for (i=FIRST_SAMPLE_INDEX+1; i<=last; i+=2) {
            avg+=activeSnps[lastActiveSnp][i];
            denom+=activeSnps[lastActiveSnp][i]*activeSnps[lastActiveSnp][i];
        }
        avg/=totalNSamples;
        denom-=totalNSamples*avg*avg;
        rSquaredSnp[lastActiveSnp][0]=avg;
        rSquaredSnp[lastActiveSnp][1]=denom;
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
     * Removes from the buffer the first few expired SNPs that end before `pos`
     * (zero-based). The procedure stops at the first SNP that is not expired, 
     * so it is just best-effort heuristic.
     * 
     * Remark: at the end of the procedure, `activeSnps` may contain some
     * expired SNPs.
     * 
     * @param pos zero-based, inclusive.
     */
    private static final void removeInactiveSnps(int pos) {
        boolean lastReached;

        if (firstActiveSnp==-1) return;
        lastReached=false;
        while (true) {
            if (firstActiveSnp==lastActiveSnp) lastReached=true;
            if (activeSnps[firstActiveSnp][1]>=pos) break;
            if (lastReached) { firstActiveSnp=-1; lastActiveSnp=-1; break; }
            else firstActiveSnp=(firstActiveSnp+1)%activeSnps.length;
        }
    }


    /**
     * Removes from the buffer the first few expired SVs that end before `pos`
     * (zero-based). The procedure stops at the first SV that is not expired, 
     * so it is just a best-effort heuristic.
     * 
     * Remark: at the end of the procedure, `activeSvs` may contain some
     * expired SVs.
     * 
     * @param pos zero-based, inclusive.
     */
    private static final void removeInactiveSvs(int pos, BufferedWriter bw) throws IOException {
        boolean lastReached;

        if (firstActiveSv==-1) return;
        lastReached=false;
        while (true) {
            if (firstActiveSv==lastActiveSv) lastReached=true;
            if (activeSvs[firstActiveSv][1]>=pos) break;
            bw.write(activeSvs[firstActiveSv][2]+"\t"+activeSvs[firstActiveSv][3]+"\t"+rSquaredSv[firstActiveSv][2]+"\t"+rSquaredSv[firstActiveSv][3]+"\t"+rSquaredSv[firstActiveSv][4]+"\n");
            if (lastReached) { firstActiveSv=-1; lastActiveSv=-1; break; }
            else firstActiveSv=(firstActiveSv+1)%activeSvs.length;   
        }
    }


    /**
     * Initializes `rSquaredSv[lastActiveSv][2]` using all the active SNPs that 
     * are at distance `<=threshold` from `lastActiveSv` and do not overlap with
     * it, if any.
     */
    private static final void initializeMaxRsquare(int threshold) {
        boolean lastReached;
        int i, distance;
        double r, max;
        
        if (firstActiveSnp==-1) return;
        i=firstActiveSnp; lastReached=false; max=rSquaredSv[lastActiveSv][2];
        while (true) {
            if (i==lastActiveSnp) lastReached=true;
            distance=getDistance(lastActiveSv,i);
            if (distance<=threshold && distance>=0) {
                r=rSquare(lastActiveSv,i);
                if (r>=max) max=r;
            }
            if (lastReached) break;
            i=(i+1)%activeSnps.length;
        }
        rSquaredSv[lastActiveSv][2]=max;
    }


    /**
     * Uses the last active SNP to update the `rSquared` value of every active 
     * SV that is at distance `<=threshold` from it and does not overlap with 
     * it, if any.
     */
    private static final void updateMaxRsquare(int threshold) {
        boolean lastReached;
        int i, distance;
        double r;
        
        if (firstActiveSv==-1) return;
        i=firstActiveSv; lastReached=false;
        while (true) {
            if (i==lastActiveSv) lastReached=true;
            distance=getDistance(i,lastActiveSnp);
            if (distance<=threshold && distance>=0) {
                r=rSquare(i,lastActiveSnp);
                if (r>=rSquaredSv[i][2]) rSquaredSv[i][2]=r;
            }
            if (lastReached) break;
            i=(i+1)%activeSvs.length;
        }
    }


    /**
     * @return the number of basepairs between two intervals, in any respective
     * order (<0 means that the intervals overlap).
     */
    private static final int getDistance(int svIndex, int snpIndex) {
        return Math.max(activeSnps[snpIndex][0]-activeSvs[svIndex][1]-1,activeSvs[svIndex][0]-activeSnps[snpIndex][1]-1);
    }


    /**
     * Empties the SV buffer to disk and closes `bw`.
     */
    private static final void printRemainingSvs(BufferedWriter bw) throws IOException {
        int i;

        if (firstActiveSv!=-1) { 
            i=firstActiveSv;
            while (true) {
                bw.write(activeSvs[i][2]+"\t"+activeSvs[i][3]+"\t"+rSquaredSv[i][2]+"\t"+rSquaredSv[i][3]+"\t"+rSquaredSv[i][4]+"\n");
                if (i==lastActiveSv) break;
                i=(i+1)%activeSvs.length;
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


    private static final void resizeActiveSvs() {
        final int NEW_SIZE = (int)(activeSvs.length*RESIZE_FACTOR);

        int[][] newArray = new int[NEW_SIZE][];
        System.arraycopy(activeSvs,firstActiveSv,newArray,0,activeSvs.length-firstActiveSv);
        if (lastActiveSv<activeSvs.length-1) System.arraycopy(activeSvs,0,newArray,activeSvs.length-firstActiveSv,lastActiveSv+1);
        int[] newArrayLast = new int[NEW_SIZE];
        System.arraycopy(activeSvsLast,firstActiveSv,newArrayLast,0,activeSvsLast.length-firstActiveSv);
        if (lastActiveSv<activeSvsLast.length-1) System.arraycopy(activeSvsLast,0,newArrayLast,activeSvsLast.length-firstActiveSv,lastActiveSv+1);
        double[][] newArrayRsq = new double[NEW_SIZE][];
        System.arraycopy(rSquaredSv,firstActiveSv,newArrayRsq,0,rSquaredSv.length-firstActiveSv);
        if (lastActiveSv<rSquaredSv.length-1) System.arraycopy(rSquaredSv,0,newArrayRsq,rSquaredSv.length-firstActiveSv,lastActiveSv+1);
        firstActiveSv=0; lastActiveSv=activeSvs.length-1;
        activeSvs=newArray; activeSvsLast=newArrayLast; rSquaredSv=newArrayRsq;
    }


    private static final void resizeActiveSnps() {
        final int NEW_SIZE = (int)(activeSnps.length*RESIZE_FACTOR);

        int[][] newArray = new int[NEW_SIZE][];
        System.arraycopy(activeSnps,firstActiveSnp,newArray,0,activeSnps.length-firstActiveSnp);
        if (lastActiveSnp<activeSnps.length-1) System.arraycopy(activeSnps,0,newArray,activeSnps.length-firstActiveSnp,lastActiveSnp+1);
        int[] newArrayLast = new int[NEW_SIZE];
        System.arraycopy(activeSnpsLast,firstActiveSnp,newArrayLast,0,activeSnpsLast.length-firstActiveSnp);
        if (lastActiveSnp<activeSnpsLast.length-1) System.arraycopy(activeSnpsLast,0,newArrayLast,activeSnpsLast.length-firstActiveSnp,lastActiveSnp+1);
        double[][] newArrayRsq = new double[NEW_SIZE][];
        System.arraycopy(rSquaredSnp,firstActiveSnp,newArrayRsq,0,rSquaredSnp.length-firstActiveSnp);
        if (lastActiveSnp<rSquaredSnp.length-1) System.arraycopy(rSquaredSnp,0,newArrayRsq,rSquaredSnp.length-firstActiveSnp,lastActiveSnp+1);
        firstActiveSnp=0; lastActiveSnp=activeSnps.length-1;
        activeSnps=newArray; activeSnpsLast=newArrayLast; rSquaredSnp=newArrayRsq;
    }


    /**
     * Remark: this procedure essentially just computes the numerator of the R^2
     * coefficient, since the denominator's quantities are already cached in 
     * `rSquaredSv` and `rSquaredSnp`.
     * 
     * @return the R^2 coefficient between `activeSvs[svIndex]` and 
     * `activeSnps[snpIndex]`.
     */
    private static final double rSquare(int svIndex, int snpIndex) {
        final int LAST_SV = activeSvsLast[svIndex];
        final int LAST_SNP = activeSnpsLast[snpIndex];
        final double AVG_SV = rSquaredSv[svIndex][0];
        final double AVG_SNP = rSquaredSnp[snpIndex][0];
        int i, j;
        double n, numerator;

        numerator=0.0;
        i=FIRST_SAMPLE_INDEX; j=FIRST_SAMPLE_INDEX;
        while (i<=LAST_SV && j<=LAST_SNP) {
            if (activeSvs[svIndex][i]<activeSnps[snpIndex][j]) i+=2;
            else if (activeSvs[svIndex][i]>activeSnps[snpIndex][j]) j+=2;
            else {
                numerator+=activeSvs[svIndex][i+1]*activeSnps[snpIndex][j+1];
                i+=2; j+=2;
            }
        }
        numerator-=totalNSamples*AVG_SV*AVG_SNP;
        n=numerator/Math.sqrt(rSquaredSv[svIndex][1]*rSquaredSnp[snpIndex][1]);
        return n*n;
    }
    
}