import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * 
 */
public class Rsquare {

    /**
     * Reused buffers
     */
    private static int firstActiveSv, firstActiveSnps, lastActiveSv, lastActiveSnp;
    private static double[] rSquared;

    /**
     * Row: pos,len,sample1,count1,sample2,counts2,...,sampleK,countK
     */
    private static int[][] activeSvs, activeSnps;
    private static int[] activeSvsLast, activeSnpsLast;

    /**
     * Scratch space
     */
    private static boolean[] tmpArray1, tmpArray2;


    /**
     * @param args 0 format: POS,REF,ALT,sample=GT,sample=GT,...,sample=GT
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV = args[0];
        final int N_SAMPLES = Integer.parseInt(args[1]);
        final int MIN_SV_LENGTH = Integer.parseInt(args[2]);
        final int MAX_DISTANCE_BP = Integer.parseInt(args[3]);
        
        
        int i, p, q, r, s, t;
        int chunkID, pos, refLength, altLength, nRecords;
        int first, last, maxLast;  // One-based, inclusive.
        String str, headerStr;
        StringBuilder header;
        BufferedReader br;
        BufferedWriter bw;
        


        
        
        boolean isDel, isIns, isSv;
        int length, firstActiveSv, firstActiveSnps, lastActiveSv, lastActiveSnp;
        long nRecords;
        int[] tmpArray;
        
        
        
        
        // Allocating reused memory
        activeSvs = new int[CAPACITY][2+2*N_SAMPLES];
        activeSvsLast = new int[CAPACITY];
        firstActiveSv=-1; lastActiveSv=-1;
        activeSnps = new int[CAPACITY][2+2*N_SAMPLES];
        activeSnpsLast = new int[CAPACITY];
        firstActiveSnps=-1; lastActiveSnp=-1;
        tmpArray1 = new boolean[2+2*N_SAMPLES];
        tmpArray2 = new boolean[2+2*N_SAMPLES];
        rSquared = new double[CAPACITY];


        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_TSV))));
        str=br.readLine(); nRecords=0;
        while (str!=null) {
            nRecords++;
            p1=str.indexOf('\t');
            pos=Integer.parseInt(str.substring(0,p1));
            p2=str.indexOf('\t',p1+1);
            p3=str.indexOf('\t',p2+1);
            length=(p3-p2-1)-(p2-p1-1);
            isDel=length<0; isIns=length>0; isSv=(length<=-MIN_SV_LENGTH)||(length>=MIN_SV_LENGTH);
            if (isSv) {
                // Adding a new active SV
                tmpArray[0]=pos; tmpArray[1]=length;
                last=loadRecord(str,p3+1,tmpArray);
                allocateActiveSv();
                System.arraycopy(tmpArray,0,activeSvs[lastActiveSv],0,last+1);
                activeSvsLast[lastActiveSv]=last;
                // Removing inactive SNPs
                removeInactiveSnps(pos-MAX_DISTANCE_BP);
                // Computing max R^2
                computeMaxRsquare(lastActiveSv);
            }
            else {
                // Adding a new active SNP

                // Removing and printing inactive SVs

            }
            str=br.readLine();
        }
        bw.close();
        System.err.println("chunk="+chunkID+" nRecords="+nRecords);
        br.close();
        System.err.println("Created "+(chunkID+1)+" truvari chunks");
    }


    /**
     * @param p the starting position of `str` from which to begin parsing;
     * @return the last cell of the output array that was loaded.
     */
    private static final int loadRecord(String str, int p, int[] tmpArray) {
        final int STR_LENGTH = str.length();
        char c;
        int i, q;
        int last;

        last=1; q=p;
        for (i=p; i<STR_LENGTH; i++) {
            c=str.charAt(i);
            if (c=='=' || c=='\t') { 
                tmpArray[++last]=Integer.parseInt(str.substring(q,i));
                q=i+1;
            }
        }
        tmpArray[++last]=Integer.parseInt(str.substring(q,STR_LENGTH));
        return last;
    }


    /**
     * Makes room for a new SV in the reused buffers.
     */
    private static final void allocateActiveSv() {
        if (lastActiveSv==-1) { firstActiveSv=0; lastActiveSv=0; }
        else {
            if ((lastActiveSv+1)%activeSvs.length==firstActiveSv) resizeActiveSvs();
            lastActiveSv=(lastActiveSv+1)%activeSvs.length;
        }
    }


    /**
     * Makes room for a new SNP in the reused buffers.
     */
    private static final void allocateActiveSnp() {
        
    }



    private static final void resizeActiveSvs() {


    }


    private static final void resizeActiveSnps() {

        
    }


    /**
     * Removes from the buffer every SNP that occurs before `oldestPos`.
     */
    private static final void removeInactiveSnps(int oldestPos) {
        boolean lastReached;
        int i;

        i=firstActiveSnp; lastReached=false;
        while (true) {
            if (i==lastActiveSnp) lastReached=true;
            if (activeSnps[i][0]>=oldestPos) break;
            firstActiveSnp=(firstActiveSnp+1)%activeSnps.length;
            if (lastReached) break;
            i=(i+1)%activeSnps.length;
        }
        if (firstActiveSnp>lastActiveSnp) { firstActiveSnp=-1; lastActiveSnp=-1; }
    }


    private static final void computeMaxRsquare(int svIndex) {
        boolean lastReached;
        int i;
        double r, max;
        
        i=firstActiveSnp; lastReached=false; max=0.0;
        while (true) {
            if (i==lastActiveSnp) lastReached=true;
            r=rSquare(svIndex,i,tmpArray1,tmpArray2);
            if (r>=max) max=r;
            if (lastReached) break;
            i=(i+1)%activeSnps.length;
        }
        rSquared[svIndex]=max;
    }


    /**
     * @return the R^2 coefficient between `activeSvs[svIndex]` and 
     * `activeSnps[snpIndex]`.
     */
    private static final double rSquare(int svIndex, int snpIndex, boolean[] svMatches, boolean[] snpMatches) {
        final int LAST_SV = activeSvsLast[svIndex];
        final int LAST_SNP = activeSnpsLast[snpIndex];
        int i, j;
        double n, d1, d2, denominator, numerator, avgSv, avgSnp;

        // 1. Averages
        avgSv=0.0;
        for (i=3; i<=LAST_SV; i+=2) avgSv+=activeSvs[svIndex][i];
        avgSv/=((LAST_SV+1-2)/2);
        avgSnp=0.0;
        for (j=3; j<=LAST_SNP; j+=2) avgSnp+=activeSnps[snpIndex][j];
        avgSnp/=((LAST_SNP+1-2)/2);

        // 2. Denominator
        d1=0.0; d2=0.0;
        for (i=3; i<=LAST_SV; i+=2) { n=activeSvs[svIndex][i+1]-avgSv; d1+=n*n; }
        for (i=3; i<=LAST_SNP; i+=2) { n=activeSnps[snpIndex][i+1]-avgSnp; d2+=n*n; }
        denominator=Math.sqrt(d1*d2);

        // 3. Numerator
        numerator=0.0;

        // 3.1 Matches (merge sort)
        for (i=0; i<=LAST_SV; i++) svMatches[i]=false;
        for (i=0; i<=LAST_SNP; i++) snpMatches[i]=false;
        i=2; j=2;
        while (i<=LAST_SV && j<=LAST_SNP) {
            if (activeSvs[svIndex][i]<activeSnps[snpIndex][j]) i+=2;
            else if (activeSvs[svIndex][i]>activeSnps[snpIndex][j]) j+=2;
            else {
                svMatches[i]=true; snpMatches[j]=true;
                numerator+=(activeSvs[svIndex][i+1]-avgSv)*(activeSnps[snpIndex][j+1]-avgSnp);
                i+=2; j+=2;
            }
        }

        // 3.2 Mismatches
        for (i=2; i<=LAST_SV; i+=2) {
            if (!svMatches[i]) numerator+=(activeSvs[svIndex][i+1]-avgSv)*(0.0-avgSnp);
        }
        for (i=2; i<=LAST_SNP; i+=2) {
            if (!snpMatches[i]) numerator+=(0.0-avgSv)*(activeSnps[snpIndex][i+1]-avgSnp);
        }

        n=numerator/denominator;
        return n*n;
    }
    
}