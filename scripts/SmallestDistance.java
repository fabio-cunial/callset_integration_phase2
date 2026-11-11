import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Computes the smallest distance between every long call and every other long
 * call located within a max distance from its endpoints. This is essentially
 * the same algorithm as in `Rsquare.java`.
 * 
 * Remark: for speed, the program works on a small TSV projection of the
 * original VCF, where every record comes from the same chromosome, and
 * records are sorted by POS.
 */
public class SmallestDistance {
    /**
     * Variant types
     */
    private static final int TYPE_DEL = 0;
    private static final int TYPE_INS = 1;
    private static final int TYPE_DUP = 2;
    private static final int TYPE_INV = 3;
    private static final int N_TYPES = 4;

    /**
     * Distance types:
     * 0,1: distance between first positions or last positions;
     *   2: gap between two non-overlapping intervals; `Integer.MAX_VALUE` if
     *      the intervals overlap;
     */
    private static final int DISTANCE_FIRST = 0;
    private static final int DISTANCE_LAST = 1;
    private static final int DISTANCE_GAP = 2;

    /**
     * Array size parameters
     */
    private static final int CAPACITY = 100;  // Arbitrary
    private static final double RESIZE_FACTOR = 1.5;  // Arbitrary

    /**
     * Format of each row:
     * 
     * pos,first,last,type,length
     * 
     * where:
     * - `pos` is the one-based POS from the VCF;
     * - `[first..last]` is the interval of zero-based, inclusive reference 
     *   positions that are affected by the variant.
     */
    private static int[][] activeLong;
    private static int[] nActiveLongs;  // By SVTYPE
    private static int firstActiveLong, lastActiveLong;
    private static String[][] activeLongIds;  // The ID of every active long, and the ID of its closest long of each type (any,del,ins,dup,inv).
    private static String chrId;

    /**
     * Format of each row:
     * 
     * to_any,to_any_overlaps,  to_del,to_del_overlaps,  to_ins,to_ins_overlaps,  to_dup,to_dup_overlaps,  to_inv,to_inv_overlaps
     */
    private static int[][] smallestDistanceLong;


    /**
     * @param args 0 format:
     * 
     * POS, SVLEN, AF, ID, SVTYPE 
     * 
     * where `POS` is sorted and comes from the VCF;
     * @param args 1 the chromosome all records come from; must be one of the
     * autosomes (no chrX, chrY, chrM);
     * @param args 2 only variants of length >=this are used;
     * @param args 3 only variants with AF >=this are used;
     * @param args 4 max bp distance to compare two variants;
     * @param args 5 one of the DISTANCE_* constants;
     * @param args 6 output BED; see `writeToBed()` for format.
     * 
     * Remark: output records are not necessarily sorted by their first 
     * position.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV_GZ = args[0];
        chrId=args[1];
        final int MIN_LONG_LENGTH = Integer.parseInt(args[2]);
        final double MIN_AF = Double.parseDouble(args[3]);
        final int MAX_DISTANCE_BP = Integer.parseInt(args[4]);
        final int DISTANCE_TYPE = Integer.parseInt(args[5]);
        final String OUTPUT_BED = args[6];

        final int QUANTUM = 10000;  // Arbitrary

        int p1, p2, p3, p4;
        int pos, svtype, svlen;
        long nRecords, nLong, nLongWithAF;
        double af;
        String str, id, svtypeStr;
        BufferedReader br;
        BufferedWriter bw;
        
        // Allocating reused space
        activeLong = new int[CAPACITY][5];
        nActiveLongs = new int[N_TYPES];
        firstActiveLong=-1; lastActiveLong=-1;
        activeLongIds = new String[CAPACITY][N_TYPES+2];
        smallestDistanceLong = new int[CAPACITY][(N_TYPES+1)*2];

        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_TSV_GZ))));
        bw = new BufferedWriter(new FileWriter(OUTPUT_BED));
        str=br.readLine(); nRecords=0; nLong=0; nLongWithAF=0;
        while (str!=null) {
            p1=str.indexOf('\t');
            pos=Integer.parseInt(str.substring(0,p1));
            p2=str.indexOf('\t',p1+1);
            svlen=Math.abs(Integer.parseInt(str.substring(p1+1,p2)));
            if (svlen>=MIN_LONG_LENGTH) {
                nLong++;
                p3=str.indexOf('\t',p2+1);
                af=Double.parseDouble(str.substring(p2+1,p3));
                if (af>=MIN_AF) {
                    nLongWithAF++;
                    p4=str.indexOf('\t',p3+1);
                    id=str.substring(p3+1,p4);
                    svtypeStr=str.substring(p4+1);
                    if (svtypeStr.equalsIgnoreCase("DEL")) svtype=TYPE_DEL;
                    else if (svtypeStr.equalsIgnoreCase("INS")) svtype=TYPE_INS;
                    else if (svtypeStr.equalsIgnoreCase("DUP")) svtype=TYPE_DUP;
                    else if (svtypeStr.equalsIgnoreCase("INV")) svtype=TYPE_INV;
                    else {
                        System.err.println("ERROR: wrong type="+svtypeStr+" in string: "+str);
                        System.exit(1);
                        svtype=-1;
                    }
                    removeInactiveLongs(pos-1-MAX_DISTANCE_BP-1,bw);
                    addActiveLong(pos,svlen,id,svtype);
                    initializeSmallestDistance(MAX_DISTANCE_BP,DISTANCE_TYPE);
                    updateSmallestDistance(MAX_DISTANCE_BP,DISTANCE_TYPE);
                }
            }

            // Next iteration
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            str=br.readLine();
        }
        br.close();
        printRemainingLongs(bw);
        System.out.println("Processed "+nRecords+" total records, "+nLong+" long ("+nLongWithAF+" with the desired AF).");
    }


    /**
     * @param pos one-based, from the VCF.
     */
    private static final void addActiveLong(int pos, int svlen, String id, int svtype) {
        int i;

        // activeLong
        if (lastActiveLong==-1) { firstActiveLong=0; lastActiveLong=0; }
        else {
            if ((lastActiveLong+1)%activeLong.length==firstActiveLong) resizeActiveLongs();
            lastActiveLong=(lastActiveLong+1)%activeLong.length;
        }
        if (activeLong[lastActiveLong]==null) activeLong[lastActiveLong] = new int[5];        
        activeLong[lastActiveLong][0]=pos;
        getInterval(pos,svlen,svtype,activeLong[lastActiveLong],1);
        activeLong[lastActiveLong][3]=svtype;
        activeLong[lastActiveLong][4]=svlen;
        nActiveLongs[svtype]++;

        // activeLongIds
        if (activeLongIds[lastActiveLong]==null) activeLongIds[lastActiveLong] = new String[N_TYPES+2];
        activeLongIds[lastActiveLong][0]=id; 
        for (i=1; i<activeLongIds[lastActiveLong].length; i++) activeLongIds[lastActiveLong][i]=".";

        // smallestDistanceLong
        if (smallestDistanceLong[lastActiveLong]==null) smallestDistanceLong[lastActiveLong] = new int[(N_TYPES+1)*2];
        for (i=0; i<N_TYPES+1; i++) {
            smallestDistanceLong[lastActiveLong][2*i]=Integer.MAX_VALUE;
            smallestDistanceLong[lastActiveLong][2*i+1]=0;
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
        int i;

        if (firstActiveLong==-1) return;
        lastReached=false;
        while (true) {
            if (firstActiveLong==lastActiveLong) lastReached=true;
            if (activeLong[firstActiveLong][2]>=pos) break;
            writeToBed(firstActiveLong,bw);
            nActiveLongs[activeLong[firstActiveLong][3]]--;
            for (i=0; i<activeLongIds[firstActiveLong].length; i++) activeLongIds[firstActiveLong][i]=null;
            if (lastReached) { firstActiveLong=-1; lastActiveLong=-1; break; }
            else firstActiveLong=(firstActiveLong+1)%activeLong.length;   
        }
    }


    /**
     * Initializes `smallestDistanceLong[lastActiveLong]` using all the active 
     * longs that are at distance `<=threshold` from `lastActiveLong` (if any).
     * 
     * Remark: if `getDistance()` is monotonically non-decreasing with 
     * decreasing POS, the procedure stops as soon as min distances to every
     * SVTYPE have been computed.
     */
    private static final void initializeSmallestDistance(int threshold, int distanceType) {
        final int OUT_OF_RANGE = (firstActiveLong-1+activeLong.length)%activeLong.length;
        final boolean EARLY_QUIT = distanceType==DISTANCE_FIRST;
        boolean allDistancesComputed, overlapping;
        int i, j;
        int distance, svtype;
        
        i=(lastActiveLong-1+activeLong.length)%activeLong.length;
        while (true) {
            if (i==OUT_OF_RANGE) break;
            distance=getDistance(distanceType,lastActiveLong,i);
            if (EARLY_QUIT) {
                if (distance>threshold) break;
                allDistancesComputed=true;
                for (j=0; j<N_TYPES; j++) {
                    if ( ( (j==activeLong[lastActiveLong][3] && nActiveLongs[j]>1) || 
                           (j!=activeLong[lastActiveLong][3] && nActiveLongs[j]>0)
                         ) &&
                         smallestDistanceLong[lastActiveLong][2+2*j]==Integer.MAX_VALUE
                       ) { allDistancesComputed=false; break; }
                }
                if (allDistancesComputed) break;
            }
            if (distance<=threshold) {
                overlapping=areOverlapping(lastActiveLong,i);
                if (distance<smallestDistanceLong[lastActiveLong][0]) {
                    smallestDistanceLong[lastActiveLong][0]=distance;
                    smallestDistanceLong[lastActiveLong][1]=overlapping?1:0;
                    activeLongIds[lastActiveLong][1]=activeLongIds[i][0];
                }
                svtype=activeLong[i][3];
                if (distance<smallestDistanceLong[lastActiveLong][2+svtype*2]) {
                    smallestDistanceLong[lastActiveLong][2+svtype*2]=distance;
                    smallestDistanceLong[lastActiveLong][2+svtype*2+1]=overlapping?1:0;
                    activeLongIds[lastActiveLong][2+svtype]=activeLongIds[i][0];
                }
            }
            i=(i-1+activeLong.length)%activeLong.length;
        }
    }


    /**
     * Uses the last active long to update the `smallestDistanceLong` row of
     * every active long that is at distance `<=threshold` from it (if any).
     * 
     * Remark: if `getDistance()` is monotonically non-decreasing with 
     * decreasing POS, the procedure stops as soon as another call of the same
     * type is encountered.
     */
    private static final void updateSmallestDistance(int threshold, int distanceType) {
        final int OUT_OF_RANGE = (firstActiveLong-1+activeLong.length)%activeLong.length;
        final boolean EARLY_QUIT = distanceType==DISTANCE_FIRST;
        final int SVTYPE = activeLong[lastActiveLong][3];
        boolean overlapping;
        int i;
        int distance;
        
        i=(lastActiveLong-1+activeLong.length)%activeLong.length;
        while (true) {
            if (i==OUT_OF_RANGE) break;
            distance=getDistance(distanceType,i,lastActiveLong);
            if (distance<=threshold) {
                overlapping=areOverlapping(i,lastActiveLong);
                if (distance<smallestDistanceLong[i][0]) {
                    smallestDistanceLong[i][0]=distance;
                    smallestDistanceLong[i][1]=overlapping?1:0;
                    activeLongIds[i][1]=activeLongIds[lastActiveLong][0];
                }
                if (distance<smallestDistanceLong[i][2+SVTYPE*2]) {
                    smallestDistanceLong[i][2+SVTYPE*2]=distance;
                    smallestDistanceLong[i][2+SVTYPE*2+1]=overlapping?1:0;
                    activeLongIds[i][2+SVTYPE]=activeLongIds[lastActiveLong][0];
                }
                if (activeLong[i][3]==SVTYPE && EARLY_QUIT) break;
            }
            else if (EARLY_QUIT) break;
            i=(i-1+activeLong.length)%activeLong.length;
        }
    }


    private static final int getDistance(int distanceType, int index1, int index2) {
        int distance;

        if (distanceType==DISTANCE_FIRST) return Math.abs(activeLong[index1][1]-activeLong[index2][1]);
        else if (distanceType==DISTANCE_LAST) return Math.abs(activeLong[index1][2]-activeLong[index2][2]);
        else if (distanceType==DISTANCE_GAP) {
            // Special case: two INS at the same POS are considered overlapping.
            if (activeLong[index1][2]<activeLong[index1][1] && activeLong[index2][2]<activeLong[index2][1] && activeLong[index1][1]==activeLong[index2][1]) return Integer.MAX_VALUE;
            // Every other pair of calls.
            // Remark: an INS overlaps with any interval call that contains or is 
            // identical to its two adjacent positions.
            distance=Math.max(activeLong[index2][1]-activeLong[index1][2]-1,activeLong[index1][1]-activeLong[index2][2]-1);
            return distance<0?Integer.MAX_VALUE:distance;
        }
        else return Integer.MAX_VALUE;
    }


    public static final boolean areOverlapping(int index1, int index2) {
        // Special case: two INS at the same POS are considered overlapping.
        if (activeLong[index1][2]<activeLong[index1][1] && activeLong[index2][2]<activeLong[index2][1] && activeLong[index1][1]==activeLong[index2][1]) return true;
        // Every other pair of calls.
        // Remark: an INS overlaps with any interval call that contains or is 
        // identical to its two adjacent positions.
        final int distance = Math.max(activeLong[index2][1]-activeLong[index1][2]-1,activeLong[index1][1]-activeLong[index2][2]-1);
        return distance<0;
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
        int i;

        bw.write( chrId+"\t" +                              // CHROM
                  // first
                  (activeLong[longIndex][3]==TYPE_INS?activeLong[longIndex][1]-1:activeLong[longIndex][1])+"\t" +
                  // last+1, avoiding empty intervals that may interfere with bedtools intersect.
                  (activeLong[longIndex][3]==TYPE_INS?activeLong[longIndex][1]:activeLong[longIndex][2]+1)+"\t" +
                  activeLongIds[longIndex][0]+"\t" +        // ID
                  activeLong[longIndex][3]+"\t" +           // type
                  activeLong[longIndex][4]                  // length
                );
        for (i=0; i<N_TYPES+1; i++) bw.write("\t"+(smallestDistanceLong[longIndex][2*i]==Integer.MAX_VALUE?"-1":smallestDistanceLong[longIndex][2*i])+"\t"+smallestDistanceLong[longIndex][2*i+1]);
        for (i=0; i<N_TYPES+1; i++) bw.write("\t"+activeLongIds[longIndex][1+i]);
        bw.newLine();
    }


    /**
     * Loads in `out[outFirst..outFirst+1]` the interval [first..last] (zero-
     * based, inclusive) that corresponds to all and only the positions affected
     * by the variant.
     * 
     * Remark: the interval is empty for INS. This is represented by setting
     * `last=first-1`.
     */
    private static final void getInterval(int pos, int svlen, int type, int[] out, int outFirst) {
        if (type==TYPE_DEL || type==TYPE_INV || type==TYPE_DUP) { out[outFirst]=(pos-1)+1; out[outFirst+1]=(pos-1)+svlen; }
        else if (type==TYPE_INS) { out[outFirst]=pos; out[outFirst+1]=pos-1; }
        else {
            System.err.println("ERROR: wrong type="+type+" for length "+svlen);
            System.exit(1);
        }
    }


    private static final void resizeActiveLongs() {
        final int NEW_SIZE = (int)(activeLong.length*RESIZE_FACTOR);

        int[][] newArray = new int[NEW_SIZE][];
        System.arraycopy(activeLong,firstActiveLong,newArray,0,activeLong.length-firstActiveLong);
        if (lastActiveLong<activeLong.length-1) System.arraycopy(activeLong,0,newArray,activeLong.length-firstActiveLong,lastActiveLong+1);
        int[][] newArrayDistance = new int[NEW_SIZE][];
        System.arraycopy(smallestDistanceLong,firstActiveLong,newArrayDistance,0,smallestDistanceLong.length-firstActiveLong);
        if (lastActiveLong<smallestDistanceLong.length-1) System.arraycopy(smallestDistanceLong,0,newArrayDistance,smallestDistanceLong.length-firstActiveLong,lastActiveLong+1);
        String[][] newArrayIds = new String[NEW_SIZE][];
        System.arraycopy(activeLongIds,firstActiveLong,newArrayIds,0,activeLongIds.length-firstActiveLong);
        if (lastActiveLong<activeLongIds.length-1) System.arraycopy(activeLongIds,0,newArrayIds,activeLongIds.length-firstActiveLong,lastActiveLong+1);
        firstActiveLong=0; lastActiveLong=activeLong.length-1;
        activeLong=newArray; smallestDistanceLong=newArrayDistance; activeLongIds=newArrayIds;
    }
    
}