import java.util.ArrayList;
import java.util.Collections;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Keeps the records of a BND-only VCF that are close enough to a given set of
 * breakpoints, without being a reference gap.
 * 
 * Remarks:
 * - the orientation of the BNDs is not taken into account;
 * - for simplicity, BNDs are assumed to follow the simple form without 
 *   inserted sequence;
 * - the implementation could be made much faster.
 */
public class BndFilterWithAssemblyBreakpoints {

    /**
     * @param args
     * 0: not necessarily sorted;
     * 1: not necessarily sorted;
     * 3: reference gaps file (AGP); not necessarily sorted;
     * 4: 2=both sides of the BND must be close to a breakpoint; 1=at least one
     *    side of the BND must be close to a breakpoint.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String BREAKPOINTS_CSV = args[1];
        final int MAX_DISTANCE = Integer.parseInt(args[2]);
        final String REFERENCE_AGP = args[3];
        final int MODE = Integer.parseInt(args[4]);
    
        int p, q, r;
        int nRecordsIn, nRecordsOut;
        String str, chrom;
        BufferedReader br;
        String[] tokens;
        Breakpoint tmp;
        ArrayList<Breakpoint> gaps;
        ArrayList<Breakpoint> breakpoints;

        // Loading all breakpoints
        breakpoints = new ArrayList<Breakpoint>();
        br = new BufferedReader(new FileReader(BREAKPOINTS_CSV));
        str=br.readLine();
        while (str!=null) {
            p=str.indexOf(',');
            q=str.indexOf(',',p+1);
            if (q<0) q=str.length();
            breakpoints.add(new Breakpoint(str.substring(0,p),Integer.parseInt(str.substring(p+1,q))));
            str=br.readLine();
        }
        br.close();
        breakpoints.sort(null);

        // Loading all gaps
        gaps = new ArrayList<Breakpoint>();
        br = new BufferedReader(new FileReader(REFERENCE_AGP));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') { str=br.readLine(); continue; }
            p=str.indexOf('\t');
            q=str.indexOf('\t',p+1);
            r=str.indexOf('\t',q+1);
            if (r<0) r=str.length();
            chrom=str.substring(0,p);
            gaps.add(new Breakpoint(chrom,Integer.parseInt(str.substring(p+1,q))));
            gaps.add(new Breakpoint(chrom,Integer.parseInt(str.substring(q+1,r))));
            str=br.readLine();
        }
        br.close();
        gaps.sort(null);

        // Filtering the VCF
        tmp = new Breakpoint("",0);
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecordsIn=0; nRecordsOut=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecordsIn++;
            tokens=str.split("\t");
            if (keepVcfRecord(tokens,breakpoints,gaps,MAX_DISTANCE,MODE,tmp)) {
                System.out.println(str);
                nRecordsOut++;
            }
            str=br.readLine();
        }
        br.close();
        System.err.println("Kept "+nRecordsOut+" records out of "+nRecordsIn);
    }


    /**
     * Keeps a VCF record iff its REF or ALT is close enough to a breakpoint.
     * 
     * @param mode 2=both sides of the BND must be close to a breakpoint; 
     * 1=at least one side of the BND must be close to a breakpoint;
     * @param tmp temporary space.
     */
    private static final boolean keepVcfRecord(String[] tokens, ArrayList<Breakpoint> breakpoints, ArrayList<Breakpoint> gaps, int maxDistance, int mode, Breakpoint tmp) {
        boolean refFound;
        char separator;
        int p, q;
        int first, pos, altPos;
        String alt, altChrom;

        // REF
        pos=Integer.parseInt(tokens[1]);
        refFound = !find(tokens[0],pos,gaps,maxDistance,tmp) && find(tokens[0],pos,breakpoints,maxDistance,tmp);
        if (mode==1 && refFound) return true;
        else if (mode==2 && !refFound) return false;

        // ALT
        alt=tokens[4];
        p=alt.indexOf('['); q=alt.indexOf(']'); first=-1; separator='_';
        if (p>=0) { separator='['; first=p; }
        else if (q>=0) { separator=']'; first=q; }
        else {
            System.err.println("ERROR: unrecognized ALT = "+alt);
            System.exit(1);
        }
        p=alt.indexOf(':',first+1);
        altChrom=alt.substring(first+1,p);
        q=alt.indexOf(separator,p+1);
        altPos=Integer.parseInt(alt.substring(p+1,q));
        return !find(altChrom,altPos,gaps,maxDistance,tmp) && find(altChrom,altPos,breakpoints,maxDistance,tmp);
    }


    /**
     * @param tmp temporary space;
     * @return TRUE iff the input is at <=maxDistance from a breakpoint.
     */
    private static final boolean find(String chrom, int pos, ArrayList<Breakpoint> breakpoints, int maxDistance, Breakpoint tmp) {
        final int N_BREAKPOINTS = breakpoints.size();
        int i, p;
        Breakpoint breakpoint;

        tmp.chr=chrom; tmp.pos=pos;
        p=Collections.binarySearch(breakpoints,tmp);
        if (p>=0) return true;
        p=-p-1;
        for (i=p-1; i>=0; i--) {
            breakpoint=breakpoints.get(i);
            if (!breakpoint.chr.equals(chrom)) break;
            if (breakpoint.pos<pos-maxDistance) break;
            return true;
        }
        for (i=p; i<N_BREAKPOINTS; i++) {
            breakpoint=breakpoints.get(i);
            if (!breakpoint.chr.equals(chrom)) break;
            if (breakpoint.pos>pos+maxDistance) break;
            return true;
        }
        return false;
    }


    private static final class Breakpoint implements Comparable<Breakpoint> {
        String chr;
        int pos;

        Breakpoint(String chr, int pos) { 
            this.chr=chr; this.pos=pos; 
        }

        public int compareTo(Breakpoint other) {
            int c = this.chr.compareTo(other.chr);
            if (c!=0) return c;
            return this.pos-other.pos;
        }

        @Override
        public boolean equals(Object o) {
            Breakpoint other = (Breakpoint)o;
            return this.pos==other.pos && this.chr.equals(other.chr);
        }
    }

}