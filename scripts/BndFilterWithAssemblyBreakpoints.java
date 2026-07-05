import java.util.Vector;
import java.util.Collections;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Keeps the records of a BND-only VCF that are close enough to a given set of
 * breakpoints.
 * 
 * Remarks: 
 * - for simplicity, BNDs are assumed to follow the simple form without 
 *   inserted sequence;
 * - the implementation could be made much faster.
 */
public class BndFilterWithAssemblyBreakpoints {

    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String BREAKPOINTS_CSV = args[1];
        final int MAX_DISTANCE = Integer.parseInt(args[2]);
    
        int p;
        String str;
        BufferedReader br;
        String[] tokens;
        Breakpoint tmp;
        Vector<Breakpoint> breakpoints;

        // Loading all breakpoints
        breakpoints = new Vector<Breakpoint>();
        br = new BufferedReader(new FileReader(BREAKPOINTS_CSV));
        str=br.readLine();
        while (str!=null) {
            p=str.indexOf(',');
            breakpoints.add(new Breakpoint(str.substring(0,p),Integer.parseInt(str.substring(p+1))));
            str=br.readLine();
        }
        br.close();
        breakpoints.sort(null);

        // Filtering the VCF
        tmp = new Breakpoint("",0);
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            if (keepVcfRecord(tokens,breakpoints,MAX_DISTANCE,tmp)) System.out.println(str);
            str=br.readLine();
        }
        br.close();
    }


    /**
     * Keeps a VCF record iff its REF or ALT is close enough to a breakpoint.
     * 
     * @param tmp temporary space.
     */
    private static final boolean keepVcfRecord(String[] tokens, Vector<Breakpoint> breakpoints, int maxDistance, Breakpoint tmp) {
        char separator;
        int p, q;
        int first;
        String alt, altChrom, altPos;

        // REF
        if (isBreakpoint(tokens[0],Integer.parseInt(tokens[1]),breakpoints,maxDistance,tmp)) return true;

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
        altPos=alt.substring(p+1,q);
        return isBreakpoint(altChrom,Integer.parseInt(altPos),breakpoints,maxDistance,tmp);
    }


    /**
     * @param tmp temporary space;
     * @return TRUE iff the input is at <=maxDistance from a breakpoint.
     */
    private static final boolean isBreakpoint(String chrom, int pos, Vector<Breakpoint> breakpoints, int maxDistance, Breakpoint tmp) {
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