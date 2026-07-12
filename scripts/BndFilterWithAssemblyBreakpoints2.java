import java.util.ArrayList;
import java.util.Collections;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Keeps the records of a BND-only VCF that are close enough to a given set of
 * true breakpoints, without being close to reference gaps.
 * 
 * Remarks:
 * - the orientation of the BNDs is not taken into account; this would be easy
 *   to implement;
 * - for simplicity, BNDs are assumed to follow the simple form without
 *   inserted sequence;
 * - the implementation could be made much faster;
 * - all the internal matching logic works on coordinates in 0-based form.
 */
public class BndFilterWithAssemblyBreakpoints2 {

    /**
     * @param args
     * 0: not necessarily sorted; every breakend if assumed to have two sides
     *    and to be represented by exactly one record (rather than by two
     *    symmetric records);
     * 1: not necessarily sorted; every breakend with two sides is assumed to 
     *    be represented by two symmetric records;
     * 3: reference gaps file (AGP); not necessarily sorted.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String BREAKPOINTS_CSV = args[1];
        final int MAX_DISTANCE = Integer.parseInt(args[2]);
        final String REFERENCE_AGP = args[3];
    
        int i, p;
        int nRecordsIn, nRecordsOut;
        String str;
        BufferedReader br;
        Breakpoint query, queryPrime;
        String[] tokens;
        ArrayList<Breakpoint> gaps;
        ArrayList<Breakpoint> breakpoints;

        // Loading all breakpoints.
        // Remark: from `AssemblySam2Breakpoints2.java`, breakpoints without a 
        // mate have otherChr=-1, otherPos=-1, type=4 or 5.
        breakpoints = new ArrayList<Breakpoint>();
        br = new BufferedReader(new FileReader(BREAKPOINTS_CSV));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            breakpoints.add(new Breakpoint(tokens[0],Integer.parseInt(tokens[1]),tokens[2],Integer.parseInt(tokens[3]),Integer.parseInt(tokens[4])));  // Already zero-based
            str=br.readLine();
        }
        br.close();
        breakpoints.sort(null);

        // Loading all gaps
        // Remarks: 
        // - AGP files are 1-based inclusive;
        // - we mark gaps with otherChr=-2, otherPos=-2, type=-1.
        gaps = new ArrayList<Breakpoint>();
        br = new BufferedReader(new FileReader(REFERENCE_AGP));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') { str=br.readLine(); continue; }
            tokens=str.split("\t");
            if (tokens[4].equals("N") || tokens[4].equals("U")) {
                gaps.add(new Breakpoint(tokens[0],Integer.parseInt(tokens[1])-1,"-2",-2,-1));
                gaps.add(new Breakpoint(tokens[0],Integer.parseInt(tokens[2])-1,"-2",-2,-1));
            }
            str=br.readLine();
        }
        br.close();
        gaps.sort(null);

        // Filtering the VCF
        query = new Breakpoint("",0,"",0,-2);
        queryPrime = new Breakpoint("",0,"",0,-2);
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecordsIn=0; nRecordsOut=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.startsWith("#CHROM")) System.out.println("##INFO=<ID=TRUTH_BND_TYPE,Number=1,Type=Integer,Description=\"Type of assembly-derived BND that matches with this BND\">");
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecordsIn++;
            tokens=str.split("\t");
            query.loadFromVcf(tokens); query.symmetrize(queryPrime);
            if (find(query,gaps,MAX_DISTANCE)==-1 && find(queryPrime,gaps,MAX_DISTANCE)==-1) {
                p=find(query,breakpoints,MAX_DISTANCE);
                if (p<0) p=find(queryPrime,breakpoints,MAX_DISTANCE);
                if (p>=0) {
                    tokens[7]+=";TRUTH_BND_TYPE="+breakpoints.get(p).type;
                    System.out.println(String.join("\t",tokens));
                    nRecordsOut++;
                }
            }
            str=br.readLine();
        }
        br.close();
        System.err.println("Kept "+nRecordsOut+" records out of "+nRecordsIn);
    }


    /**
     * Remark: this procedure focuses on the first side of `query`. If no match
     * is found, the caller should try again with a symmetrized version of 
     * `query`.
     * 
     * @param query assumed to always be of mated type;
     * @return the index of a best matching breakpoint in `breakpoints`, or -1
     * if no match is found; mated matches whose mate also matches the query are
     * prioritized over unmated matches that match the query only on one side. 
     * No ranking is enforced over the breakpoint types defined in
     * `AssemblySam2Breakpoints2.java`.
     */
    private static final int find(Breakpoint query, ArrayList<Breakpoint> breakpoints, int maxDistance) {
        final int N_BREAKPOINTS = breakpoints.size();
        int i, p;
        int out;
        Breakpoint current;

        p=Collections.binarySearch(breakpoints,query);
        if (p<0) p=-p-1;
        out=-1;
        for (i=p-1; i>=0; i--) {
            current=breakpoints.get(i);
            if (!current.chr.equals(query.chr)) break;
            if (current.pos<query.pos-maxDistance) break;
            if (current.otherPos<0) {
                if (out==-1) out=i;
            }
            else if (Math.abs(current.otherPos-query.otherPos)<=maxDistance && current.otherChr.equals(query.otherChr)) {
                if (out==-1 || breakpoints.get(out).otherPos<0) out=i;
            }
        }
        for (i=p; i<N_BREAKPOINTS; i++) {
            current=breakpoints.get(i);
            if (!current.chr.equals(query.chr)) break;
            if (current.pos>query.pos+maxDistance) break;
            if (current.otherPos<0) {
                if (out==-1) out=i;
            }
            else if (Math.abs(current.otherPos-query.otherPos)<=maxDistance && current.otherChr.equals(query.otherChr)) {
                if (out==-1 || breakpoints.get(out).otherPos<0) out=i;
            }
        }
        return out;
    }


    private static final class Breakpoint implements Comparable<Breakpoint> {
        public String chr, otherChr;
        public int pos, otherPos;

        /**
         * Decided by `AssemblySam2Breakpoints2.java`:
         * 
         * 1=violation1, 2=violation2, 3=violation3, 4=first/last, 5=internal.
         * 
         * For completeness we also set:
         * 
         * 0=input VCF breakpoint, -1=reference gap.
         */
        public int type;

        Breakpoint(String chr, int pos, String otherChr, int otherPos, int type) { 
            this.chr=chr; this.pos=pos; this.otherChr=otherChr; this.otherPos=otherPos; this.type=type;
        }

        public void loadFromVcf(String[] tokens) {
            int p, q, r;
            String alt;

            this.chr=tokens[0];
            this.pos=Integer.parseInt(tokens[1])-1;
            alt=tokens[4];
            p=alt.indexOf('['); q=alt.indexOf(']'); 
            if (p>=0) {
                r=alt.indexOf(':',p+1);
                this.otherChr=alt.substring(p+1,r);
                this.otherPos=Integer.parseInt(alt.substring(r+1,alt.indexOf('[',r+1)))-1;
            }
            else if (q>=0) {
                r=alt.indexOf(':',q+1);
                this.otherChr=alt.substring(q+1,r);
                this.otherPos=Integer.parseInt(alt.substring(r+1,alt.indexOf(']',r+1)))-1;
            }
            else {
                System.err.println("ERROR: unrecognized ALT: "+alt);
                System.exit(1);
            }
            this.type=0;
        }
        
        /**
         * Stores in `other` the symmetrized version of this breakpoint.
         */
        public void symmetrize(Breakpoint other) {
            other.chr=this.otherChr; other.pos=this.otherPos;
            other.otherChr=this.chr; other.otherPos=this.pos;
            other.type=this.type;
        }

        /**
         * Sorts by `chr,pos` only, disregarding the mate.
         */
        public int compareTo(Breakpoint other) {
            int c = this.chr.compareTo(other.chr);
            if (c!=0) return c;
            return this.pos-other.pos;
        }
    }

}