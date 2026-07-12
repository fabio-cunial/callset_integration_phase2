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
public class BndFilterWithAssemblyBreakpoints2 {

    /**
     * @param args
     * 0: not necessarily sorted;
     * 1: not necessarily sorted; every breakpoint with a mate is assumed to 
     *    have a symmetric record;
     * 3: reference gaps file (AGP); not necessarily sorted.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String BREAKPOINTS_CSV = args[1];
        final int MAX_DISTANCE = Integer.parseInt(args[2]);
        final String REFERENCE_AGP = args[3];
    
        int p, q, r;
        int nRecordsIn, nRecordsOut;
        String str, chrom;
        BufferedReader br;
        String[] tokens;
        Breakpoint query, queryPrime;
        ArrayList<Breakpoint> gaps;
        ArrayList<Breakpoint> breakpoints;

        // Loading all breakpoints.
        // Remark: from `AssemblySam2Breakpoints2.java`, breakpoints without a 
        // mate have otherChr=-1, otherPos=-1.
        breakpoints = new ArrayList<Breakpoint>();
        br = new BufferedReader(new FileReader(BREAKPOINTS_CSV));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            breakpoints.add(new Breakpoint(tokens[0],Integer.parseInt(tokens[1]),tokens[2],Integer.parseInt(tokens[3])));
            str=br.readLine();
        }
        br.close();
        breakpoints.sort(null);

        // Loading all gaps.
        // Remark: gaps are marked with otherChr=-2, otherPos=-2.
        gaps = new ArrayList<Breakpoint>();
        br = new BufferedReader(new FileReader(REFERENCE_AGP));

------> Only load the correct records from this file!!!!!!!! Or pre-filter the file in the WDL.

        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') { str=br.readLine(); continue; }
            p=str.indexOf('\t');
            q=str.indexOf('\t',p+1);
            r=str.indexOf('\t',q+1);
            if (r<0) r=str.length();
            chrom=str.substring(0,p);
            gaps.add(new Breakpoint(chrom,Integer.parseInt(str.substring(p+1,q)),"-2",-2));
            gaps.add(new Breakpoint(chrom,Integer.parseInt(str.substring(q+1,r)),"-2",-2));
            str=br.readLine();
        }
        br.close();
        gaps.sort(null);

        // Filtering the VCF
        query = new Breakpoint("",0,"",0);
        queryPrime = new Breakpoint("",0,"",0);
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecordsIn=0; nRecordsOut=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecordsIn++;
            query.loadFromVcf(str.split("\t")); query.symmetrize(queryPrime);
            if (!find(query,gaps,MAX_DISTANCE) && (query.otherPos>0 && !find(queryPrime,gaps,MAX_DISTANCE)) && find(query,breakpoints,MAX_DISTANCE)) {
                System.out.println(str);
                nRecordsOut++;
            }
            str=br.readLine();
        }
        br.close();
        System.err.println("Kept "+nRecordsOut+" records out of "+nRecordsIn);
    }


    /**
     * A breakpoint without a mate can match anything.
     * A breakpoint with a mate can match only gaps and breakpoints with a mate.
     */
    private static final boolean find(Breakpoint query, ArrayList<Breakpoint> breakpoints, int maxDistance) {
        final int N_BREAKPOINTS = breakpoints.size();
        int i, p;
        Breakpoint current;

        p=Collections.binarySearch(breakpoints,query);
        if (p>=0) {
            if (query.otherPos<0) {
                // Query does not have a mate
                return true;
            }
        }
        else p=-p-1;
        for (i=p-1; i>=0; i--) {
            current=breakpoints.get(i);
            if (!current.chr.equals(query.chr)) break;
            if (current.pos<query.pos-maxDistance) break;
            if (query.otherPos<0) {
                // Query does not have a mate
                return true;
            }
            else {
                // Query has a mate
                if (current.otherPos==-2) {
                    // The current breakpoint is a gap
                    return true;
                }
                else if (current.otherPos>0) {
                    // The current breakpoint has a mate
                    if (Math.abs(current.otherPos-query.otherPos)<=maxDistance && current.otherChr.equals(query.otherChr)) return true;
                }
            }
        }
        for (i=p; i<N_BREAKPOINTS; i++) {
            current=breakpoints.get(i);
            if (!current.chr.equals(query.chr)) break;
            if (current.pos>query.pos+maxDistance) break;
            if (query.otherPos<0) {
                // Query does not have a mate
                return true;
            }
            else {
                // Query has a mate
                if (current.otherPos==-2) {
                    // The current breakpoint is a gap
                    return true;
                }
                else if (current.otherPos>0) {
                    // The current breakpoint has a mate
                    if (Math.abs(current.otherPos-query.otherPos)<=maxDistance && current.otherChr.equals(query.otherChr)) return true;
                }
            }
        }
        return false;
    }


    private static final class Breakpoint implements Comparable<Breakpoint> {
        public String chr, otherChr;
        public int pos, otherPos;

        Breakpoint(String chr, int pos, String otherChr, int otherPos) { 
            this.chr=chr; this.pos=pos; this.otherChr=otherChr; this.otherPos=otherPos;
        }

        public void loadFromVcf(String[] tokens) {
            int p, q, r;
            String alt;

            this.chr=tokens[0];
            this.pos=Integer.parseInt(tokens[1]);
            alt=tokens[4];
            p=alt.indexOf('['); q=alt.indexOf(']'); 
            if (p>=0) {
                r=alt.indexOf(':',p+1);
                this.otherChr=alt.substring(p+1,r);
                this.otherPos=Integer.parseInt(alt.substring(r+1,alt.indexOf('[',r+1)));
            }
            else if (q>=0) {
                r=alt.indexOf(':',q+1);
                this.otherChr=alt.substring(q+1,r);
                this.otherPos=Integer.parseInt(alt.substring(r+1,alt.indexOf(']',r+1)));
            }
            else {
                System.err.println("ERROR: unrecognized ALT: "+alt);
                System.exit(1);
            }
        }
        
        /**
         * Stores in `other` the symmetrized version of this breakpoint.
         */
        public void symmetrize(Breakpoint other) {
            other.chr=this.otherChr; other.pos=this.otherPos;
            other.otherChr=this.chr; other.otherPos=this.pos;
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