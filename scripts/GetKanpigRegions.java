import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Uses the NE field printed by kanpig to build a BED of kanpig regions.
 */
public class GetKanpigRegions {
    
    /**
     * @param args 0 a VCF whose FORMAT field is assumed to be:
     * GT:FT:SQ:GQ:PS:NE:DP:AD:KS
     */
    public static void main(String[] args) throws IOException {
        final String KANPIG_VCF_GZ = args[0];
        final int MAX_REGION_LENGTH = Integer.parseInt(args[1]);
        
        final int QUANTUM = 10000;  // Arbitrary
        
        int pos, region, refLength, altLength, start, end, currentRegion, currentStart, currentEnd, nRecords;
        String str, chr, currentChr;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(KANPIG_VCF_GZ))));
        str=br.readLine(); nRecords=0; 
        currentRegion=-1; currentChr=""; currentStart=-1; currentEnd=-1;
        while (str!=null) { 
            if (str.charAt(0)=='#') { str=br.readLine(); continue; }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records...");
            tokens=str.split("\t");
            tokensPrime=tokens[9].split(":");
            region=Integer.parseInt(tokensPrime[5]);
            chr=tokens[0];
            pos=Integer.parseInt(tokens[1]);
            refLength=tokens[3].length();
            altLength=tokens[4].length();
            if (refLength==1 && altLength>1) {
                // INS
                start=pos-1;  // Zero-based, inclusive.
                end=(pos+1)-1;  // Zero-based, inclusive.
            }
            else if (refLength>1 && altLength==1) {
                // DEL
                start=(pos+1)-1;
                end=pos+(refLength-1)-1;
            }
            else {
                // Replacement
                start=(pos+1)-1;
                end=pos+(refLength-1)-1;
            }
            if (region!=currentRegion) {
                if (currentRegion!=-1 && currentEnd-currentStart+1<=MAX_REGION_LENGTH) System.out.println(currentChr+"\t"+currentStart+"\t"+(currentEnd+1));
                currentRegion=region; currentChr=chr; currentStart=start; currentEnd=end;
            }
            else {
                if (end>currentEnd) currentEnd=end;
            }
            str=br.readLine();
        }
        System.out.println(currentChr+"\t"+currentStart+"\t"+(currentEnd+1));
        br.close();
    }
    
}