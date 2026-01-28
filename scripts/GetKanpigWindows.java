import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Uses the PS field printed by kanpig to build a BED of kanpig windows.
 * The output BED contains two additional fields: PS, and number of records with
 * that PS value.
 */
public class GetKanpigWindows {
    
    /**
     * @param args 0 a single-sample VCF whose FORMAT field is:
     *
     * GT:FT:SQ:GQ:PS:DP:AD:KS
     */
    public static void main(String[] args) throws IOException {
        final String KANPIG_VCF_GZ = args[0];
        
        final int QUANTUM = 10000;  // Arbitrary
        
        int p, q;
        int pos, refLength, altLength, start, end, currentStart, currentEnd, nRecords, total;
        String str, chr, currentChr, region, currentRegion;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(KANPIG_VCF_GZ))));
        str=br.readLine(); nRecords=0; 
        currentRegion=""; currentChr=""; currentStart=-1; currentEnd=-1; total=0;
        while (str!=null) { 
            if (str.charAt(0)=='#') { str=br.readLine(); continue; }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records...");
            tokens=str.split("\t");
            tokensPrime=tokens[9].split(":");
            region=tokensPrime[4];
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
            if (!region.equals(currentRegion)) {
                if (currentRegion.length()!=0) System.out.println(currentChr+"\t"+currentStart+"\t"+(currentEnd+1)+"\t"+currentRegion+"\t"+total);
                currentRegion=region; currentChr=chr; currentStart=start; currentEnd=end; total=1;
            }
            else {
                if (end>currentEnd) currentEnd=end;
                total++;
            }
            str=br.readLine();
        }
        System.out.println(currentChr+"\t"+currentStart+"\t"+(currentEnd+1)+"\t"+currentRegion+"\t"+total);
        br.close();
    }
    
}