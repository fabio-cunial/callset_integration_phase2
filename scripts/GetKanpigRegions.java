import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Uses the NE field printed by kanpig to build a BED of kanpig regions.
 */
public class GetKanpigRegions {
    
    /**
     * Remark: the program handles calls that do not occur in any individual of
     * the first trio.
     *
     * @param args 0 a VCF with >=3 samples and whose FORMAT field is:
     *
     * GT:FT:SQ:GQ:PS:NE:DP:AD:KS
     *
     * If the VCF was annotated with `bcftools +mendelian2 -m a`, the output BED
     * contains 3 additional columns for each region: 
     * 4. total number of calls that are non-missing in every sample and that 
     *    occur in some sample;
     * 5. number of calls above that have a Mendelian error; 
     * 6. total number of calls in the input VCF (no constraints on their GTs).
     */
    public static void main(String[] args) throws IOException {
        final String KANPIG_VCF_GZ = args[0];
        final int MAX_REGION_LENGTH = Integer.parseInt(args[1]);
        
        final int QUANTUM = 10000;  // Arbitrary
        final String MERR_STR = "MERR=";
        final int MERR_STR_LENGTH = MERR_STR.length();
        
        boolean merr, presentChild, presentFather, presentMother, missingChild, missingFather, missingMother;
        int p, q;
        int pos, region, refLength, altLength, start, end, currentRegion, currentStart, currentEnd, nRecords, numerator, denominator, total;
        String str, chr, currentChr;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(KANPIG_VCF_GZ))));
        str=br.readLine(); nRecords=0; 
        currentRegion=-1; currentChr=""; currentStart=-1; currentEnd=-1; numerator=0; denominator=0; total=0;
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
            p=tokens[7].indexOf(MERR_STR);
            if (p>=0) {
                merr=tokens[7].charAt(p+MERR_STR_LENGTH)=='1'?true:false;
            }
            else merr=false;
            p=tokens[9].indexOf(":"); q=tokens[9].indexOf("1"); presentChild=q>=0&&q<p;
            p=tokens[10].indexOf(":"); q=tokens[10].indexOf("1"); presentFather=q>=0&&q<p;
            p=tokens[11].indexOf(":"); q=tokens[11].indexOf("1"); presentMother=q>=0&&q<p;
            p=tokens[9].indexOf(":"); q=tokens[9].indexOf("."); missingChild=q>=0&&q<p;
            p=tokens[10].indexOf(":"); q=tokens[10].indexOf("."); missingFather=q>=0&&q<p;
            p=tokens[11].indexOf(":"); q=tokens[11].indexOf("."); missingMother=q>=0&&q<p;
            if (region!=currentRegion) {
                if (currentRegion!=-1 && currentEnd-currentStart+1<=MAX_REGION_LENGTH && denominator>0) System.out.println(currentChr+"\t"+currentStart+"\t"+(currentEnd+1)+"\t"+denominator+"\t"+numerator+"\t"+total);
                currentRegion=region; currentChr=chr; currentStart=start; currentEnd=end; total=1;
                if ( !missingChild && !missingFather && !missingMother && (presentChild||presentFather||presentMother) ) { numerator=merr?1:0; denominator=1; }
                else { numerator=0; denominator=0; }
            }
            else {
                if (end>currentEnd) currentEnd=end;
                total++;
                if ( !missingChild && !missingFather && !missingMother && (presentChild||presentFather||presentMother) ) {
                    if (merr) numerator++;
                    denominator++;
                }
                else {
                    // NOP
                }
            }
            str=br.readLine();
        }
        if (denominator>0) System.out.println(currentChr+"\t"+currentStart+"\t"+(currentEnd+1)+"\t"+denominator+"\t"+numerator+"\t"+total);
        br.close();
    }
    
}