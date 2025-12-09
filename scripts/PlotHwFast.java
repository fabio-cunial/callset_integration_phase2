import java.io.*;
import java.util.zip.*;


/**
 * Remark: the program prints only records whose GT has two alleles and that
 * are ALT in some sample. Missing characters in a GT are assumed to be zeros.
 */
public class PlotHwFast {
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String VCF_GZ = args[0];
        final String OUT_HWE = args[1];
        
        final char FIELD_SEPARATOR = '\t';
        final char GT_SEPARATOR = ':';
        final int QUANTUM = 1000;  // Arbitrary
        
        boolean onLeft, onRight;
        int i, p;
        int length, nLines, start, end;
        int gt00, gt01, gt11;
        long ms;
        String str;
        BufferedReader br;
        BufferedWriter bwHwe;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(VCF_GZ))));
        bwHwe = new BufferedWriter(new FileWriter(OUT_HWE));
        bwHwe.write("AA,AB,BB\n");
        str=br.readLine(); nLines=0;
        ms=System.currentTimeMillis();
        while (str!=null) {
            gt00=0; gt01=0; gt11=0;
            length=str.length();
            start=0;
            for (i=0; i<=8; i++) start=str.indexOf(FIELD_SEPARATOR,start)+1;
            do {
                end=str.indexOf(FIELD_SEPARATOR,start);
                if (end<0) end=length;
                p=str.indexOf(GT_SEPARATOR,start);
                if (p==start+3) {  // Only GTs with two alleles
                    onLeft=str.charAt(start)=='1';
                    onRight=str.charAt(start+2)=='1';
                    if (onLeft==onRight) {
                        if (onLeft) gt11++;
                        else gt00++;
                    }
                    else gt01++;
                }
                start=end+1;
            } while (start<length);
            if (gt01+gt11>0) {  // Only records that occur in some sample
                bwHwe.write(gt00+","+gt01+","+gt11+"\n");
            }
            nLines++;
            if (nLines%QUANTUM==0) System.err.println("Processed "+nLines+" lines in "+((double)(System.currentTimeMillis()-ms)/1000)+"s");
            str=br.readLine();
        }
        br.close(); bwHwe.close();
    }
    
}