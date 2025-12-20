import java.io.*;
import java.util.zip.*;


/**
 * Prints $DP,AD_ref,AD_alt$ for every 0/0 cell of a cohort VCF.
 */
public class ReadsOfHomAlts {
    
    /**
     * Format: GT:FT:SQ:GQ:PS:NE:DP:AD:KS
     * Example: 0|1:22:53:5:.:73810:3:1,2:97
     */
    public static void main(String[] args) throws IOException {
        final String VCF_GZ = args[0];
        
        final char FIELD_SEPARATOR = '\t';
        final char GT_SEPARATOR = ':';
        final int QUANTUM = 1000;  // Arbitrary
        
        int i, p, q;
        int length, nLines, start, end;
        long ms;
        String str;
        BufferedReader br;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(VCF_GZ))));
        str=br.readLine(); 
        while (str.charAt(0)=='#') str=br.readLine();
        nLines=0;
        ms=System.currentTimeMillis();
        while (str!=null) {
            length=str.length();
            start=0;
            for (i=0; i<=8; i++) start=str.indexOf(FIELD_SEPARATOR,start)+1;
            do {
                end=str.indexOf(FIELD_SEPARATOR,start);
                if (end<0) end=length;
                if (str.indexOf(GT_SEPARATOR,start)==start+3) {
                    // Only GTs with two alleles
                    if (str.charAt(start)=='0' && str.charAt(start+2)=='0') {
                        p=start;
                        for (i=0; i<=5; i++) p=str.indexOf(GT_SEPARATOR,p+1);
                        q=str.indexOf(GT_SEPARATOR,p+1);
                        System.out.print(str.substring(p+1,q)+",");
                        p=q+1;
                        q=str.indexOf(GT_SEPARATOR,p+1);
                        if (q>end) q=end;
                        System.out.println(str.substring(p,q));
                    }
                }
                start=end+1;
            } while (start<length);
            nLines++;
            if (nLines%QUANTUM==0) System.err.println("Processed "+nLines+" lines in "+((double)(System.currentTimeMillis()-ms)/1000)+"s");
            str=br.readLine();
        }
        br.close();
    }
    
}
