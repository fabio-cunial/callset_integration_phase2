import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Resets ALT to <DEL-ID>, regardless of the actual SVTYPE. This is necessary to
 * make `ExtractVariantAnnotations --resource-matching-strategy START_POSITION_
 * AND_GIVEN_REPRESENTATION` work correctly downstream.
 * 
 * Remark: to make `ExtractVariantAnnotations` work, it is not necessary to set 
 * SVTYPE=DEL in INFO.
 * 
 * Remark: simply setting ALT=ID makes `ExtractVariantAnnotations` fail, because
 * it assumes the record is a BND and it claims that this type is not supported.
 */
public class SetAltToID {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        boolean isEven1, isEven2;
        int nRecords, parityType;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader( new InputStreamReader( INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz") ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.startsWith("#CHROM")) System.out.println("##INFO=<ID=BND_PARITY_TYPE,Number=1,Type=String,Description=\"0=even-even, 1=even-odd, 2=odd-odd.\">");
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");

            // Computing parity
            isEven1=isChromEven(tokens[0]);
            isEven2=isChromEven(getAltChrom(tokens[4]));
            if (isEven1 && isEven2) parityType=0;
            else if ((!isEven1) && (!isEven2)) parityType=2;
            else parityType=1;
            tokens[7]+=";BND_PARITY_TYPE="+parityType;

            // Updating ALT
            tokens[4]="<DEL-"+tokens[2]+">";

            System.out.println(String.join("\t",tokens));
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
    }


    private static String getAltChrom(String alt) {
        char separator;
        int p, q;
        int first;
        
        first=-1; separator='_';
        p=alt.indexOf('['); q=alt.indexOf(']');
        if (p>=0) { separator='['; first=p; }
        else if (q>=0) { separator=']'; first=q; }
        else {
            System.err.println("ERROR: Unknown ALT format: "+alt);
            System.exit(1);
        }
        if (p>1 || q>1) {
            System.err.println("ERROR: Unknown ALT format: "+alt);
            System.exit(1);
        }
        p=alt.indexOf(':',first+1);
        return alt.substring(first+1,p);
    }


    private static final boolean isChromEven(String chrom) {
        if (chrom.charAt(3)=='X' || chrom.charAt(3)=='x') return false;
        else if (chrom.charAt(3)=='Y' || chrom.charAt(3)=='y') return true;
        else if (chrom.charAt(3)=='M' || chrom.charAt(3)=='m') return false;
        else return Integer.parseInt(chrom.substring(3))%2==0;
    }

}