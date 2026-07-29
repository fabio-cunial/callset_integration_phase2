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

        int nRecords, parityType, isEven1, isEven2, isBndVcf;
        String str, svtype;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader( new InputStreamReader( INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz") ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0; isBndVcf=-1;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.startsWith("#CHROM")) System.out.println("##INFO=<ID=BND_PARITY_TYPE,Number=1,Type=String,Description=\"0=even-even, 1=even-odd, 2=odd-odd, 3=non-canonical chr.\">");
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");

            // Deciding the type of input calls
            if (isBndVcf==-1) {
                svtype=getInfoField(tokens[7],"SVTYPE");
                if (svtype==null) {
                    System.err.println("ERROR: Unknown SVTYPE: "+str);
                    System.exit(1);
                }
                if (svtype.equals("BND")) isBndVcf=1;
                else isBndVcf=0;
            }

            // Computing BND parity
            if (isBndVcf==1) {
                isEven1=isChromEven(tokens[0]);
                isEven2=isChromEven(getAltChrom(tokens[4]));
                if (isEven1==-1 || isEven2==-1) parityType=3;
                else if (isEven1==1 && isEven2==1) parityType=0;
                else if (isEven1==0 && isEven2==0) parityType=2;
                else parityType=1;
                tokens[7]+=";BND_PARITY_TYPE="+parityType;
            }

            // Updating ALT
            tokens[4]="<DEL-"+tokens[2]+">";

            System.out.println(String.join("\t",tokens));
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
    }


    /**
	 * @return NULL if $field$ does not occur in $info$.
	 */
	private static final String getInfoField(String info, String field) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
        
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return null;
		q=info.indexOf(";",p+FIELD_LENGTH);
		return info.substring(p+FIELD_LENGTH,q<0?info.length():q);
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


    /**
     * @return
     * -1=non-canonical chr;
     *  0=odd;
     *  1=even.
     */
    private static final int isChromEven(String chrom) {
        final int length = chrom.length();
        if (length==4 && (chrom.charAt(3)=='X' || chrom.charAt(3)=='x')) return 0;
        else if (length==4 && (chrom.charAt(3)=='Y' || chrom.charAt(3)=='y')) return 1;
        else if (length==4 && (chrom.charAt(3)=='M' || chrom.charAt(3)=='m')) return 0;
        else if ( (length==4 && chrom.charAt(3)>='0' && chrom.charAt(3)<='9') ||
                  (length==5 && chrom.charAt(3)>='0' && chrom.charAt(3)<='9' && chrom.charAt(4)>='0' && chrom.charAt(4)<='9') 
                ) return Integer.parseInt(chrom.substring(3))%2==0?1:0;
        else return -1;
    }

}