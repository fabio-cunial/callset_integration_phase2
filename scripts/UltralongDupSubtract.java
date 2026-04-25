import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given an input VCF with only ultralong DUPs, the program subtracts SVLEN from
 * each POS.
 */
public class UltralongDupSubtract {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        int i;
        int pos, newPos, svlen;
        String str, info;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            pos=Integer.parseInt(tokens[1]);  // 1-based, exclusive.
            info=tokens[7];
            svlen=Integer.parseInt(getInfoField(info,"SVLEN"));
            newPos=pos-svlen;
            if (newPos<1) newPos=1;
            tokens[1]=newPos+"";

            // Outputting
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) System.out.print("\t"+tokens[i]);
            System.out.println();
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
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

}