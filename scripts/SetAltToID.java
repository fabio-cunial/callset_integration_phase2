import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * 
 */
public class SetAltToID {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        int i;
        int nRecords;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader( new InputStreamReader( INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz") ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");
            tokens[4]=tokens[2];
            tokens[7]=addOrReplaceInfoField(tokens[7],"SVTYPE","DEL");
            System.out.println(String.join("\t",tokens));
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
    }


    private static final String addOrReplaceInfoField(String info, String field, String newValue) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
        
        if (info.equals(".")) return field+"="+newValue;
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return info+";"+field+"="+newValue;
		q=info.indexOf(";",p+FIELD_LENGTH);
        return info.substring(0,p+FIELD_LENGTH)+newValue+(q>=0?info.substring(q):"");
	}

}