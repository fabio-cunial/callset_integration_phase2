import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * Assumes that every symbolic record has SVLEN in INFO.
 */
public class AddSvlenToSymbolic {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String OUTPUT_VCF_GZ = args[1];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        int i;
        int nRecords, nEdited, svlen;
        String str, alt, info;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        
        bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(OUTPUT_VCF_GZ))));
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine(); nRecords=0; nEdited=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                bw.write(str); bw.newLine();
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");
            alt=tokens[4]; info=tokens[7];
            if (alt.charAt(0)!='<') {
                bw.write(str); bw.newLine();
                str=br.readLine();
                continue;
            }
            svlen=Integer.parseInt(getInfoField(info,"SVLEN"));
            if (info.indexOf("SVTYPE=DEL")>=0) { 
                alt="<DEL-"+svlen+">";
                nEdited++;
            }
            else if (info.indexOf("SVTYPE=INV")>=0) { 
                alt="<INV-"+svlen+">";
                nEdited++;
            }
            else if (info.indexOf("SVTYPE=DUP")>=0) { 
                alt="<DUP-"+svlen+">";
                nEdited++;
            }
            else if (info.indexOf("SVTYPE=CNV")>=0) {
                alt="<CNV-"+svlen+">";
                nEdited++;
            }
            tokens[4]=alt;
            
            // Outputting
            bw.write(tokens[0]);
            for (i=1; i<tokens.length; i++) { bw.write('\t'); bw.write(tokens[i]); }
            bw.newLine();
            
            // Next iteration
            str=br.readLine();
        }
        br.close(); bw.close();
        System.err.println("nRecords="+nRecords);
        System.err.println("nEdited="+nEdited+" ("+((100.0*nEdited)/nRecords)+"% of nRecords)");
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