import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * 
 */
public class FixUltralongRecords {
    
    private static HashMap<String,Integer> fai;
    
    
    /**
     * @param args
     * 0: assumed not to contain ultra-long calls or BND calls.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String INPUT_FAI = args[1];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        int i;
        int chromLength, pos, svlen, nRecords, nDiscarded_pos, nDiscarded_svlen;
        String chrom, ref, alt, info, str, svlenStr, endStr;
        BufferedReader br;
        String[] tokens;
        

        loadFai(INPUT_FAI);
        br = new BufferedReader( new InputStreamReader( INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz") ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0; nDiscarded_pos=0; nDiscarded_svlen=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");
            chrom=tokens[0];
            chromLength=fai.get(chrom).intValue();
            pos=Integer.parseInt(tokens[1]);  // 1-based, previous base in REF.
            if (pos<1 || pos>chromLength) {
                nDiscarded_pos++;
                str=br.readLine();
                continue;
            }
            alt=tokens[4];
            if (alt.charAt(0)!='<') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            info=tokens[7];
            svlen=Integer.parseInt(getInfoField(info,"SVLEN"));
            if (pos+svlen>chromLength) {
                nDiscarded_svlen++;
                System.err.println("Call ends beyond "+chrom+": "+(pos+svlen)+">"+chromLength+" (excess: "+(pos+svlen-chromLength)+")");
                str=br.readLine();
                continue;
            }
            
            // Outputting
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) { System.out.print('\t'); System.out.print(tokens[i]); }
            System.out.println();
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
        System.err.println("nDiscarded_pos="+nDiscarded_pos+" ("+((100.0*nDiscarded_pos)/nRecords)+"% of nRecords)");
        System.err.println("nDiscarded_svlen="+nDiscarded_svlen+" ("+((100.0*nDiscarded_svlen)/nRecords)+"% of nRecords)");
    }
    
    
    private static final void loadFai(String path) throws IOException {
        String str;
        BufferedReader br;
        String[] tokens;
        
        fai = new HashMap<String,Integer>();
        br = new BufferedReader(new FileReader(path));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            fai.put(tokens[0],Integer.valueOf(tokens[1]));
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