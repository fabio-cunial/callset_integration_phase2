import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given an input VCF with only ultralong DUPs, the program adds to each POS its
 * SVLEN.
 */
public class UltralongDupAdd {
    
    private static HashMap<String,Integer> fai;
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String INPUT_FAI = args[1];
        
        int i;
        int chromLength, pos, newPos, svlen;
        String str, chrom, info;
        BufferedReader br;
        String[] tokens;
        
        loadFai(INPUT_FAI);
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            chrom=tokens[0];
            chromLength=fai.get(chrom).intValue();
            pos=Integer.parseInt(tokens[1]);  // 1-based, exclusive.
            info=tokens[7];
            svlen=Integer.parseInt(getInfoField(info,"SVLEN"));
            newPos=pos+svlen;
            if (newPos>chromLength) newPos=chromLength;
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