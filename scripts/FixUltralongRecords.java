import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * - Discards records with first or last position outside the chromosome.
 * - Ensures that the correct value of END is present in every record.
 * - Removes SVLEN from symbolic ALTs.
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
        
        boolean endFound;
        int i;
        int chromLength, pos, svlen, end, nRecords, nDiscarded_pos, nDiscarded_svlen, nNoEnd, nWrongEnd;
        String chrom, ref, alt, info, str, svlenStr, endStr;
        BufferedReader br;
        String[] tokens;
        

        loadFai(INPUT_FAI);
        br = new BufferedReader( new InputStreamReader( INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz") ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0; endFound=false; nDiscarded_pos=0; nDiscarded_svlen=0; nNoEnd=0; nWrongEnd=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.indexOf("INFO=<ID=END,")>=0) endFound=true;
                if (!endFound && str.substring(0,6).equals("#CHROM")) System.out.println("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">");
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
            alt=tokens[4]; info=tokens[7];
            svlen=Integer.parseInt(getInfoField(info,"SVLEN"));
            if (alt.charAt(0)=='<') {
                if (pos+svlen>chromLength) {
                    nDiscarded_svlen++;
                    System.err.println("Call ends after "+chrom+": "+(pos+svlen)+">"+chromLength+" (excess: "+(pos+svlen-chromLength)+")");
                    str=br.readLine();
                    continue;
                }
                alt=alt.substring(0,4)+'>';
                end=pos+svlen;
            }
            else end=pos;
            endStr=getInfoField(info,"END");
            if (endStr==null) {
                nNoEnd++;
                info+=";END="+end;
            }
            else if (Integer.parseInt(endStr)!=end) {
                System.err.println("Call "+tokens[2]+" has wrong END: "+endStr+" != "+end);
                nWrongEnd++;
                info=replaceInfoField(info,"END",end+"");
            }
            tokens[4]=alt; tokens[7]=info;
            
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
        System.err.println("nNoEnd="+nNoEnd+" ("+((100.0*nNoEnd)/nRecords)+"% of nRecords)");
        System.err.println("nWrongEnd="+nWrongEnd+" ("+((100.0*nWrongEnd)/nRecords)+"% of nRecords)");
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
    
    
    /**
     * @return `ìnfo` if `field` does not occur in `info`.
     */
	private static final String replaceInfoField(String info, String field, String newValue) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
    
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return info;
		q=info.indexOf(";",p+FIELD_LENGTH);
		return info.substring(0,p+FIELD_LENGTH)+newValue+info.substring(q);
	}

}