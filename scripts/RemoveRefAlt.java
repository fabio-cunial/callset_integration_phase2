import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * Removes all the REF sequences that can be reconstructed from INFO, and stores
 * their length in the symbolic ALT to prevent overcollapse by `bcftools merge` 
 * downstream. Forces a QUAL value on every record.
 *
 * This program is typically used to compress a VCF that contains only ultralong
 * calls.
 *
 * Remark: POS values are not checked to be consistent with the reference.
 */
public class RemoveRefAlt {
    
    private static HashMap<String,Integer> fai;
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String FORCE_QUAL = args[1];
        final String INPUT_FAI = args[2];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        boolean refIsDna, altIsDna;
        char c;
        int i;
        int nRecords, nEdited, pos, svlen, nDiscarded;
        String str, ref, alt, info, svlenStr, endStr, chrom;
        BufferedReader br;
        String[] tokens;
        
        
        loadFai(INPUT_FAI);
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine(); nRecords=0; nEdited=0; nDiscarded=0;
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
            pos=Integer.parseInt(tokens[1]);  // 1-based, previous base in REF.
            if (pos<1 || pos>fai.get(chrom).intValue()) {
                // Discarded: wrong POS.
                nDiscarded++;
                str=br.readLine();
                continue;
            }
            ref=tokens[3]; alt=tokens[4]; info=tokens[7];
            c=ref.charAt(0);
            refIsDna = ref.length()>=1 && (c=='A' || c=='C' || c=='G' || c=='T' || c=='N' || c=='a' || c=='c' || c=='g' || c=='t' || c=='n');
            c=alt.charAt(0);
            altIsDna = alt.length()>=1 && (c=='A' || c=='C' || c=='G' || c=='T' || c=='N' || c=='a' || c=='c' || c=='g' || c=='t' || c=='n');
            
            // Forcing QUAL
            tokens[5]=FORCE_QUAL;
            
            // Computing SVLEN
            svlen=-1;
            if (refIsDna && altIsDna) {
                svlen=alt.length()-ref.length();
                if (svlen<0) svlen=-svlen;
            }
            else {
                svlenStr=getInfoField(info,"SVLEN");
                if (svlenStr!=null) {
                    svlen=Integer.parseInt(svlenStr);
                    if (svlen<0) svlen=-svlen;
                }
                else {
                    // 1-based, last base in REF.
                    endStr=getInfoField(info,"END");
                    if (endStr!=null) svlen=Integer.parseInt(endStr)-pos;
                    else {
                        // NOP: unknown SVLEN.
                        System.out.print(tokens[0]);
                        for (i=1; i<tokens.length; i++) { System.out.print('\t'); System.out.print(tokens[i]); }
                        System.out.println();
                        str=br.readLine();
                        continue;
                    }
                }
            }
            
            // Compressing the record
            if (info.indexOf("SVTYPE=DEL")>=0) { 
                if (info.indexOf("SVLEN=")<0) info+=";SVLEN="+svlen;
                alt="<DEL-"+svlen+">";
                ref="N";
                nEdited++;
            }
            else if (info.indexOf("SVTYPE=INV")>=0) { 
                if (info.indexOf("SVLEN=")<0) info+=";SVLEN="+svlen;
                alt="<INV-"+svlen+">";
                ref="N";
                nEdited++;
            }
            else if (info.indexOf("SVTYPE=DUP")>=0) { 
                if (info.indexOf("SVLEN=")<0) info+=";SVLEN="+svlen;
                alt="<DUP-"+svlen+">";
                ref="N";
                nEdited++;
            }
            else if (info.indexOf("SVTYPE=CNV")>=0) {
                if (info.indexOf("SVLEN=")<0) info+=";SVLEN="+svlen;
                alt="<CNV-"+svlen+">";
                ref="N";
                nEdited++;
            }
            tokens[3]=ref; tokens[4]=alt;
            
            // Outputting
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) { System.out.print('\t'); System.out.print(tokens[i]); }
            System.out.println();
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
        System.err.println("nEdited="+nEdited+" ("+((100.0*nEdited)/nRecords)+"% of nRecords)");
        System.err.println("nDiscarded="+nDiscarded+" ("+((100.0*nDiscarded)/nRecords)+"% of nRecords)");
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

}