import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * Makes sure that SVTYPE and SVLEN are present and computed consistently.
 * This is essentially a reimplementation of `truvari anno svinfo`, which is 
 * necessary since:
 *
 * 1. That program does not preserve the input SVTYPE, so e.g. a raw call with 
 * SVTYPE=INV and explicit REF and ALT, but such that ALT is shorter than REF by
 * a few bases, is reassigned an SVTYPE=DEL. Instead, we want to preserve the
 * original SVTYPE to guide merging and for downstream interpretation.
 *
 * 2. In case of a substitution with explicit REF and ALT, that program sets
 * SVLEN=| |ALT|-|REF| |. If downstream we use SVLEN to e.g. filter for large
 * SVs, this is incorrect. Instead, we want SVLEN=|REF| for substitutions.
 */
public class AddSvtypeSvlen {
    
    /**
     * @param args 0: a set of raw calls from an SV caller.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        boolean edited;
        int i, p;
        int nRecords, refType, altType, svlen, pos, value;
        String str, ref, alt, info, svtype;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
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
            pos=Integer.parseInt(tokens[1]); ref=tokens[3]; alt=tokens[4]; info=tokens[7];
            refType=getAltType(ref); altType=getAltType(alt);
            
            // Deciding SVTYPE
            edited=false;
            svtype=getInfoField(info,"SVTYPE");
            if (svtype==null) {
                if (altType==1) svtype=alt.substring(1,4);
                else if (altType==2) svtype="BND";
                else if (altType==0 && refType==0) {
                    if (ref.length()==1 && alt.length()>1) svtype="INS";
                    else if (ref.length()>1 && alt.length()==1) svtype="DEL";
                    else if (ref.length()>1 && alt.length()>1) svtype="SUB";
                    else svtype="UNK";
                }
                else svtype="UNK";
                edited=true;
            }
            else if (svtype.length()>3) { svtype=svtype.substring(0,3); edited=true; }
            if (edited) info=addOrReplaceInfoField(info,"SVTYPE",svtype);
            
            // Computing SVLEN
            svlen=0;
            if (svtype.equalsIgnoreCase("BND")) {
                // NOP
            }
            else {
                if (refType==0 && altType==0) {
                    if (svtype.equalsIgnoreCase("INS")) svlen=alt.length()-1;
                    else svlen=ref.length()-1;
                }
                else if (altType==1) {
                    value=getSymbolicSvlen(pos,info);
                    if (value!=-1) svlen=value;
                }
            }
            if (svlen>0 || isInfoFieldPresent(info,"SVLEN")) info=addOrReplaceInfoField(info,"SVLEN",svlen+"");
            
            tokens[7]=info;
            
            // Outputting
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) { System.out.print('\t'); System.out.print(tokens[i]); }
            System.out.println();
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
    }
    
    
    /**
     * Since the input is assumed to come from a raw call, the procedure allows
     * the entire IUPAC alphabet, see e.g.:
     *
     * https://en.wikipedia.org/wiki/Nucleic_acid_sequence
     *
     * @return
     * -1: unknown;
     *  0: fully DNA;
     *  1: symbolic;
     *  2: BND.
     */
    private static final int getAltType(String alt) {
        char a, z, c;
        boolean found;
        int i;
        final String ALPHABET = "acgtuwsmkrybdhvnz";
        final int ALT_LENGTH = alt.length();
        
        a=alt.charAt(0); z=alt.charAt(ALT_LENGTH-1);
        if (a=='<' && z=='>') return 1;
        else if (a=='.') {
            found=true;
            for (i=1; i<ALT_LENGTH; i++) {
                if (ALPHABET.indexOf(Character.toLowerCase(alt.charAt(i)))<0) { found=false; break; }
            }
            if (found) return 2;
        }
        else if (z=='.') {
            found=true;
            for (i=0; i<ALT_LENGTH-1; i++) {
                if (ALPHABET.indexOf(Character.toLowerCase(alt.charAt(i)))<0) { found=false; break; }
            }
            if (found) return 2;
        }
        else if ( (a==']' && alt.indexOf("]",1)>=0) ||
                  (a=='[' && alt.indexOf("[",1)>=0) ||
                  (z=='[' && alt.indexOf("[")<ALT_LENGTH-1) ||
                  (z==']' && alt.indexOf("]")<ALT_LENGTH-1)
                ) return 2;
        else {
            found=true;
            for (i=0; i<ALT_LENGTH; i++) {
                if (ALPHABET.indexOf(Character.toLowerCase(alt.charAt(i)))<0) { found=false; break; }
            }
            if (found) return 0;
        }
        return -1;
    }
    
    
    /**
     * @return -1 if SVLEN cannot be determined from `info`. Otherwise, a 
     * non-negative integer.
     */
    private static int getSymbolicSvlen(int pos, String info) {
        int end, out;
        String value;
        
        out=-1;
        value=getInfoField(info,"SVLEN");
        if (value!=null) {
            out=Integer.parseInt(value);
            if (out<0) out=-out;
        }
        else {
            value=getInfoField(info,"END");
            if (value!=null) {
                end=Integer.parseInt(value);
                if (end>pos) out=end-pos;
            }
        }
        return out;
    }
    
    
	private static final boolean isInfoFieldPresent(String info, String field) {
		final int FIELD_LENGTH = field.length()+1;
        int p;
        
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
        return p>=0;
	}
    
    
	/**
	 * @return NULL if `field` does not occur in `info`.
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
    
    
	private static final String addOrReplaceInfoField(String info, String field, String newValue) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
    
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return info+";"+field+"="+newValue;
		q=info.indexOf(";",p+FIELD_LENGTH);
		return info.substring(0,p+FIELD_LENGTH)+newValue+info.substring(q);
	}

}