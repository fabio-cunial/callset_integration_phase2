import java.io.*;


/**
 * Resets SVLEN to $|REF|-1$ for every record where REF and ALT are DNA 
 * sequences of length >1. This is important for INV and replacement records,
 * since `truvari anno svinfo` sets SVLEN to $| |ALT|-|REF| |$ in this case.
 */
public class FixTruvariAnnoSvlen {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        boolean refIsDna, altIsDna;
        char c;
        int i;
        int nRecords, nEdited;
        String str, ref, alt, info;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_VCF)));
        str=br.readLine(); nRecords=0; nEdited=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");
            ref=tokens[3]; alt=tokens[4]; info=tokens[7];
            
            // Fixing SVLEN
            if (isDna(ref) && isDna(alt) && ref.length()>1 && alt.length()>1) {
                info=replaceInfoField(info,"SVLEN",(ref.length()-1)+"");
                nEdited++;
            }
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
        System.err.println("nEdited="+nEdited+" ("+((100.0*nEdited)/nRecords)+"% of nRecords)");
    }
    
    
    private static final boolean isDna(String str) {
        char c;
        int i;
        final int STR_LENGTH = str.length();
        
        for (i=0; i<STR_LENGTH; i++) {
            c=str.charAt(i);
            if (c!='A' && c!='C' && c!='G' && c!='T' && c!='N' && c!='a' && c!='c' && c!='g' && c!='t' && c!='n') return false;
        }
        return true;
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