import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * Tries to transform a symbolic record into a non-symbolic record, and discards 
 * a symbolic record if the transformation is not possible.
 */
public class FixSymbolicRecords {
    
    private static HashMap<String,StringBuilder> fasta;
    
    
    /**
     * @param args
     * 0: assumed not to contain ultra-long calls or BND calls.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String INPUT_FASTA = args[1];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        int i, p, q;
        int pos, svlen, nRecords, nSymbolicRecords, nDiscarded;
        long time;
        String chrom, ref, alt, info, str, svlenStr, endStr;
        StringBuilder buffer;
        BufferedReader br;
        String[] tokens;
        
        System.err.print("Loading reference...");
        time=System.currentTimeMillis();
        loadFasta(INPUT_FASTA);
        System.err.println(" done in "+((System.currentTimeMillis()-time)/1000)+"s");
        
        buffer = new StringBuilder();
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine(); nRecords=0; nSymbolicRecords=0; nDiscarded=0;
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
            if (pos<1 || pos>fasta.get(chrom).length()) {
                // Discarded: wrong POS.
                str=br.readLine();
                continue;
            }
            alt=tokens[4];
            if (alt.charAt(0)!='<') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            nSymbolicRecords++;
            if (alt.substring(0,4).equalsIgnoreCase("<INS") || alt.substring(0,4).equalsIgnoreCase("<CNV")) {
                // Discarded: not enough information.
                nDiscarded++;
                str=br.readLine();
                continue;
            }
            
            // Computing SVLEN
            info=tokens[7];
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
                    // Discarded: not enough information.
                    nDiscarded++;
                    str=br.readLine();
                    continue;
                }
            }
            
            // Fixing the symbolic record
            ref=tokens[3];
            if (alt.substring(0,4).equalsIgnoreCase("<DEL")) {
                ref=fasta.get(chrom).substring(pos-1,pos-1+svlen+1);
                alt=ref.charAt(0)+"";
            }
            else if (alt.substring(0,4).equalsIgnoreCase("<INV")) {
                ref=fasta.get(chrom).substring(pos-1,pos-1+svlen+1);
                reverseComplement(ref.substring(1),buffer);
                alt=ref.charAt(0)+buffer.toString();
            }
            else if (alt.substring(0,4).equalsIgnoreCase("<DUP")) {
                alt=fasta.get(chrom).substring(pos-1,pos-1+svlen+1);
                ref=alt.charAt(0)+"";
                info=info.replace("SVTYPE=DUP","SVTYPE=INS");
                p=info.indexOf("END=");
                if (p>=0) {
                    q=info.indexOf(";",p+1);
                    info=info.substring(0,p)+"END="+pos+(q>=0?info.substring(q):"");
                }
            }
            else {
                // Discarded: unknown type.
                nDiscarded++;
                str=br.readLine();
                continue;
            }
            tokens[3]=ref; tokens[4]=alt; tokens[7]=info;
            
            // Outputting
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) { System.out.print('\t'); System.out.print(tokens[i]); }
            System.out.println();
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        System.err.println("nRecords="+nRecords);
        System.err.println("nSymbolicRecords="+nSymbolicRecords+" ("+((100.0*nSymbolicRecords)/nRecords)+"% of nRecords)");
        System.err.println("nDiscarded="+nDiscarded+" ("+((100.0*nDiscarded)/nSymbolicRecords)+"% of nSymbolicRecords)");
    }
    
    
    private static final void loadFasta(String path) throws IOException {
        String currentChr, str;
        StringBuilder buffer;
        BufferedReader br;
        
        fasta = new HashMap<String,StringBuilder>();
        br = new BufferedReader(new FileReader(path));
        str=br.readLine(); currentChr=""; buffer=null;
        while (str!=null) {
            if (str.charAt(0)=='>') {
                if (currentChr.length()>0) fasta.put(currentChr,buffer);
                currentChr=str.substring(1,str.indexOf(" "));
                buffer = new StringBuilder();
            }
            else buffer.append(str);
            str=br.readLine();
        }
        if (currentChr.length()>0) fasta.put(currentChr,buffer);
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
    
    
	private static final void reverseComplement(String str, StringBuilder out) {
        final int LENGTH = str.length();
        
        out.delete(0,out.length());
        for (int i=LENGTH-1; i>=0; i--) {
            switch (str.charAt(i)) {
                case 'A': out.append('T'); break;
                case 'C': out.append('G'); break;
                case 'G': out.append('C'); break;
                case 'T': out.append('A'); break;
                case 'a': out.append('t'); break;
                case 'c': out.append('g'); break;
                case 'g': out.append('c'); break;
                case 't': out.append('a'); break;
                default: out.append('N'); break;
            }
        }
	}    

}