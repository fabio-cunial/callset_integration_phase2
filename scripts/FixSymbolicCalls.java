import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * 
 */
public class FixSymbolicCalls {
    
    private static final char COMMENT = '#';
    private static HashMap<String,String> fasta;
    
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String OUTPUT_VCF_GZ = args[0];
        
        final int MIN_SV_LENGTH=Integer.parseInt(args[1]);
        final String OUTPUT_SV_VCF = args[2];
        final String OUTPUT_SNP_VCF = args[3];
        
        int row, length, nCalls, pos;
        String str, endStr, svlenStr;
        BufferedReader br;
        BufferedWriter bw1, bw2;
        String[] tokens;
        
        
        loadFasta(FASTA_PATH);
        
        bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(OUTPUT_VCF_GZ))));
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                bw.write(str); bw.newLine();
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            if (tokens[5].charAt(0)!='<') {
                bw.write(str); bw.newLine();
                str=br.readLine();
                continue;
            }
            if (tokens[5].substring(0,4).equalsIgnoreCase("<INS") || tokens[5].substring(0,4).equalsIgnoreCase("<CNV")) {
                // Discarded because of not enough information
                str=br.readLine();
                continue;
            }
            
            // Computing SVLEN
            pos=Integer.parseInt(tokens[1]);  // 1-based
            endStr=getInfoField(tokens[7],"END");  // 1-based, last base in REF.
            if (endStr!=null) length=Integer.parseInt(endStr)-pos;
            else {
                svlenStr=getInfoField(tokens[7],"SVLEN");
                if (svlenStr!=null) length=Integer.parseInt(svlenStr);
                else {
                    str=br.readLine();
                    continue;
                }
            }
            
            // Fixing the symbolic call
            if (tokens[5].substring(0,4).equalsIgnoreCase("<DEL")) {
                tokens[3]=fasta.get(tokens[0]).substring(pos-1,pos-1+length+1);
                tokens[4]=tokens[3].charAt(0)+"";
            }
            else if (tokens[5].substring(0,4).equalsIgnoreCase("<INV")) {
                tokens[3]=fasta.get(tokens[0]).substring(pos-1,pos-1+length+1);
                tokens[4]=tokens[3].charAt(0)+reverseComplement(tokens[3].substring(1));
            }
            else if (tokens[5].equalsIgnoreCase("<DUP>")) {
                tokens[4]=fasta.get(tokens[0]).substring(pos-1,pos-1+length+1);
                tokens[3]=tokens[4].charAt(0);
                tokens[7]=tokens[7].replace("SVTYPE=DUP;","SVTYPE=INS;");
            }
            if (tokens[5].charAt(0)=='<') tokens[5]=1;
            
            bw.write(tokens[0]);
            for (i=1; i<tokens.length; i++) { bw.write('\t'); bw.write(tokens[i]); }
            bw.newLine();
            str=br.readLine();
        }
        br.close(); bw.close();
        
        
        
        
        
    }
    
    
    
    
    
    
    
    
    private static final void loadFasta(String path) throws IOException {
        String currentChr, str;
        StringBuilder buffer;
        BufferedReader br;
        
        fasta = new HashMap<String,String>();
        buffer = new StringBuilder();
        br = new BufferedReader(new FileReader(path));
        str=br.readLine(); currentChr="";
        while (str!=null) {
            if (str.charAt(0)=='>') {
                if (currentChr.length()>0) fasta.put(currentChr,buffer.toString());
                currentChr=str.substring(1,str.indexOf(" "));
                buffer.delete(0,buffer.length());
            }
            else buffer.append(str);
        }
        if (currentChr.length()>0) fasta.put(currentChr,buffer.toString());
        br.close();
    }
    
    
	/**
	 * @return NULL if $field$ does not occur in $info$.
	 */
	private static final String getInfoField(String info, String field) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
        
        do { p=info.indexOf(field+"="); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return null;
		q=info.indexOf(";",p+FIELD_LENGTH);
		return info.substring(p+FIELD_LENGTH,q<0?info.length():q);
	}
    
    
    /**
     * @param str assumed to be uppercase.
     */
	private static final void replaceNonstandardChars(String str, StringBuilder out) {
        final int LENGTH = str.length();
        
        out.delete(0,out.length());
        for (int i=0; i<LENGTH; i++) {
            c=str.charAt(i);
            if (c=='A' || c=='C' || c=='G' || c=='T') out.append(c);
            else out.append('N');
        }
	}
    
    
    /**
     * @param str assumed to be uppercase.
     */
	private static final void reverseComplement(String str, StringBuilder out) {
        final int LENGTH = str.length();
        
        out.delete(0,out.length());
        for (int i=LENGTH-1; i>=0; i--) {
            switch (str.charAt(i)) {
                case 'A': out.append('T'); break;
                case 'C': out.append('G'); break;
                case 'G': out.append('C'); break;
                case 'T': out.append('A'); break;
                default: 'N'; break;
            }
        }
	}


}