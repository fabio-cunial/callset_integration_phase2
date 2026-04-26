import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a VCF with only ultralong DUPs and a VCF with only ultralong INS, the 
 * program finds every DUP that has exactly one corresponding INS (i.e. there is
 * only one INS that is fully contained in the DUP and that has similar length).
 * 
 * Remark: length similarity is computed as in `truvari bench`.
 * 
 * Remark: this is quadratic just for simplicity, could be made much faster.
 */
public class UltralongBenchDupIns {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String DUP_VCF_GZ = args[0];
        final String INS_VCF_GZ = args[1];
        final int INS_VCF_NRECORDS = Integer.parseInt(args[2]);
        final double PCTSIZE = Double.parseDouble(args[3]);
;
        int i;
        int found;
        double pos, svlen;
        String str, chr, info;
        BufferedReader br;
        double[][] insPosLen;
        String[] insChr, tokens;

        // Loading the INS file in memory
        insPosLen = new double[INS_VCF_NRECORDS][2];
        insChr = new String[INS_VCF_NRECORDS];
        br = new BufferedReader( new InputStreamReader( (INS_VCF_GZ.length()>=7&&INS_VCF_GZ.substring(INS_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INS_VCF_GZ)) : new FileInputStream(INS_VCF_GZ) ) );
        str=br.readLine(); i=-1;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            i++;
            tokens=str.split("\t");
            chr=tokens[0];
            pos=Double.parseDouble(tokens[1]);  // 1-based, exclusive.
            info=tokens[7];
            svlen=Double.parseDouble(getInfoField(info,"SVLEN"));
            insChr[i]=chr; insPosLen[i][0]=pos; insPosLen[i][1]=svlen;
            str=br.readLine();
        }
        br.close();

        // Scanning the DUP file
        br = new BufferedReader( new InputStreamReader( (DUP_VCF_GZ.length()>=7&&DUP_VCF_GZ.substring(DUP_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(DUP_VCF_GZ)) : new FileInputStream(DUP_VCF_GZ) ) );
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            chr=tokens[0];
            pos=Double.parseDouble(tokens[1]);  // 1-based, exclusive.
            info=tokens[7];
            svlen=Double.parseDouble(getInfoField(info,"SVLEN"));
            found=0;
            for (i=0; i<insChr.length; i++) {
                if (insChr[i].equals(chr) && insPosLen[i][0]>=pos && insPosLen[i][0]<=pos+svlen && Math.min(svlen,insPosLen[i][1])/Math.max(svlen,insPosLen[i][1])>=PCTSIZE) found++;
            } 
            if (found==1) {
                System.out.print(tokens[0]);
                for (i=1; i<tokens.length; i++) System.out.print("\t"+tokens[i]);
                System.out.println();
            }
            
            // Next iteration
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