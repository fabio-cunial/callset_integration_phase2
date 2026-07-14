import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a BND-only VCF, the program counts the number of records of each type.
 */
public class BndCount {

    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
    
        String str;
        BufferedReader br;
        int[][] counts = new int[4][3];
        String[] tokens;

        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            count(tokens);
            str=br.readLine();
        }
        br.close();
        System.out.println("TYPE | SAME CHR, alt>ref | SAME CHR, alt<ref | OTHER CHR");
        System.out.println("t[p[ | "+counts[0][0]+" | "+counts[0][1]+" | "+counts[0][2]);
        System.out.println("t]p] | "+counts[1][0]+" | "+counts[1][1]+" | "+counts[1][2]);
        System.out.println("]p]t | "+counts[2][0]+" | "+counts[2][1]+" | "+counts[2][2]);
        System.out.println("[p[t | "+counts[3][0]+" | "+counts[3][1]+" | "+counts[3][2]);
    }




    private static final void count(String[] tokens) {
        boolean sameChr, greaterAlt;
        int refPos, altPos;
        String str, alt;
        
        alt=tokens[4];
        refPos=Integer.parseInt(tokens[1]);
        altPos=getAltPos(alt);
        sameChr=tokens[0].equalsIgnoreCase(getAltChr(alt));
        greaterAlt=altPos>refPos;
        if (alt.charAt(alt.length()-1)=='[') counts[0][sameChr?(greaterAlt?0:1):2]++;
        else if (alt.charAt(alt.length()-1)==']') counts[1][sameChr?(greaterAlt?0:1):2]++;
        else if (alt.charAt(0)==']') counts[2][sameChr?(greaterAlt?0:1):2]++;
        else if (alt.charAt(0)=='[') counts[3][sameChr?(greaterAlt?0:1):2]++;
        else {
            System.err.println("ERROR: unrecognized ALT: "+alt);
            System.exit(1);
        }
    }


    private static String getAltChr(String alt) {
        int p, q, r;

        p=alt.indexOf('['); q=alt.indexOf(']'); 
        if (p>=0) {
            r=alt.indexOf(':',p+1);
            return alt.substring(p+1,r);
        }
        else if (q>=0) {
            r=alt.indexOf(':',q+1);
            return alt.substring(q+1,r);
        }
        else {
            System.err.println("ERROR: unrecognized ALT: "+alt);
            System.exit(1);
        }
        return null;
    }


    private static int getAltPos(String alt) {
        int p, q, r;

        p=alt.indexOf('['); q=alt.indexOf(']'); 
        if (p>=0) {
            r=alt.indexOf(':',p+1);
            p=alt.indexOf('[',r+1);
            return Integer.parseInt(alt.substring(r+1,p));
        }
        else if (q>=0) {
            r=alt.indexOf(':',q+1);
            q=alt.indexOf(']',r+1);
            return Integer.parseInt(alt.substring(r+1,q));
        }
        else {
            System.err.println("ERROR: unrecognized ALT: "+alt);
            System.exit(1);
        }
        return -1;
    }

}