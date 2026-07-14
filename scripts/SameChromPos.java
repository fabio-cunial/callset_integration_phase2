import java.util.zip.GZIPInputStream;
import java.io.*;


public class SameChromPos {
    /**
     * The program prints a CSV row with the following fields to STDOUT:
     * 
     * 0: total number of records;
     * 1: number of records with the same CHROM,POS as another record;
     * 
     * 2: number of distinct CHROM,POS pairs;
     * 3: number of CHROM,POS pairs with an INS and a DEL;
     * 4: number of CHROM,POS pairs with multiple INS;
     * 5: number of CHROM,POS pairs with multiple DEL.
     */
    public static int[] counts = new int[6];

    /**
     * @param args 0: assumed to be sorted by CHROM,POS.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];

        int i;
        int pos, currentPos, currentRecords, currentIns, currentDel;
        String str, svtype, currentChrom;
        BufferedReader br;
        String[] tokens;

        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        currentChrom=""; currentPos=-1; currentRecords=0; currentIns=0; currentDel=0;
        tokens=null; pos=-1;
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                str=br.readLine();
                continue;
            }
            counts[0]++;
            tokens=str.split("\t");
            pos=Integer.parseInt(tokens[1]);
            if (!tokens[0].equals(currentChrom) || pos!=currentPos) {
                // End of a block
                if (currentChrom.length()>0) {
                    counts[2]++;
                    if (currentRecords>1) {
                        counts[1]+=currentRecords;
                        if (currentIns>0 && currentDel>0) counts[3]++;
                        if (currentIns>1) counts[4]++;
                        if (currentDel>1) counts[5]++;
                    }
                }
                currentChrom=tokens[0]; currentPos=pos;
                currentRecords=0; currentIns=0; currentDel=0;
            }
            currentRecords++;
            svtype=getInfoField(tokens[7],"SVTYPE");
            if (svtype.equalsIgnoreCase("INS")) currentIns++;
            else if (svtype.equalsIgnoreCase("DEL")) currentDel++;
            str=br.readLine();
        }
        // End of last block
        if (currentChrom.length()>0) {
            counts[2]++;
            if (currentRecords>1) {
                counts[1]+=currentRecords;
                if (currentIns>0 && currentDel>0) counts[3]++;
                if (currentIns>1) counts[4]++;
                if (currentDel>1) counts[5]++;
            }
        }
        br.close();

        // Outputting
        System.out.print(counts[0]+"");
        for (i=1; i<counts.length; i++) System.out.print(","+counts[i]);
        System.out.println();
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