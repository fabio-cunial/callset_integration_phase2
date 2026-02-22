import java.io.*;


/**
 * A version of `TruvariDivide2.java` for the ultralong and BND VCFs.
 */
public class TruvariDivide2Ultralong {
    /**
     * @param args 5: 'ultralong' or 'bnd'.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV = args[0];
        final int BUFFER = Integer.parseInt(args[1]);  // Assumed >> 10
        final int MIN_RECORDS_PER_VCF = Integer.parseInt(args[2]);
        final String CHROM = args[3];
        final int EXPECTED_N_RECORDS_TOTAL = Integer.parseInt(args[4]);
        final int MODE = args[5].equalsIgnoreCase("ultralong")?0:1;
        
        final int SLACK = 5;  // Arbitrary
        
        int i, p, q;
        int chunkID, pos, svlen, nRecords, nRecordsTotal;
        int first, last, minFirst, maxLast;  // One-based, inclusive.
        String str;
        BufferedReader br;
        
        chunkID=0; nRecordsTotal=0;
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_TSV)));
        str=br.readLine(); nRecords=0; first=0; last=0; minFirst=0; maxLast=0;
        while (str!=null) {
            p=str.indexOf('\t');
            q=str.indexOf('\t',p+1);
            pos=Integer.parseInt(str.substring(0,p));            
            if (MODE==0) {
                if (str.charAt(q+1)=='<') {
                    svlen=Integer.parseInt(str.substring(q+6,str.length()-1));
                    first=pos+1;
                    last=pos+svlen-1;
                }
                else {  // INS
                    first=pos;
                    last=pos+1;
                }
            }
            else {  // BND
                first=pos;
                last=pos+1;
            }
            if (first>maxLast+BUFFER && nRecords>=MIN_RECORDS_PER_VCF) {
                System.out.println(CHROM+":"+(minFirst-SLACK)+"-"+(maxLast+SLACK)+" "+chunkID);
                nRecordsTotal+=nRecords;
                System.err.println("chunk="+chunkID+" nRecords="+nRecords);
                chunkID++;
                nRecords=1;
                minFirst=first;
                maxLast=last;
            }
            else {
                nRecords++;
                if (minFirst==0) minFirst=first;
                if (last>maxLast) maxLast=last;
            }
            str=br.readLine();
        }
        br.close();
        System.out.println(CHROM+":"+(minFirst-SLACK)+"-"+(maxLast+SLACK)+" "+chunkID);
        nRecordsTotal+=nRecords;
        System.err.println("chunk="+chunkID+" nRecords="+nRecords);
        if (nRecordsTotal!=EXPECTED_N_RECORDS_TOTAL) {
            System.err.println("ERROR: Expected "+EXPECTED_N_RECORDS_TOTAL+" records but created "+nRecords+" records.");
            System.exit(1);
        }
        else System.err.println("Created "+(chunkID+1)+" truvari chunks");
    }
    
}