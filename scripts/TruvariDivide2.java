import java.io.*;


/**
 * A faster, smaller-memory implementation of `truvari divide`, which assumes
 * the input VCF contains only one chromosome, is sorted, and has sequence-
 * resolved REF/ALT.
 * This program takes in input a TSV with fields:
 *
 * POS, REF, ALT
 *
 * and it prints to STDOUT regions for `bcftools view`.
 */
public class TruvariDivide2 {
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV = args[0];
        final int BUFFER = Integer.parseInt(args[1]);  // Assumed >> 10
        final int MIN_RECORDS_PER_VCF = Integer.parseInt(args[2]);
        final String CHROM = args[3];
        
        final int SLACK = 5;  // Arbitrary
        
        int i, p, q;
        int chunkID, pos, refLength, altLength, nRecords;
        int first, last, minFirst, maxLast;  // One-based, inclusive.
        String str;
        BufferedReader br;
        
        chunkID=0;
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_TSV)));
        str=br.readLine(); nRecords=0; first=0; last=0; minFirst=0; maxLast=0;
        while (str!=null) {
            p=str.indexOf('\t');
            q=str.indexOf('\t',p+1);
            pos=Integer.parseInt(str.substring(0,p));
            refLength=q-p-1; altLength=str.length()-q-1;  // Length of fields
            if (refLength>1 && altLength==1) {
                // DEL
                first=pos+1;
                last=pos+refLength-1;
            }
            else if (refLength==1 && altLength>1) {
                // INS
                first=pos;
                last=pos+1;
            }
            else if (refLength>1 && altLength>1) {
                // Substitution
                first=pos+1;
                last=pos+refLength-1;
            }
            else {
                System.err.println("ERROR: Record type not recognized:");
                System.err.println(str);
                System.exit(1);
            }
            if (first>maxLast+BUFFER && nRecords>=MIN_RECORDS_PER_VCF) {
                System.out.println(CHROM+":"+(minFirst-SLACK)+"-"+(maxLast+SLACK)+" "+chunkID);
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
        System.err.println("chunk="+chunkID+" nRecords="+nRecords);
        System.err.println("Created "+(chunkID+1)+" truvari chunks");
    }
    
}