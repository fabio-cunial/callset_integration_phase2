import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * Collapses all the SAMPLE columns of a VCF into the first column. For each
 * SAMPLE field, the program picks the first non-missing value in the order of 
 * the SAMPLE columns, as described in:
 *
 * https://github.com/acenglish/truvari/wiki/collapse#--intra
 *
 * The program also adds a truvari-collapse-like SUPP field to SAMPLE.
 */
public class CollapseSamples {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String OUTPUT_VCF_GZ = args[1];
        
        final int QUANTUM = 5000;  // Arbitrary
        
        int i, j;
        int nRecords, nFields, gtIndex, gtCount, count, supp;
        String str, alt, info;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens, tokensPrime, outputFields;
        
        outputFields = new String[100];  // Arbitrary
        bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(OUTPUT_VCF_GZ))));
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine(); nRecords=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.startsWith("#CHROM")) {
                    bw.write("##FORMAT=<ID=SUPP,Number=1,Type=Integer,Description=\"Truvari-collapse-like support flag\">\n");
                    tokens=str.split("\t");
                    bw.write(tokens[0]);
                    for (i=1; i<=9; i++) { bw.write('\t'); bw.write(tokens[i]); }
                    bw.newLine();
                }
                else { bw.write(str); bw.newLine(); }
                str=br.readLine();
                continue;
            }
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records");
            tokens=str.split("\t");
            
            // Initializing GT fields
            tokensPrime=tokens[8].split(":");
            nFields=tokensPrime.length; gtIndex=-1;
            for (i=0; i<nFields; i++) {
                if (tokensPrime[i].equalsIgnoreCase("GT")) gtIndex=i;
                outputFields[i]=".";
            }
            gtCount=0; supp=0;
            
            // Collapsing GT fields and computing SUPP
            for (j=9; j<tokens.length; j++) {
                tokensPrime=tokens[j].split(":");
                count=(tokensPrime[gtIndex].charAt(0)!='.'?1:0)+(tokensPrime[gtIndex].charAt(2)!='.'?1:0);
                if (count>gtCount) outputFields[gtIndex]=tokensPrime[gtIndex];
                if (tokensPrime[gtIndex].charAt(0)=='1' || tokensPrime[gtIndex].charAt(2)=='1') supp|=1<<(j-9);
                for (i=0; i<nFields; i++) {
                    if (i!=gtIndex && outputFields[i].equals(".")) outputFields[i]=tokensPrime[i];
                }
            }
            
            // Outputting
            bw.write(tokens[0]);
            for (i=1; i<=8; i++) { bw.write('\t'); bw.write(tokens[i]); }
            bw.write(":SUPP\t"); bw.write(outputFields[0]);
            for (i=1; i<nFields; i++) { bw.write(':'); bw.write(outputFields[i]); }
            bw.write(':'); bw.write(""+supp);
            bw.newLine();
            
            // Next iteration
            str=br.readLine();
        }
        br.close(); bw.close();
        System.err.println("nRecords="+nRecords);
    }

}