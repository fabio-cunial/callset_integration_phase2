import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * 
 */
public class HapsVcf2Fasta {
    
    /**
     * @param 
     */
    public static void main(String[] args) throws IOException {
        final String HAPS_VCF_GZ = args[0];
        final String CHROMOSOME_FA = args[1];
        final int PADDING_BP = Integer.parseInt(args[2]);
        final String OUTPUT_DIR = args[3];
        
        final int QUANTUM = 10000;  // Arbitrary
        
        int p, q, r;
        int pos, currentPos, nRecords, chrLength;
        String str, chrom, currentChrom, ref, alt;
        BufferedReader br;
        BufferedWriter bw;
        StringBuilder chromosome;
        String[] tokens;
        
        System.err.print("Loading the reference... ");
        chromosome = new StringBuilder();
        br = new BufferedReader(new FileReader(CHROMOSOME_FA));
        str=br.readLine(); str=br.readLine();  // Skipping header
        while (str!=null) {
            chromosome.append(str);
            str=br.readLine();
        }
        br.close();
        chrLength=chromosome.length();
        System.err.println("done");
        
        System.err.print("Creating FASTAs... ");
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(HAPS_VCF_GZ))));
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/null.fa"));
        str=br.readLine(); currentChrom=""; currentPos=-1; nRecords=0;
        while (str!=null) {
            if (str.charAt(0)=='#') { str=br.readLine(); continue; }
            tokens=str.split("\t");
            chrom=tokens[0]; 
            pos=Integer.parseInt(tokens[1]);
            ref=tokens[3];
            p=(pos-1)-PADDING_BP;
            if (p<0) p=0;
            q=(pos-1)+ref.length();
            r=q+PADDING_BP;
            if (r>chrLength) r=chrLength;
            ref=chromosome.substring(p,pos-1)+ref+chromosome.substring(q,r);
            alt=tokens[4];
            alt=chromosome.substring(p,pos-1)+alt+chromosome.substring(q,r);
            if (!chrom.equals(currentChrom) || pos!=currentPos) {
                System.out.println(nRecords+"");
                bw.close();
                bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+chrom+"_"+pos+".fa"));
                bw.write(">ref\n");
                bw.write(ref); bw.newLine();
                bw.write(">alt_1\n");
                bw.write(alt); bw.newLine();
                currentChrom=chrom; currentPos=pos; nRecords=1;
            }
            else {
                nRecords++;
                bw.write(">alt_"+nRecords); bw.newLine();
                bw.write(alt); bw.newLine();
            }
            str=br.readLine();
        }
        br.close(); bw.close();
    }
    
}