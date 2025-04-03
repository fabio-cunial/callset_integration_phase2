import java.io.*;


/**
 * 
 */
public class FillTablePhase2 {
    
    /**
     * @args 
     * 0: TSV with format: 
     *    0: id
     *    1: 02_aligned_bai
     *    2: 02_aligned_bam
     *    3: 03_pav_bed
     *    4: 03_pav_tbi
     *    5: 03_pav_vcf
     *    6: 03_pbsv_tbi
     *    7: 03_pbsv_vcf
     *    8: 03_sniffles_tbi
     *    9: 03_sniffles_vcf
     * 1: gs://fc-secure-8f7d6a20-04ce-40d7-8c88-aececeac3e09/CCS/terra-f6a059a8/outputs/T2T
     */
    public static void main(String[] args) throws Exception {
        final String INPUT_FILE = args[0];
        final String ROOT_DIR = args[1];
        
        int i;
        String str, sampleId;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t",-1);
            sampleId=tokens[0];
            if (tokens[1].length()==0) tokens[1]=ROOT_DIR+"/alignments/"+sampleId+".bam.bai";
            if (tokens[2].length()==0) tokens[2]=ROOT_DIR+"/alignments/"+sampleId+".bam";
            if (tokens[3].length()==0) tokens[3]=ROOT_DIR+"/PAV/"+sampleId+"/"+sampleId+".pav.final_bed.tgz";
            if (tokens[4].length()==0) tokens[4]=ROOT_DIR+"/PAV/"+sampleId+"/pav_"+sampleId+".vcf.gz.tbi";
            if (tokens[5].length()==0) tokens[5]=ROOT_DIR+"/PAV/"+sampleId+"/pav_"+sampleId+".vcf.gz";
            if (tokens[6].length()==0) tokens[6]=ROOT_DIR+"/variants/sv/"+sampleId+".pbsv.vcf.gz.tbi";
            if (tokens[7].length()==0) tokens[7]=ROOT_DIR+"/variants/sv/"+sampleId+".pbsv.vcf.gz";
            if (tokens[8].length()==0) tokens[8]=ROOT_DIR+"/variants/sv/"+sampleId+".sniffles.vcf.gz.tbi";
            if (tokens[9].length()==0) tokens[9]=ROOT_DIR+"/variants/sv/"+sampleId+".sniffles.vcf.gz";
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) System.out.print("\t"+tokens[i]);
            System.out.println();
            str=br.readLine();
        }
        br.close();
    }

}