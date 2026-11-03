import java.util.Arrays;
import java.io.*;


/**
 * Prints histograms of `SUPP_*` and `SCORE` values of a VCF record over all of
 * its samples. Specifically, for every call, the program prints:
 * 
 * ID \t SUPP_HISTOGRAM \t SCORE_HISTOGRAM
 * 
 * where `SUPP_HISTOGRAM` is a comma-separated list of counts for each possible 
 * value of the bitmap where PBSV=1, SNIFFLES=2, PAV=4 (0 is the number of
 * samples with no support), and `SCORE_HISTOGRAM` is a comma-separated list of
 * counts for each equal-size bin in `[SCORE_MIN..SCORE_MAX]` (the first and 
 * last bin contain scores outside this range).
 */
public class GetDiscoveryAnnotations {
    
    /**
     * @args 0 format: `ID \t SAMPLE1 \t SAMPLE2 \t ...`, where every tab-
     * separated cell contains: `SUPP_PBSV,SUPP_SNIFFLES,SUPP_PAV,SCORE`. This
     * TSV was created e.g. with:
     * 
     * bcftools query --format '%ID[\t%SUPP_PBSV,%SUPP_SNIFFLES,%SUPP_PAV,%SCORE]\n'
     * 
     * over the truvari-collapsed inter-sample VCF before re-genotyping.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_TSV = args[0];
        final int N_SCORE_QUANTA = Integer.parseInt(args[1]);
        final double SCORE_MIN = Double.parseDouble(args[2]);  // >=0
        final double SCORE_MAX = Double.parseDouble(args[3]);  // >=0
        final String OUTPUT_TSV = args[4];
        
        final char SAMPLE_SEPARATOR = '\t';
        final int QUANTUM = 1000;  // Arbitrary
        
        boolean suppPbsv, suppSniffles, suppPav;
        int i;
        int length, nLines, start, end, suppID;
        long ms;
        double score;
        String str, id;
        BufferedReader br;
        BufferedWriter bw;
        int[] suppArray = new int[8];
        int[] scoreArray = new int[N_SCORE_QUANTA+2];
        
        br = new BufferedReader(new InputStreamReader(INPUT_TSV.equals("-")?System.in:new FileInputStream(INPUT_TSV)));
        bw = new BufferedWriter(new FileWriter(OUTPUT_TSV));
        str=br.readLine(); nLines=0;
        ms=System.currentTimeMillis();
        while (str!=null) {
            length=str.length();
            Arrays.fill(suppArray,0); Arrays.fill(scoreArray,0);
            start=str.indexOf(SAMPLE_SEPARATOR,0)+1;
            id=str.substring(0,start-1);
            do {
                end=str.indexOf(SAMPLE_SEPARATOR,start);
                if (end<0) end=length;
                suppPbsv=false; suppSniffles=false; suppPav=false; score=-1;
                if (str.charAt(start)!='.') suppPbsv=str.charAt(start)=='1';
                if (str.charAt(start+2)!='.') suppSniffles=str.charAt(start+2)=='1';
                if (str.charAt(start+4)!='.') suppPav=str.charAt(start+4)=='1';
                if (str.charAt(start+6)!='.') score=Double.parseDouble(str.substring(start+6,end));
                suppID=(suppPbsv?1:0)+(suppSniffles?2:0)+(suppPav?4:0);
                suppArray[suppID]++;
                if (score>=0 && score<SCORE_MIN) scoreArray[0]++;
                else if (score>=SCORE_MIN && score<=SCORE_MAX) scoreArray[1+Math.min((int)((score-SCORE_MIN)/(SCORE_MAX-SCORE_MIN)*N_SCORE_QUANTA),N_SCORE_QUANTA-1)]++;
                else if (score>SCORE_MAX) scoreArray[N_SCORE_QUANTA+1]++;
                start=end+1;
            } while (start<length);
            bw.write(id); bw.write('\t');
            bw.write(Integer.toString(suppArray[0]));
            for (i=1; i<8; i++) { bw.write(','); bw.write(Integer.toString(suppArray[i])); }
            bw.write('\t');
            bw.write(Integer.toString(scoreArray[0]));
            for (i=1; i<scoreArray.length; i++) { bw.write(','); bw.write(Integer.toString(scoreArray[i])); }
            bw.write("\n");
            nLines++;
            if (nLines%QUANTUM==0) System.err.println("Processed "+nLines+" lines in "+((double)(System.currentTimeMillis()-ms)/1000)+"s");
            str=br.readLine();
        }
        br.close(); bw.close();
    }
    
}