import java.io.*;

/**
 * 
 */
public class CountDeNovoSimple {
    
    /**
     * @param args a CSV with format `GT_child,GT_father,GT_mother`, built with
     * e.g. `bcftools query --format '[%GT,]\n' trio.vcf.gz > trio.csv`.
     */
    public static void main(String[] args) throws IOException {
        final String TRIO_MATRIX_CSV = args[0];
        
        final int QUANTUM = 10000;  // Arbitrary
        
        int numerator, denominator;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(TRIO_MATRIX_CSV));
        numerator=0; denominator=0;
        str=br.readLine();
        while (str!=null) { 
            tokens=str.split(",");
            if (tokens[0].indexOf("1")<0) {
                str=br.readLine();
                continue;
            }
            if (tokens[0].indexOf(".")>=0 || tokens[1].indexOf(".")>=0 || tokens[2].indexOf(".")>=0) {
                str=br.readLine();
                continue;
            }
            denominator++;
            if (tokens[1].indexOf("1")<0 && tokens[2].indexOf("1")<0) numerator++;
            str=br.readLine();
        }
        br.close();
        System.out.println(numerator+","+denominator);
    }
    
}