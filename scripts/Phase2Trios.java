import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class Phase2Trios {
    
    /**
     * @param args
     * 0: CSV with columns: $pedigree_id,child_id,father_id,mother_id,child_sex$;
     * 1: TSV from Terra, with columns: $id,coverage,sex$. Rows must be sorted
     *    by $id$.
     */
    public static void main(String[] args) throws IOException {
        final String PEDIGREES_CSV = args[0];
        final String COVERAGES_TSV = args[1];  // ID, coverage, sex.
        
        int i;
        int child, father, mother, sexChild;
        double coverageChild, coverageFather, coverageMother;
        String str;
        BufferedReader br;
        int[] ids, sex;
        double[] coverages;
        String[] tokens;
        
        // Loading coverages and sex from the Terra table
        br = new BufferedReader(new FileReader(COVERAGES_TSV));
        str=br.readLine(); i=0;
        while (str!=null) { i++; str=br.readLine(); }
        br.close();
        ids = new int[i]; coverages = new double[i]; sex = new int[i];
        br = new BufferedReader(new FileReader(COVERAGES_TSV));
        str=br.readLine(); str=br.readLine();  // Skipping header
        i=0;
        while (str!=null) {
            tokens=str.split("\t");
            ids[i]=Integer.parseInt(tokens[0]);
            coverages[i]=Double.parseDouble(tokens[1]);
            sex[i]=tokens[2].equalsIgnoreCase("M")?1:2;
            i++;
            str=br.readLine();
        }
        br.close();
        
        // Parsing pedigrees
        br = new BufferedReader(new FileReader(PEDIGREES_CSV));
        str=br.readLine(); str=br.readLine();  // Skipping header
        while (str!=null) {            
            tokens=str.split(",");
            child=Integer.parseInt(tokens[1]);
            father=Integer.parseInt(tokens[2]);
            mother=Integer.parseInt(tokens[3]);
            if (father==0 || mother==0) {
                str=br.readLine();
                continue;
            }
            i=Arrays.binarySearch(ids,child);
            if (i<0) {
                System.err.println("Not found: child "+child);
                str=br.readLine();
                continue;
            }
            coverageChild=coverages[i]; sexChild=sex[i];
            i=Arrays.binarySearch(ids,father);
            if (i<0) {
                System.err.println("Not found: father "+father);
                str=br.readLine();
                continue;
            }
            coverageFather=coverages[i];
            i=Arrays.binarySearch(ids,mother);
            if (i<0) {
                System.err.println("Not found: mother "+mother);
                str=br.readLine();
                continue;
            }
            coverageMother=coverages[i];
            System.out.println(tokens[0]+","+child+","+father+","+mother+","+sexChild+","+coverageChild+","+coverageFather+","+coverageMother);
            str=br.readLine();
        }
        br.close();
    }
    
}