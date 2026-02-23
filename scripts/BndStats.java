import java.io.*;


/**
 * 
 */
public class BndStats {
    
    /**
     * @args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        
        int i;
        int from, to;
        String str;
        BufferedReader br;
        String[] tokens;
        int[][] matrix;
        
        // Building matrix
        matrix = new int[24][25];
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            if (tokens[0].charAt(3)=='X') from=23;
            else if (tokens[0].charAt(3)=='Y') from=24;
            else from=Integer.parseInt(tokens[0].substring(3));
            to=0;
            for (i=22; i>=1; i--) {
                if (tokens[1].indexOf("chr"+i)>=0) { to=i; break; }
            }
            if (to==0) {
                if (tokens[1].indexOf("chrX")>=0) to=23;
                else if (tokens[1].indexOf("chrY")>=0) to=24;
                else to=25;  // Non-standard chromosome
            }
            matrix[from-1][to-1]++;
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        printMatrix(matrix,"matrix_bnd.csv");
    }
    
    
    private static final void printMatrix(int[][] matrix, String path) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(path));
        for (int i=0; i<matrix.length; i++) {
            for (int j=0; j<matrix[i].length; j++) bw.write(matrix[i][j]+",");
            bw.newLine();
        }
        bw.close();
    }

}