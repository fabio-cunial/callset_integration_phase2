import java.io.*;


/**
 * 
 */
public class UltralongStats {
    
    /**
     * @args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final int MAX_SAMPLES = Integer.parseInt(args[1]);
        
        final int QUANTUM_SAMPLES = 10;
        final int MAX_LENGTH = 250000000;
        
        int length, nSamples;
        String str;
        BufferedReader br;
        String[] tokens;
        int[][] matrix_del, matrix_inv, matrix_dup, matrix_ins;
        
        // Building matrices
        matrix_del = new int[MAX_SAMPLES+1][(int)Math.log10(MAX_LENGTH)+1];
        matrix_inv = new int[MAX_SAMPLES+1][(int)Math.log10(MAX_LENGTH)+1];
        matrix_dup = new int[MAX_SAMPLES+1][(int)Math.log10(MAX_LENGTH)+1];
        matrix_ins = new int[MAX_SAMPLES+1][(int)Math.log10(MAX_LENGTH)+1];
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            length=(int)Math.log10(Integer.parseInt(tokens[1]));
            if (length>=matrix_del[0].length) length=matrix_del[0].length;
            nSamples=Integer.parseInt(tokens[2]);
            if (tokens[0].equalsIgnoreCase("DEL")) matrix_del[nSamples][length]++;
            else if (tokens[0].equalsIgnoreCase("INV")) matrix_inv[nSamples][length]++;
            else if (tokens[0].equalsIgnoreCase("DUP")) matrix_dup[nSamples][length]++;
            else if (tokens[0].equalsIgnoreCase("INS")) matrix_ins[nSamples][length]++;
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        printMatrix(matrix_del,"matrix_del.csv");
        printMatrix(matrix_inv,"matrix_inv.csv");
        printMatrix(matrix_dup,"matrix_dup.csv");
        printMatrix(matrix_ins,"matrix_ins.csv");
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