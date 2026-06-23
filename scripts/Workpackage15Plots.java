import java.io.*;


/**
 * 
 */
public class Workpackage15Plots {
    
    /**
     * @args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final int MAX_FREQUENCY = Integer.parseInt(args[1]);
        final int LENGTH_QUANTUM = Integer.parseInt(args[2]);
        
        final int MAX_LENGTH = 2000000;
        final int N_LENGTH_BINS = MAX_LENGTH / LENGTH_QUANTUM + 1;
        
        int length, frequency;
        String str;
        BufferedReader br;
        String[] tokens;
        int[][] matrix_del, matrix_inv, matrix_dup, matrix_ins;
        
        // Building matrices
        matrix_del = new int[MAX_FREQUENCY+1][N_LENGTH_BINS];
        matrix_inv = new int[MAX_FREQUENCY+1][N_LENGTH_BINS];
        matrix_dup = new int[MAX_FREQUENCY+1][N_LENGTH_BINS];
        matrix_ins = new int[MAX_FREQUENCY+1][N_LENGTH_BINS];
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            length=Integer.parseInt(tokens[1])/LENGTH_QUANTUM;
            if (length>=N_LENGTH_BINS) length=N_LENGTH_BINS-1;
            frequency=Integer.parseInt(tokens[2]);
            if (tokens[0].equalsIgnoreCase("DEL")) matrix_del[frequency][length]++;
            else if (tokens[0].equalsIgnoreCase("INV")) matrix_inv[frequency][length]++;
            else if (tokens[0].equalsIgnoreCase("DUP")) matrix_dup[frequency][length]++;
            else if (tokens[0].equalsIgnoreCase("INS")) matrix_ins[frequency][length]++;
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