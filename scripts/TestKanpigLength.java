import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class TestKanpigLength {
    
    /**
     * @args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final String OUTPUT_PREFIX = args[1];
        
        final int[] LENGTH_SEPARATORS = new int[] {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
        final int N_LENGTH_SEPARATORS = LENGTH_SEPARATORS.length;
        final int MAX_LENGTH_SEPARATOR = LENGTH_SEPARATORS[N_LENGTH_SEPARATORS-1];
        
        char a, b;
        int i, j, k;
        int length, from, to, bin;
        String str, svtype;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        int[] recordsAtLength_all, recordsAtLength_del, recordsAtLength_ins;
        int[][][] matrix_all, matrix_del, matrix_ins;
        
        // Building matrix
        recordsAtLength_all = new int[N_LENGTH_SEPARATORS];
        recordsAtLength_del = new int[N_LENGTH_SEPARATORS];
        recordsAtLength_ins = new int[N_LENGTH_SEPARATORS];
        matrix_all = new int[3][3][N_LENGTH_SEPARATORS];
        matrix_del = new int[3][3][N_LENGTH_SEPARATORS];
        matrix_ins = new int[3][3][N_LENGTH_SEPARATORS];
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            svtype=tokens[0];
            length=Integer.parseInt(tokens[1]);
            if (length>MAX_LENGTH_SEPARATOR) {
                str=br.readLine();
                continue;
            }
            a=tokens[2].charAt(0); b=tokens[2].charAt(2); 
            from=(a=='1'?1:0)+(b=='1'?1:0);
            a=tokens[3].charAt(0); b=tokens[3].charAt(2); 
            to=(a=='1'?1:0)+(b=='1'?1:0);
            bin=Arrays.binarySearch(LENGTH_SEPARATORS,length);
            if (bin<0) bin=-bin-1;
            recordsAtLength_all[bin]++;
            if (svtype.equalsIgnoreCase("DEL")) recordsAtLength_del[bin]++;
            else if (svtype.equalsIgnoreCase("INS")) recordsAtLength_ins[bin]++;
            matrix_all[from][to][bin]++;
            if (svtype.equalsIgnoreCase("DEL")) matrix_del[from][to][bin]++;
            else if (svtype.equalsIgnoreCase("INS")) matrix_ins[from][to][bin]++;
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) {
                bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_all_"+i+"_"+j+".csv"));
                for (k=0; k<N_LENGTH_SEPARATORS; k++) bw.write((((double)matrix_all[i][j][k])/recordsAtLength_all[k])+"\n");
                bw.close();
            }
        }
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) {
                bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_del_"+i+"_"+j+".csv"));
                for (k=0; k<N_LENGTH_SEPARATORS; k++) bw.write((((double)matrix_del[i][j][k])/recordsAtLength_del[k])+"\n");
                bw.close();
            }
        }
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) {
                bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"_ins_"+i+"_"+j+".csv"));
                for (k=0; k<N_LENGTH_SEPARATORS; k++) bw.write((((double)matrix_ins[i][j][k])/recordsAtLength_ins[k])+"\n");
                bw.close();
            }
        }
    }

}