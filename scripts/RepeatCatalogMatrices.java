import java.io.*;


/**
 * 
 */
public class RepeatCatalogMatrices {
    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_BED = args[0];
        final double PURITY_QUANTUM = Double.parseDouble(args[1]);
        
        final int INTERVAL_LENGTH_MAX = 1000;
        final int MOTIF_LENGTH_MAX = 200;
        final int N_REPEATS_MAX = 200;
        
        int i, j;
        int intervalLength, motifLength, nRepeatsInRef;
        double purity;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        int[][] intervalLength_motifLength, intervalLength_nRepeatsInRef, motifLength_nRepeatsInRef;
        
        intervalLength_motifLength = new int[INTERVAL_LENGTH_MAX+1][MOTIF_LENGTH_MAX+1];
        intervalLength_nRepeatsInRef = new int[INTERVAL_LENGTH_MAX+1][N_REPEATS_MAX+1];
        motifLength_nRepeatsInRef = new int[MOTIF_LENGTH_MAX+1][N_REPEATS_MAX+1];
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_BED)));
        str=br.readLine(); str=br.readLine();  // Skipping header
        while (str!=null) { 
            tokens=str.split("\t");
            intervalLength=Integer.parseInt(tokens[2])-Integer.parseInt(tokens[1]);
            motifLength=Integer.parseInt(tokens[3]);
            nRepeatsInRef=Integer.parseInt(tokens[4]);
            purity=Double.parseDouble(tokens[5]);
            if (intervalLength>INTERVAL_LENGTH_MAX) intervalLength=INTERVAL_LENGTH_MAX;
            if (motifLength>MOTIF_LENGTH_MAX) motifLength=MOTIF_LENGTH_MAX;
            if (nRepeatsInRef>N_REPEATS_MAX) nRepeatsInRef=N_REPEATS_MAX;
            intervalLength_motifLength[intervalLength][motifLength]++;
            intervalLength_nRepeatsInRef[intervalLength][nRepeatsInRef]++;
            motifLength_nRepeatsInRef[motifLength][nRepeatsInRef]++;
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        bw = new BufferedWriter(new FileWriter("intervalLength_motifLength.csv"));
        for (i=0; i<INTERVAL_LENGTH_MAX; i++) {
            for (j=0; j<MOTIF_LENGTH_MAX; j++) bw.write(intervalLength_motifLength[i][j]+",");
            bw.newLine();
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter("intervalLength_nRepeatsInRef.csv"));
        for (i=0; i<INTERVAL_LENGTH_MAX; i++) {
            for (j=0; j<N_REPEATS_MAX; j++) bw.write(intervalLength_nRepeatsInRef[i][j]+",");
            bw.newLine();
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter("motifLength_nRepeatsInRef.csv"));
        for (i=0; i<MOTIF_LENGTH_MAX; i++) {
            for (j=0; j<N_REPEATS_MAX; j++) bw.write(motifLength_nRepeatsInRef[i][j]+",");
            bw.newLine();
        }
        bw.close();
    }
    
}