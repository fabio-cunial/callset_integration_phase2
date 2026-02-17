import java.io.*;


/**
 * Creates a (sample,posBin) matrix for chrX and chrY that stores how many
 * dipcall BED intervals belong to each bin in each sample.
 */
public class PlotDipcallBeds {

	public static void main(String[] args) throws IOException {
	    final int CHRX_LENGTH = Integer.parseInt(args[0]);
        final int CHRY_LENGTH = Integer.parseInt(args[1]);
        final int QUANTUM = Integer.parseInt(args[2]);
        final String FILE_LIST = args[3];
        final int N_FILES = Integer.parseInt(args[4]);
        
        final int N_CHUNKS_X = (CHRX_LENGTH+QUANTUM-1)/QUANTUM;
        final int N_CHUNKS_Y = (CHRY_LENGTH+QUANTUM-1)/QUANTUM;
        
		int i, j, p, q;
        int first, last;
        String str, strPrime, chrom;
        BufferedReader br, brPrime;
        BufferedWriter bw;
		int[] histogramX, histogramY;
        int[][] matrixX, matrixY;
        
        // Allocating memory
        histogramX = new int[N_CHUNKS_X];
        histogramY = new int[N_CHUNKS_Y];
        matrixX = new int[N_FILES][N_CHUNKS_X];
        matrixY = new int[N_FILES][N_CHUNKS_Y];
        
        // Loading matrices
        br = new BufferedReader(new FileReader(FILE_LIST));
        str=br.readLine(); i=0;
        while (str!=null) {
            brPrime = new BufferedReader(new FileReader(str));
            strPrime=brPrime.readLine();
            while (strPrime!=null) {
                p=strPrime.indexOf('\t'); q=strPrime.indexOf('\t',p+1);
                chrom=strPrime.substring(0,p);
                first=Integer.parseInt(strPrime.substring(p+1,q));
                last=Integer.parseInt(strPrime.substring(q+1))-1;
                if (chrom.equalsIgnoreCase("chrX")) {
                    for (j=first/QUANTUM; j<=last/QUANTUM; j++) matrixX[i][j]++;
                }
                else if (chrom.equalsIgnoreCase("chrY")) {
                    for (j=first/QUANTUM; j<=last/QUANTUM; j++) matrixY[i][j]++;
                }
                strPrime=brPrime.readLine();
            }
            brPrime.close();
            str=br.readLine(); i++;
        }
        br.close();
		for (i=0; i<N_FILES; i++) {
		    for (j=0; j<N_CHUNKS_X; j++) histogramX[j]+=matrixX[i][j];
		}
		for (i=0; i<N_FILES; i++) {
		    for (j=0; j<N_CHUNKS_Y; j++) histogramY[j]+=matrixY[i][j];
		}
        
        // Outputting
		bw = new BufferedWriter(new FileWriter("histogram_x.txt"));
        for (i=0; i<N_CHUNKS_X; i++) bw.write(histogramX[i]+"\n");
        bw.close();
		bw = new BufferedWriter(new FileWriter("histogram_y.txt"));
        for (i=0; i<N_CHUNKS_Y; i++) bw.write(histogramY[i]+"\n");
        bw.close();
		bw = new BufferedWriter(new FileWriter("matrix_x.txt"));
        for (i=0; i<N_FILES; i++) {
            for (j=0; j<N_CHUNKS_X; j++) bw.write(matrixX[i][j]+",");
            bw.newLine();
        }
        bw.close();
		bw = new BufferedWriter(new FileWriter("matrix_y.txt"));
        for (i=0; i<N_FILES; i++) {
            for (j=0; j<N_CHUNKS_Y; j++) bw.write(matrixY[i][j]+",");
            bw.newLine();
        }
        bw.close();
	}

}