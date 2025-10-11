import java.util.*;
import java.io.*;


/**
 * 
 */
public class PrintRepeatCatalogPairsBash {
    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_BED = args[0];
        final int N_THREADS = Integer.parseInt(args[1]);
        
        final int[] MOTIF_LENGTH = new int[] {0,1,2,4,8,16,32,64,128};
        final int[] N_REPEATS = new int[] {0,1,2,4,8,16,32,64};
        final int[] INTERVAL_LENGTH = {0,2,4,6,8,10,12,14,16,18,20,40,80,160,320};
        
        int i, j;
        int quantum;
        BufferedWriter bw;
        Vector<String> instructions = new Vector<String>();
        
        // Collecting instructions
        for (i=1; i<MOTIF_LENGTH.length; i++) {
            for (j=1; j<N_REPEATS.length; j++) instructions.add("awk 'BEGIN{}{ if ($4>"+MOTIF_LENGTH[i-1]+" && $4<="+MOTIF_LENGTH[i]+" && $5>"+N_REPEATS[j-1]+" && $5<="+N_REPEATS[j]+") print $0 }' "+INPUT_BED+" > ml_"+MOTIF_LENGTH[i]+"_nr_"+N_REPEATS[j]+".bed \n");
        }
        for (i=1; i<MOTIF_LENGTH.length; i++) {
            for (j=1; j<INTERVAL_LENGTH.length; j++) instructions.add("awk 'BEGIN{}{ if ($4>"+MOTIF_LENGTH[i-1]+" && $4<="+MOTIF_LENGTH[i]+" && ($3-$2)>"+INTERVAL_LENGTH[j-1]+" && ($3-$2)<="+INTERVAL_LENGTH[j]+") print $0 }' "+INPUT_BED+" > ml_"+MOTIF_LENGTH[i]+"_il_"+INTERVAL_LENGTH[j]+".bed \n");
        }
        for (i=1; i<N_REPEATS.length; i++) {
            for (j=1; j<INTERVAL_LENGTH.length; j++) instructions.add("awk 'BEGIN{}{ if ($5>"+N_REPEATS[i-1]+" && $5<="+N_REPEATS[i]+" && ($3-$2)>"+INTERVAL_LENGTH[j-1]+" && ($3-$2)<="+INTERVAL_LENGTH[j]+") print $0 }' "+INPUT_BED+" > nr_"+N_REPEATS[i]+"_il_"+INTERVAL_LENGTH[j]+".bed \n");
        }
        
        // Printing bash scripts
        quantum=(int)Math.ceil(((double)instructions.size())/N_THREADS);
        for (i=0; i<N_THREADS; i++) {
            bw = new BufferedWriter(new FileWriter("thread_"+i+".sh"));
            bw.write("#!/bin/bash \n");
            bw.write("# \n");
            for (j=i*quantum; j<(i+1)*quantum && j<instructions.size(); j++) bw.write(instructions.elementAt(j));
            bw.close();
        }
    }
    
}