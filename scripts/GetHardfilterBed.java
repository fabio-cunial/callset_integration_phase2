import java.util.*;
import java.io.*;


/**
 * Creates BED files that limit chromosomes to their parts outside telomeres 
 * and centromeres.
 */
public class GetHardfilterBed {
    
    /**
     * @param 
     * 1: centromeres BED downloaded from the UCSC Genome Browser (Mapping and 
     *    Sequencing > Centromeres > BED);
     * 2: assembly gaps TSV downloaded from the UCSC Genome Browser (Mapping and
     *    Sequencing > Gap > All fields from selected table > tsv).
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FAI = args[0];
        final String INPUT_CENTROMERES_BED = args[1];
        final String INPUT_GAPS_TSV = args[2];
        
        final int BUFFER_BPS = Integer.parseInt(args[3]);
        final String BUFFER_BPS_BED = args[4];
        final double BUFFER_FRACTION = Double.parseDouble(args[5]);
        final String BUFFER_FRACTION_BED = args[6];
        
        final int MIN_OUTPUT_LENGTH = 1000;
        
        int pos, currentStart, currentEnd, currentLength, cStart, cEnd, tStart, tEnd, buffer, previousStart, previousEnd;
        String str, currentChr;
        BufferedReader br;
        BufferedWriter bw1, bw2;
        HashMap<String,Integer> chrLengths, centromereStart, centromereEnd, telomereStart, telomereEnd;
        String[] tokens;
        
        // Loading chromosome lengths
        chrLengths = new HashMap();
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_FAI)));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            chrLengths.put(tokens[0],Integer.valueOf(Integer.parseInt(tokens[1])));
            str=br.readLine();
        }
        br.close();
        
        // Loading centromere intervals
        centromereStart = new HashMap(); centromereEnd = new HashMap();
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_CENTROMERES_BED)));
        str=br.readLine(); str=br.readLine();  // Skipping header
        currentChr=""; currentStart=Integer.MAX_VALUE; currentEnd=-1;
        while (str!=null) {
            tokens=str.split("\t");
            if (!tokens[0].equalsIgnoreCase(currentChr)) {
                if (currentChr.length()>0) { centromereStart.put(currentChr,Integer.valueOf(currentStart)); centromereEnd.put(currentChr,Integer.valueOf(currentEnd)); }
                currentChr=tokens[0]; currentStart=Integer.parseInt(tokens[1]); currentEnd=Integer.parseInt(tokens[2]);
            }
            else {
                pos=Integer.parseInt(tokens[1]);
                if (pos<currentStart) currentStart=pos;
                pos=Integer.parseInt(tokens[2]);
                if (pos>currentEnd) currentEnd=pos;
            }
            str=br.readLine();
        }
        if (currentChr.length()>0) { centromereStart.put(currentChr,Integer.valueOf(currentStart)); centromereEnd.put(currentChr,Integer.valueOf(currentEnd)); }
        br.close();
        
        // Loading telomere intervals, and updating centromere intervals.
        telomereStart = new HashMap(); telomereEnd = new HashMap();
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_GAPS_TSV)));
        str=br.readLine(); str=br.readLine();  // Skipping header
        while (str!=null) {
            tokens=str.split("\t");
            currentChr=tokens[1];
            if (!chrLengths.containsKey(currentChr)) { str=br.readLine(); continue; }
            currentLength=chrLengths.get(currentChr).intValue();
            currentStart=Integer.parseInt(tokens[2]);
            currentEnd=Integer.parseInt(tokens[3]);
            if (tokens[7].equalsIgnoreCase("telomere")) {
                if (currentStart==0) {
                    previousEnd=telomereEnd.containsKey(currentChr)?telomereEnd.get(currentChr).intValue():-1;
                    if (currentEnd>previousEnd) telomereEnd.put(currentChr,Integer.valueOf(currentEnd));
                }
                else if (currentEnd==currentLength) {
                    previousStart=telomereStart.containsKey(currentChr)?telomereStart.get(currentChr).intValue():Integer.MAX_VALUE;
                    if (currentStart<previousStart) telomereStart.put(currentChr,Integer.valueOf(currentStart));
                }
                else {
                    System.err.println("ERROR: telomere interval not aligned to the beginning/end of the chromosome:");
                    System.err.println(str);
                    System.exit(1);
                }
            }
            else if (tokens[7].equalsIgnoreCase("heterochromatin")) {
                if (currentChr.equalsIgnoreCase("chrY")) {
                    previousStart=telomereStart.containsKey(currentChr)?telomereStart.get(currentChr).intValue():Integer.MAX_VALUE;
                    if (currentStart<previousStart) telomereStart.put(currentChr,Integer.valueOf(currentStart));
                }
                else {
                    previousStart=centromereStart.containsKey(currentChr)?centromereStart.get(currentChr).intValue():Integer.MAX_VALUE;
                    if (currentStart<previousStart) centromereStart.put(currentChr,Integer.valueOf(currentStart));
                    previousEnd=centromereEnd.containsKey(currentChr)?centromereEnd.get(currentChr).intValue():-1;
                    if (currentEnd>previousEnd) centromereEnd.put(currentChr,Integer.valueOf(currentEnd));
                }
            }
            else if (tokens[7].equalsIgnoreCase("short_arm")) {
                // chr13,14,15,21,22
                previousStart=centromereStart.containsKey(currentChr)?centromereStart.get(currentChr).intValue():Integer.MAX_VALUE;
                if (currentStart<previousStart) centromereStart.put(currentChr,Integer.valueOf(currentStart));
            }
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        bw1 = new BufferedWriter(new FileWriter(BUFFER_BPS_BED));
        bw2 = new BufferedWriter(new FileWriter(BUFFER_FRACTION_BED));
        br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_FAI)));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            currentChr=tokens[0];
            currentLength=chrLengths.get(currentChr).intValue();
            if (!telomereStart.containsKey(currentChr) || !telomereEnd.containsKey(currentChr) || !centromereStart.containsKey(currentChr) || !centromereEnd.containsKey(currentChr)) {
                System.err.println("WARNING: "+currentChr+" does not have full telomere/centromere information.");
                bw1.write(currentChr+"\t0\t"+currentLength+"\n");
                bw2.write(currentChr+"\t0\t"+currentLength+"\n");
                str=br.readLine();
                continue;
            }
            tStart=telomereStart.get(currentChr).intValue();
            tEnd=telomereEnd.get(currentChr).intValue();
            cStart=centromereStart.get(currentChr).intValue();
            cEnd=centromereEnd.get(currentChr).intValue();
            buffer=BUFFER_BPS;
            if (currentLength<4*buffer) bw1.write(currentChr+"\t0\t"+currentLength+"\n");
            else {
                if ((cStart-buffer)-(tEnd+buffer)>=MIN_OUTPUT_LENGTH) bw1.write(currentChr+"\t"+(tEnd+buffer)+"\t"+(cStart-buffer)+"\n");
                if ((tStart-buffer)-(cEnd+buffer)>=MIN_OUTPUT_LENGTH) bw1.write(currentChr+"\t"+(cEnd+buffer)+"\t"+(tStart-buffer)+"\n");
            }
            buffer=(int)Math.ceil(currentLength*BUFFER_FRACTION);
            if (currentLength<4*buffer) bw2.write(currentChr+"\t0\t"+currentLength+"\n");
            else {
                if ((cStart-buffer)-(tEnd+buffer)>=MIN_OUTPUT_LENGTH) bw2.write(currentChr+"\t"+(tEnd+buffer)+"\t"+(cStart-buffer)+"\n");
                if ((tStart-buffer)-(cEnd+buffer)>=MIN_OUTPUT_LENGTH) bw2.write(currentChr+"\t"+(cEnd+buffer)+"\t"+(tStart-buffer)+"\n");
            }
            str=br.readLine();
        }
        br.close(); bw1.close(); bw2.close();
    }
    
}
