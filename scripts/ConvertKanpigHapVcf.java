import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.*;


/**
 * Transforms a single-sample, re-genotyped haps VCF, into the corresponding 
 * records VCF. Every record is printed in output, including those that are not
 * genotyped as present.
 */
public class ConvertKanpigHapVcf {
    
    /**
     * Remark: the program prints the output VCF to STDOUT.
     *
     * @param args 
     * 0: the output of kanpig when given in input the cohort haps VCF and one 
     *    BAM; assumed to contain only one sample;
     * 1: every row corresponds to a chunk and has format: 
     *    `records_file_path \t map_file_path`. The rows of every 
     *    `map_file_path`, in order, are assumed to correspond to the rows of 
     *    the haps VCF.
     */
    public static void main(String[] args) throws IOException {
        final String HAPS_VCF_GZ = args[0];
        final String CHUNKS_TSV = args[1];

        final int CAPACITY = 10;  // Arbitrary
        final int QUANTUM = 100;  // Arbitrary
        
        boolean found;
        int i, j, k, p;
        int nRecords, nHaps, nChunks;
        String str, strHaps, strChunks, recordsPath, mapPath;
        BufferedReader br, brHaps, brChunks;
        int[] mapLast;
        String[] tokens, records;
        boolean[][] recordGTs;
        int[][] map;
        
        map = new int[CAPACITY][0];
        mapLast = new int[CAPACITY];
        records = new String[CAPACITY];
        recordGTs = new boolean[CAPACITY][2];
        brHaps = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(HAPS_VCF_GZ))));
        strHaps=brHaps.readLine();
        while (strHaps!=null) {
            if (strHaps.charAt(0)!='#') break;
            System.out.println(strHaps);
            strHaps=brHaps.readLine();
        }
        brChunks = new BufferedReader(new FileReader(CHUNKS_TSV));
        strChunks=brChunks.readLine(); nChunks=0;
        while (strChunks!=null) {  // For each chunk
            p=strChunks.indexOf("\t");
            recordsPath=strChunks.substring(0,p);
            mapPath=strChunks.substring(p+1);
            
            // Loading the map file of the chunk
            br = new BufferedReader(new FileReader(mapPath));
            str=br.readLine(); i=-1;  
            while (str!=null) {
                tokens=str.split(",");
                i++;
                if (i==map.length) {
                    int[][] newMap = new int[map.length*2][0];
                    System.arraycopy(map,0,newMap,0,map.length);
                    map=newMap;
                    int[] newMapLast = new int[mapLast.length*2];
                    System.arraycopy(mapLast,0,newMapLast,0,mapLast.length);
                    mapLast=newMapLast;
                }
                if (map[i]==null || map[i].length<tokens.length) map[i] = new int[tokens.length];
                for (j=0; j<tokens.length; j++) map[i][j]=Integer.parseInt(tokens[j]);
                mapLast[i]=tokens.length-1;
                str=br.readLine();
            }
            br.close();
            nHaps=i+1;
            
            // Loading the records file of the chunk
            br = new BufferedReader(new FileReader(recordsPath));
            str=br.readLine(); i=-1;
            while (str!=null) {
                i++;
                if (i==records.length) {
                    String[] newRecords = new String[records.length*2];
                    System.arraycopy(records,0,newRecords,0,records.length);
                    records=newRecords;
                }
                records[i]=str;
                str=br.readLine();
            }
            br.close();
            nRecords=i+1;
            
            // Converting hap GTs to record GTs
            if (recordGTs.length<nRecords) recordGTs = new boolean[nRecords][2];
            for (i=0; i<nRecords; i++) { recordGTs[i][0]=false; recordGTs[i][1]=false; }
            i=0;
            while (i<nHaps) {
                tokens=strHaps.split("\t");
                if (tokens[9].charAt(0)=='1') {
                    for (k=0; k<=mapLast[i]; k++) recordGTs[map[i][k]][0]=true;
                }
                if (tokens[9].charAt(2)=='1') {
                    for (k=0; k<=mapLast[i]; k++) recordGTs[map[i][k]][1]=true;
                }
                strHaps=brHaps.readLine();
                i++;
            }
            
            // Outputting all records
            for (i=0; i<nRecords; i++) {
                System.out.print(records[i]);
                System.out.print("\tGT\t");
                System.out.print(recordGTs[i][0]?"1|":"0|");
                System.out.println(recordGTs[i][1]?"1":"0");
            }
            
            // Next iteration
            nChunks++;
            if (nChunks%QUANTUM==0) System.err.println("Processed "+nChunks+" chunks...");
            strChunks=brChunks.readLine();
        }
        brChunks.close();
        System.err.println("Processed all "+nChunks+" chunks");
    }
    
}