import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Prints to STDOUT a BED representation of file:
 *
 * `repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz`
 *
 * from:
 *
 * https://github.com/broadinstitute/tandem-repeat-catalog/releases/tag/v1.0
 */
public class FilterRepeatCatalog {
    
    private static final String REFERENCE_REGION_STR = "\"ReferenceRegion\": \"";
    private static final int REFERENCE_REGION_LENGTH = REFERENCE_REGION_STR.length();
    
    private static final String CANONICAL_MOTIF_STR = "\"CanonicalMotif\": \"";
    private static final int CANONICAL_MOTIF_LENGTH = CANONICAL_MOTIF_STR.length();
    
    private static final String NUM_REPEATS_IN_REFERENCE_STR = "\"NumRepeatsInReference\": ";
    private static final int NUM_REPEATS_IN_REFERENCE_LENGTH = NUM_REPEATS_IN_REFERENCE_STR.length();
    
    private static final String REFERENCE_REPEAT_PURITY_STR = "\"ReferenceRepeatPurity\": ";
    private static final int REFERENCE_REPEAT_PURITY_LENGTH = REFERENCE_REPEAT_PURITY_STR.length();
    
    private static final String TRS_IN_REGION_STR = "\"TRsInRegion\": ";
    private static final int TRS_IN_REGION_LENGTH = TRS_IN_REGION_STR.length();
    
    private static final String VARIATION_CLUSTER_SIZE_DIFF_STR = "\"VariationClusterSizeDiff\": ";
    private static final int VARIATION_CLUSTER_SIZE_DIFF_LENGTH = VARIATION_CLUSTER_SIZE_DIFF_STR.length();
    
    private static final String LPS_LENGTH_STR = "\"LPSLengthStdevFromHPRC100\": ";
    private static final int LPS_LENGTH_LENGTH = LPS_LENGTH_STR.length();
    
    private static final String LPS_MOTIF_STR = "\"LPSMotifFractionFromHPRC100\": \"";
    private static final int LPS_MOTIF_LENGTH = LPS_MOTIF_STR.length();
    
    private static final String STDEV_ILLUMINA_STR = "\"StdevFromIllumina174k\": ";
    private static final int STDEV_ILLUMINA_LENGTH = STDEV_ILLUMINA_STR.length();
    
    private static final String STDEV_T2T_STR = "\"StdevFromT2TAssemblies\": ";
    private static final int STDEV_T2T_LENGTH = STDEV_T2T_STR.length();

    
    /**
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_JSON_GZ = args[0];
        
        boolean open;
        int p, q, r;
        int strLength, start, end, motifLength, nRepeatsInRef, trsInRegion, clusterSizeDiff;
        double numerator, denominator, purity, lpsLength, lpsMotif, stdevIllumina, stdevT2T;
        String str, chr;
        BufferedReader br;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_JSON_GZ))));
        str=br.readLine(); open=false;
        chr=""; start=-1; end=-1; motifLength=-1; nRepeatsInRef=-1; purity=-1.0; trsInRegion=-1; clusterSizeDiff=-1; lpsLength=-1.0; lpsMotif=-1.0; stdevIllumina=-1.0; stdevT2T=-1.0;
        System.out.println("# chr start end motifLength nRepeatsInRef purity trsInRegion clusterSizeDiff lpsLength lpsMotif stdevIllumina stdevT2T");
        while (str!=null) { 
            if (str.indexOf("{")>=0) { 
                open=true;
                chr=""; start=-1; end=-1; motifLength=-1; nRepeatsInRef=-1; purity=-1.0; trsInRegion=-1; clusterSizeDiff=-1; lpsLength=-1.0; lpsMotif=-1.0; stdevIllumina=-1.0; stdevT2T=-1.0;
                str=br.readLine();
                continue;
            }
            else if (str.indexOf("}")>=0) { 
                System.out.println(chr+"\t"+start+"\t"+end+"\t"+motifLength+"\t"+nRepeatsInRef+"\t"+purity+"\t"+trsInRegion+"\t"+clusterSizeDiff+"\t"+lpsLength+"\t"+lpsMotif+"\t"+stdevIllumina+"\t"+stdevT2T);
                open=false;
                str=br.readLine();
                continue;
            }
            if (!open) { str=br.readLine(); continue; }
            strLength=str.length();
            if (str.charAt(strLength-1)==',') { str=str.substring(0,strLength-1); strLength--; }
            p=str.indexOf(REFERENCE_REGION_STR);
            if (p>=0) {
                q=str.indexOf(":",p+REFERENCE_REGION_LENGTH);
                chr=str.substring(p+REFERENCE_REGION_LENGTH,q);
                r=str.indexOf("-",q+1);
                start=Integer.parseInt(str.substring(q+1,r));
                end=Integer.parseInt(str.substring(r+1,strLength-1));
                str=br.readLine();
                continue;
            }
            p=str.indexOf(CANONICAL_MOTIF_STR);
            if (p>=0) {
                motifLength=strLength-p-CANONICAL_MOTIF_LENGTH-1;
                str=br.readLine();
                continue;
            }
            p=str.indexOf(NUM_REPEATS_IN_REFERENCE_STR);
            if (p>=0) {
                nRepeatsInRef=Integer.parseInt(str.substring(p+NUM_REPEATS_IN_REFERENCE_LENGTH,strLength));
                str=br.readLine();
                continue;
            }
            p=str.indexOf(REFERENCE_REPEAT_PURITY_STR);
            if (p>=0) {
                purity=Double.parseDouble(str.substring(p+REFERENCE_REPEAT_PURITY_LENGTH,strLength));
                str=br.readLine();
                continue;
            }
            p=str.indexOf(TRS_IN_REGION_STR);
            if (p>=0) {
                trsInRegion=Integer.parseInt(str.substring(p+TRS_IN_REGION_LENGTH,strLength));
                str=br.readLine();
                continue;
            }
            p=str.indexOf(VARIATION_CLUSTER_SIZE_DIFF_STR);
            if (p>=0) {
                clusterSizeDiff=Integer.parseInt(str.substring(p+VARIATION_CLUSTER_SIZE_DIFF_LENGTH,strLength));
                str=br.readLine();
                continue;
            }
            p=str.indexOf(LPS_LENGTH_STR);
            if (p>=0) {
                lpsLength=Double.parseDouble(str.substring(p+LPS_LENGTH_LENGTH,strLength));
                str=br.readLine();
                continue;
            }
            p=str.indexOf(LPS_MOTIF_STR);
            if (p>=0) {
                q=str.indexOf(":",p+LPS_MOTIF_LENGTH);
                r=str.indexOf("/",q+1);
                numerator=Double.parseDouble(str.substring(q+2,r));
                denominator=Double.parseDouble(str.substring(r+1,strLength-1));
                lpsMotif=numerator/denominator;
                str=br.readLine();
                continue;
            }
            p=str.indexOf(STDEV_ILLUMINA_STR);
            if (p>=0) {
                stdevIllumina=Double.parseDouble(str.substring(p+STDEV_ILLUMINA_LENGTH,strLength));
                str=br.readLine();
                continue;
            }
            p=str.indexOf(STDEV_T2T_STR);
            if (p>=0) {
                stdevT2T=Double.parseDouble(str.substring(p+STDEV_T2T_LENGTH,strLength));
                str=br.readLine();
                continue;
            }
            str=br.readLine();
        }
        br.close();
    }
    
}