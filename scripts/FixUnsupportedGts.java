import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * 
 */
public class FixUnsupportedGts {
    
    /**
     * Remark: the program prints to STDOUT.
     *
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final int MIN_AD = Integer.parseInt(args[1]);
        final int AD_INDEX = Integer.parseInt(args[2]);  // Zero-based
        
        boolean isDiploid;
        char gtPipe;
        int i, p, q;
        int record, nAlts, adRef, adAlt, nSamples;
        long n_cells, n_00_01, n_01_00, n_01_11, n_11_01, n_changes;
        String str;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine(); record=0; nSamples=0; n_cells=0;
        n_00_01=0; n_01_00=0; n_01_11=0; n_11_01=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            nSamples=tokens.length-9;
            // Changing GTs
            for (i=9; i<tokens.length; i++) {
                tokensPrime=tokens[i].split(":");
                nAlts=0; isDiploid=true; adRef=-1; adAlt=-1; gtPipe='x';
                p=tokens[i].indexOf(":");
                if (p==3) {
                    if (tokens[i].charAt(0)=='.' || tokens[i].charAt(2)=='.') continue;
                    n_cells++;
                    isDiploid=true;
                    gtPipe=tokens[i].charAt(1);
                    nAlts=(tokens[i].charAt(0)=='1'?1:0)+(tokens[i].charAt(2)=='1'?1:0);
                    q=tokensPrime[AD_INDEX].indexOf(",");
                    adRef=Integer.parseInt(tokensPrime[AD_INDEX].substring(0,q));
                    adAlt=Integer.parseInt(tokensPrime[AD_INDEX].substring(q+1));
                }
                else if (p==1) {
                    if (tokens[i].charAt(0)=='.') continue;
                    n_cells++;
                    isDiploid=false;
                    gtPipe='x';
                    nAlts=tokens[i].charAt(0)=='1'?1:0;
                    q=tokensPrime[AD_INDEX].indexOf(",");
                    adRef=Integer.parseInt(tokensPrime[AD_INDEX].substring(0,q));
                    adAlt=Integer.parseInt(tokensPrime[AD_INDEX].substring(q+1));
                }
                else {
                    System.err.println("ERROR: unknown GT: "+tokens[i]);
                    System.exit(1);
                }
                if (nAlts==0) {  // 0/0
                    if (adRef<MIN_AD) {
                        tokens[i]=(isDiploid?"0"+gtPipe+"1":"1")+tokens[i].substring(p);
                        n_00_01++;
                    }
                }
                else if (nAlts==1) {  // 0/1
                    if (adAlt<MIN_AD) {
                        tokens[i]=(isDiploid?"0"+gtPipe+"0":"0")+tokens[i].substring(p);
                        n_01_00++;
                    }
                    else if (adRef<MIN_AD) {
                        tokens[i]=(isDiploid?"1"+gtPipe+"1":"1")+tokens[i].substring(p);
                        n_01_11++;
                    }
                }
                else if (nAlts==2) {  // 1/1
                    if (adAlt<MIN_AD) {
                        tokens[i]=(isDiploid?"0"+gtPipe+"1":"0")+tokens[i].substring(p);
                        n_11_01++;
                    }
                }
            }
            record++;
            if (record%1000==0) System.err.println("Processed "+record+" records...");
            // Writing the modified record
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) System.out.print("\t"+tokens[i]);
            System.out.println();
            str=br.readLine();
        }
        br.close();
        System.err.println("Processed "+record+" total records and "+n_cells+" GTs.");
        System.err.println("Number of changes 0/0 -> 0/1: "+n_00_01);
        System.err.println("Number of changes 0/1 -> 0/0: "+n_01_00);
        System.err.println("Number of changes 0/1 -> 1/1: "+n_01_11);
        System.err.println("Number of changes 1/1 -> 0/1: "+n_11_01);
        n_changes=n_00_01+n_01_00+n_01_11+n_11_01;
        System.err.println("Total number of changes: "+n_changes+" ("+(100*((double)n_changes)/n_cells)+"% of GTs)");
    }
    
}